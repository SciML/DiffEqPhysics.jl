function gather_bodies_initial_coordinates(system::NBodySystem)
    bodies = system.bodies;
    n = length(bodies)
    u0 = zeros(3, n)
    v0 = zeros(3, n)

    for i = 1:n
        u0[:, i] = bodies[i].r 
        v0[:, i] = bodies[i].v
    end 

    (u0, v0, n)
end

function pairwise_lennard_jones_acceleration!(dv,
    rs,
    i::Int,
    n::Int,
    bodies::Vector{<:MassBody},
    p::LennardJonesParameters,
    pbc::BoundaryConditions)

    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    
    for j = 1:n
        if j != i
            rij = @MVector [ri[1] - rs[1, j], ri[2] - rs[2, j], ri[3] - rs[3, j]]

            apply_boundary_conditions!(rij, pbc)
            
            rij_2 = dot(rij, rij)
            σ_rij_6 = (p.σ2 / rij_2)^3
            σ_rij_12 = σ_rij_6^2

            if rij_2 < p.R2
                accel += (2 * σ_rij_12 - σ_rij_6 ) * rij / rij_2
            end            
        end
    end   

    dv .=  24 * p.ϵ * accel
end

function apply_boundary_conditions!(rij, pbc::PeriodicBoundaryConditions)
    for x in (1, 2, 3) 
        while (rij[x] > pbc[2x]) rij[x] -= (pbc[2x] - pbc[2x - 1]) end
        while (rij[x] < pbc[2x - 1]) rij[x] += (pbc[2x] - pbc[2x - 1]) end
    end
end

function apply_boundary_conditions!(rij, pbc::BoundaryConditions)
end

function pairwise_electrostatic_acceleration!(dv,
    rs,
    i::Int,
    n::Int,
    bodies::Vector{<:ChargedParticle},
    p::ElectrostaticParameters)

    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j = 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            accel += p.k * bodies[i].q * bodies[j].q / bodies[i].m * (ri - rj) / norm(ri - rj)^3
        end
    end    

    dv .= accel
end


function gravitational_acceleration!(dv, 
    rs,
    i::Int,
    n::Int,
    bodies::Vector{<:MassBody},
    p::GravitationalParameters)
    
    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j = 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            accel -= p.G * bodies[j].m * rij / norm(rij)^3
        end
    end
    
    dv .= accel
end

function gather_accelerations_for_potentials(simulation::NBodySimulation{<:PotentialNBodySystem})
    acceleration_functions = []
    
    for (potential, parameters) in simulation.system.potentials
        push!(acceleration_functions, get_accelerating_function(parameters, simulation))
    end

    acceleration_functions
end

function get_accelerating_function(parameters::LennardJonesParameters, simulation::NBodySimulation)
    (dv, u, v, t, i) -> pairwise_lennard_jones_acceleration!(dv, u, i, length(simulation.system.bodies), simulation.system.bodies, parameters, simulation.boundary_conditions)
end

function get_accelerating_function(parameters::ElectrostaticParameters, simulation::NBodySimulation)
    (dv, u, v, t, i) -> pairwise_electrostatic_acceleration!(dv, u, i, length(simulation.system.bodies), simulation.system.bodies, parameters)
end

function get_accelerating_function(parameters::GravitationalParameters, simulation::NBodySimulation)
    (dv, u, v, t, i) -> gravitational_acceleration!(dv, u, i, length(simulation.system.bodies), simulation.system.bodies, parameters)
end

function gather_accelerations_for_potentials(simulation::NBodySimulation{CustomAccelerationSystem})
    acceleration_functions = []
    push!(acceleration_functions, simulation.system.acceleration)
    acceleration_functions
end

function DiffEqBase.ODEProblem(simulation::NBodySimulation{<:PotentialNBodySystem})
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation.system)
    
    acceleration_functions = gather_accelerations_for_potentials(simulation)   

    function ode_system!(du, u, p, t)
        du[:, 1:n] = @view u[:, n + 1:2n];

        @inbounds for i = 1:n
            for acceleration! in acceleration_functions
                a = @view du[:, n + i]
                acceleration!(a, u[:, 1:n], u[:, n + 1:end], t, i);                
            end
        end 
    end

    return ODEProblem(ode_system!, hcat(u0, v0), simulation.tspan)
end

function DiffEqBase.SecondOrderODEProblem(simulation::NBodySimulation{<:PotentialNBodySystem})
    
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation.system)

    acceleration_functions = gather_accelerations_for_potentials(simulation)  

    function soode_system!(dv, v, u, p, t)
        @inbounds for i = 1:n
            for acceleration! in acceleration_functions
                a = @view dv[:, i] 
                acceleration!(a, u, v, t, i);       
            end
        end 
    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan)
end
