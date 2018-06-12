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
    i::Integer,
    n::Integer,
    ms::Vector{<:Real},
    p::LennardJonesParameters,
    pbc::BoundaryConditions)

    force = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    
    for j = 1:n

        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]

            rij = apply_boundary_conditions!(ri, rj, pbc, p)
            
            if rij[1] < Inf
                rij_2 = dot(rij, rij)
                σ_rij_6 = (p.σ2 / rij_2)^3
                σ_rij_12 = σ_rij_6^2
                force += (2 * σ_rij_12 - σ_rij_6 ) * rij / rij_2    
            end        
        end
    end   
    dv .+=  24 * p.ϵ * force / ms[i]
end

function apply_boundary_conditions!(ri, rj, pbc::PeriodicBoundaryConditions, p::LennardJonesParameters)
    res = @MVector [Inf, Inf, Inf]
    for dx in [0,-pbc[2],pbc[2]], dy in [0,-pbc[4],pbc[4]], dz in [0,-pbc[6],pbc[6]]
        rij = @MVector [ri[1] - rj[1] + dx, ri[2] - rj[2] + dy, ri[3] - rj[3] + dz]
        for x in (1, 2, 3) 
            rij[x] -= (pbc[2x] - pbc[2x - 1]) * div(rij[x] , (pbc[2x] - pbc[2x - 1]))
        end
        if  dot(rij, rij) < p.R2
            res = rij
            break
        end
    end
    return res
end

function apply_boundary_conditions!(ri, rj, pbc::CubicPeriodicBoundaryConditions, p::LennardJonesParameters)
    res = @MVector [Inf, Inf, Inf]
    shifts = @SVector [0,-pbc.L,pbc.L]
    for dx in shifts, dy in shifts, dz in shifts
        rij = @MVector [ri[1] - rj[1] + dx, ri[2] - rj[2] + dy, ri[3] - rj[3] + dz]
        for x in (1, 2, 3) 
            rij[x] -= pbc.L * div(rij[x] , pbc.L)
        end
        if  dot(rij, rij) < p.R2
            res = rij
            break
        end
    end
    return res
end

function apply_boundary_conditions!(ri, rj, pbc::BoundaryConditions, p::PotentialParameters)
    return ri - rj
end

function pairwise_electrostatic_acceleration!(dv,
    rs,
    i::Integer,
    n::Integer,
    qs::Vector{<:Real},
    ms::Vector{<:Real},
    p::ElectrostaticParameters)

    force = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j = 1:n
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            force += p.k * qs[i] * qs[j] * (ri - rj) / norm(ri - rj)^3
        end
    end    

    dv .+= force / ms[i]
end


function gravitational_acceleration!(dv, 
    rs,
    i::Integer,
    n::Integer,
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
    
    dv .+= accel
end

function magnetostatic_dipdip_acceleration!(dv, 
    rs,
    i::Integer,
    n::Integer,
    bodies::Vector{<:MagneticParticle},
    p::MagnetostaticParameters)

    force = @SVector [0.0, 0.0, 0.0];
    mi = bodies[i].mm
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j = 1:n
        if j != i
            mj = bodies[j].mm
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            rij4 = dot(rij, rij)^2
            r =  rij / norm(rij)
            mir = dot(mi, r)
            mij = dot(mj, r)
            force += (mi * mij + mj * mir + r * dot(mi, mj) - 5 * r * mir * mij) / rij4
        end
    end
    
    dv .+= 3 * p.μ_4π * force / bodies[i].m
end

function elastic_potential_acceleration(dv, 
    rs,
    i::Integer,
    n::Integer,
    ms::Vector{<:Real},
    neighbours::Vector{Vector{Tuple{Int,Float64}}},
    p::SPCFwParameters)
    
    force = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for (j,k) in neighbours[i]
        rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
        rij = ri - rj
        r = norm(rij)
        d = r-p.rOH
        force -= d*p.k*rij / r
    end
    
    dv .+= force / ms[i]
end

function gather_accelerations_for_potentials(simulation::NBodySimulation{<:PotentialNBodySystem})
    acceleration_functions = []
    
    for (potential, parameters) in simulation.system.potentials
        push!(acceleration_functions, get_accelerating_function(parameters, simulation))
    end

    acceleration_functions
end

function get_accelerating_function(parameters::LennardJonesParameters, simulation::NBodySimulation)
    ms = obtain_data_for_lennard_jones_interaction(simulation.system)
    (dv, u, v, t, i) -> pairwise_lennard_jones_acceleration!(dv, u, i, length(simulation.system.bodies), ms, parameters, simulation.boundary_conditions)
end

function get_accelerating_function(parameters::ElectrostaticParameters, simulation::NBodySimulation)
    (qs, ms) = obtain_data_for_electrostatic_interaction(simulation.system)
    (dv, u, v, t, i) -> pairwise_electrostatic_acceleration!(dv, u, i, length(simulation.system.bodies), qs, ms, parameters)
end

function obtain_data_for_lennard_jones_interaction(system::PotentialNBodySystem)
    bodies = system.bodies
    n = length(bodies)
    ms = zeros(Real ,n)
    for i = 1:n
        ms[i] = bodies[i].m
    end 
    return ms
end

function obtain_data_for_lennard_jones_interaction(system::WaterSPCFw)
    bodies = system.bodies
    n = length(bodies)
    ms = zeros(Real, 3*n)
    for i = 1:n
        ms[3*(i-1)+1] = system.mO
        ms[3*(i-1)+2] = system.mH
        ms[3*(i-1)+3] = system.mH
    end 
    return ms
end

function obtain_data_for_electrostatic_interaction(system::PotentialNBodySystem)
    bodies = system.bodies
    n = length(bodies)
    qs = zeros(Real, n)
    ms = zeros(Real ,n)
    for i = 1:n
        qs[i] = bodies[i].q
        ms[i] = bodies[i].m
    end 
    return (qs, ms)
end

function obtain_data_for_electrostatic_interaction(system::WaterSPCFw)
    bodies = system.bodies
    n = length(bodies)
    qs = zeros(Real, 3*n)
    ms = zeros(Real, 3*n)
    for i = 1:n
        qs[3*(i-1)+1] = system.qO
        qs[3*(i-1)+2] = system.qH
        qs[3*(i-1)+3] = system.qH
        ms[3*(i-1)+1] = system.mO
        ms[3*(i-1)+2] = system.mH
        ms[3*(i-1)+3] = system.mH
    end 
    return (qs, ms)
end

function obtain_data_for_electrostatic_interaction(system::NBodySystem)
    bodies = system.bodies
    n = length(bodies)
    qs = zeros(Real,n)
    ms = zeros(Real,n)
    return (qs, ms)
end

function get_accelerating_function(parameters::GravitationalParameters, simulation::NBodySimulation)
    (dv, u, v, t, i) -> gravitational_acceleration!(dv, u, i, length(simulation.system.bodies), simulation.system.bodies, parameters)
end

function get_accelerating_function(parameters::MagnetostaticParameters, simulation::NBodySimulation)
    (dv, u, v, t, i) -> magnetostatic_dipdip_acceleration!(dv, u, i, length(simulation.system.bodies), simulation.system.bodies, parameters)
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
            a = MVector(0.0,0.0,0.0)
            for acceleration! in acceleration_functions
                acceleration!(a, u[:, 1:n], u[:, n + 1:end], t, i);                
            end
            du[:, n + i] .= a
        end 
    end

    return ODEProblem(ode_system!, hcat(u0, v0), simulation.tspan)
end

function DiffEqBase.SecondOrderODEProblem(simulation::NBodySimulation{<:PotentialNBodySystem})
    
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation.system)

    acceleration_functions = gather_accelerations_for_potentials(simulation)  

    function soode_system!(dv, v, u, p, t)
        @inbounds for i = 1:n
            a = MVector(0.0,0.0,0.0)
            for acceleration! in acceleration_functions
                acceleration!(a, u, v, t, i);    
            end
            dv[:, i].=a
        end 
    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan)
end

function DiffEqBase.SecondOrderODEProblem(simulation::NBodySimulation{<:WaterSPCFw})
    
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation.system)

    (o_acelerations, h_acelerations) = gather_accelerations_for_potentials(simulation)   

    function soode_system!(dv, v, u, p, t)
        @inbounds for i = 1:n
            for acceleration! in o_acelerations
                a = @view dv[:, 3*(i-1)] 
                acceleration!(a, u, v, t, i);    
            end
        end
        @inbounds for i in 1:n, j in (2,3)
            for acceleration! in h_acelerations
                a = @view dv[:, 3*(i-1)+j] 
                acceleration!(a, u, v, t, 3*(i-1)+j);    
            end
        end 
    end

    SecondOrderODEProblem(soode_system!, v0, u0, simulation.tspan)
end

function gather_bodies_initial_coordinates(system::WaterSPCFw)
    bodies = system.bodies;
    n = length(bodies)
    u0 = zeros(3, 3*n)
    v0 = zeros(3, 3*n)

    for i = 1:n
        p = system.scpfw_parameters
        u0[:, 3*(i-1)+1] = bodies[i].r 
        u0[:, 3*(i-1)+2] = bodies[i].r + [p.rOH, 0.0, 0.0]
        u0[:, 3*(i-1)+3] = bodies[i].r + [cos(p.∠HOH)*p.rOH, 0.0, sin(p.∠HOH)*p.rOH]
        v0[:, 3*(i-1)+1] = bodies[i].v
        v0[:, 3*(i-1)+2] = bodies[i].v
        v0[:, 3*(i-1)+3] = bodies[i].v
    end 

    (u0, v0, n)
end

function gather_accelerations_for_potentials(simulation::NBodySimulation{<:WaterSPCFw})
    o_acelerations = []
    h_acelerations = []
    
    n = length(simulation.system.bodies)

    push!(o_acelerations, get_accelerating_function(simulation.system.e_parameters, simulation))
    push!(o_acelerations, get_accelerating_function(simulation.system.spc_parameters, simulation))
    ms = obtain_data_for_lennard_jones_interaction(simulation.system)
    
    push!(o_acelerations, (dv, u, v, t, i) -> pairwise_lennard_jones_acceleration!(dv[1:3:n], dv[1:3:n], i, length(simulation.system.bodies), ms, simulation.system.lj_parameters, simulation.boundary_conditions))


    push!(h_acelerations, get_accelerating_function(simulation.system.lj_parameters, simulation))
    push!(h_acelerations, get_accelerating_function(simulation.system.spc_parameters, simulation))

    (o_acelerations, h_acelerations)
end