# SimulationResult sould provide an interface for working with properties of a separate particle
# and with physical properties of the whole system.
struct SimulationResult{sType <: NBodySystem}
    solution::AbstractTimeseriesSolution
    simulation::NBodySimulation{sType}
end

function Base.show(stream::IO, sr::SimulationResult)
    print(stream, "N: ") 
    show(stream, length(sr.simulation.system.bodies))
    println(stream)
    show(stream, sr.simulation)
    print(stream, "Time steps: ") 
    show(stream, length(sr.solution.t))
    println(stream)
    print(stream, "t: ", minimum(sr.solution.t), ", ", maximum(sr.solution.t)) 
    println(stream)
end

(sr::SimulationResult)(args...; kwargs...) = return sr.solution(args...; kwargs...)


# This iterator interface is implemented specifically for making animation.
# Probably, there will be a wrapper for this in the future.
Base.start(::SimulationResult) = 1

Base.done(sr::SimulationResult, state) = state > length(sr.solution.t)

function Base.next(sr::SimulationResult, state) 
    positions = get_position(sr, sr.solution.t[state])

    (sr, sr.solution.t[state]), state + 1
end

function get_velocity(sr::SimulationResult, time::Real, i::Integer=0)
    if typeof(sr.solution[1]) <: RecursiveArrayTools.ArrayPartition
        velocities = sr(time).x[1]
        n = size(velocities, 2)
        if i <= 0
            return velocities[:, 1:end]
        else
            return velocities[:, i]
        end
    else
        velocities = sr(time)
        n = div(size(velocities, 2), 2)
        if i <= 0
            return velocities[:, n + 1:end]
        else
            return velocities[:, n + i]
        end
    end
end

function get_position(sr::SimulationResult, time::Real, i::Integer=0)
    if typeof(sr.solution[1]) <: RecursiveArrayTools.ArrayPartition
        positions = sr(time).x[2]
        n = size(positions, 2)
    else
        positions = sr(time)
        n = div(size(positions, 2), 2)
    end

    if i <= 0
        return positions[:, 1:n]
    else
        return positions[:, i]
    end
end

function get_masses(system::NBodySystem)
    n = length(system.bodies)
    masses = zeros(n)
    for i = 1:n
        masses[i] = system.bodies[i].m
    end
    return masses
end

function get_masses(system::WaterSPCFw)
    n = length(system.bodies)
    ms = zeros(Real, 3 * n)
    for i = 1:n
        ms[3 * (i - 1) + 1] = system.mO
        ms[3 * (i - 1) + 2] = system.mH
        ms[3 * (i - 1) + 3] = system.mH
    end 
    return ms
end

function temperature(result::SimulationResult, time::Real)
    kb = 1.38e-23
    velocities = get_velocity(result, time)
    masses = get_masses(result.simulation.system)
    temperature = mean(sum(velocities.^2, 1) .* masses) / (3kb)
    return temperature
end

function temperature(result::SimulationResult{<:WaterSPCFw}, time::Real)
    kb = 1.38e-23
    system = result.simulation.system
    n = length(system.bodies)
    vs = get_velocity(result, time)
    mH2O = system.mO+2*system.mH
    v2 = zeros(n)
    for i = 1:n
        indO, indH1, indH2 = 3 * (i - 1) + 1, 3 * (i - 1) + 2, 3 * (i - 1) + 3
        v_c = (vs[:,indO]*system.mO+vs[:,indH1]*system.mH+vs[:,indH2]*system.mH)/mH2O
        v2 = dot(v_c,v_c)
    end
    temperature = mean(v2) * mH2O / (3kb)
    return temperature
end

function kinetic_energy(velocities, masses)
    ke = sum(dot(vec(sum(velocities.^2, 1)), masses / 2))
    return ke
end

function kinetic_energy(sr::SimulationResult, time::Real)
    vs = get_velocity(sr, time)
    ms = get_masses(sr.simulation.system)
    return kinetic_energy(vs, ms)
end

function potential_energy(coordinates, simulation::NBodySimulation)
    e_potential = 0
    system = simulation.system
    n = length(system.bodies)
    if :lennard_jones ∈ keys(system.potentials)
        p = system.potentials[:lennard_jones]
        (ms, indxs) = obtain_data_for_lennard_jones_interaction(system)
        e_potential += lennard_jones_potential(p, indxs, coordinates, simulation.boundary_conditions)
    end
    e_potential
end

function potential_energy(coordinates, simulation::NBodySimulation{<:WaterSPCFw})
    e_potential = 0
    system = simulation.system
    p = system.lj_parameters
    (ms, indxs) = obtain_data_for_lennard_jones_interaction(system)
    e_potential += lennard_jones_potential(p, indxs, coordinates, simulation.boundary_conditions)

    p = system.e_parameters
    (qs, ms, indx, exclude) = obtain_data_for_electrostatic_interaction(simulation.system)
    e_potential += electrostatic_potential(p, indxs, exclude, qs, coordinates, simulation.boundary_conditions)

    p = system.scpfw_parameters
    (ms, neighbouhoods) = obtain_data_for_harmonic_bond_interaction(simulation.system, p)
    e_potential += harmonic_bonds_potential(p, coordinates, ms, neighbouhoods)

    p = system.scpfw_parameters
    (ms, bonds) = obtain_data_for_valence_angle_harmonic_interaction(simulation.system)
    e_potential += valence_angle_harmonic_potential(coordinates, bonds)
    e_potential
end

function lennard_jones_potential(p::LennardJonesParameters, indxs::Vector{<:Integer}, coordinates, pbc::BoundaryConditions)
    e_lj = 0
    n = length(indxs)
    @inbounds for ind_i = 1:n
        i = indxs[ind_i]
        ri = @SVector [coordinates[1, i], coordinates[2, i], coordinates[3, i]]
        for ind_j = ind_i + 1:n    
            j = indxs[ind_j]          
            rj = @SVector [coordinates[1, j], coordinates[2, j], coordinates[3, j]]
            
            (rij, rij_2, success) = apply_boundary_conditions!(ri, rj, pbc, p.R2)
            
            if success
                σ_rij_6 = (p.σ2 / rij_2)^3
                σ_rij_12 = σ_rij_6^2
                e_lj += (σ_rij_12 - σ_rij_6 )
            end
        end
    end 

    return 4 * p.ϵ * e_lj
end

function electrostatic_potential(p::ElectrostaticParameters, indxs::Vector{<:Integer}, exclude::Dict{Int,Vector{Int}}, qs, rs, pbc::BoundaryConditions)
    e_el = 0

    n = length(indxs)
    @inbounds for ind_i = 1:n
        i = indxs[ind_i]
        ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
        e_el_i = 0
        for ind_j = ind_i + 1:n    
            j = indxs[ind_j] 
            if !in(j, exclude[i])
                rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]

                (rij, rij_2, success) = apply_boundary_conditions!(ri, rj, pbc, p.R2)
                if success
                    e_el_i += qs[j] / norm(ri - rj)
                end
            end
        end    
        e_el += e_el_i * qs[i]
    end

    return e_el * p.k
end

function harmonic_bonds_potential(p::SPCFwParameters,
    rs,
    ms::Vector{<:Real},
    neighborhoods::Dict{Int,Vector{Tuple{Int,Float64}}})

    e_harmonic = 0
    
    @inbounds for (i, neighborhood) ∈ neighborhoods
        ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
        for (j, k) in neighborhood
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            r = norm(rij)
            d = r - p.rOH
            e_harmonic += d^2 * k
        end
    end
    return e_harmonic
end

function valence_angle_harmonic_potential(
    rs,
    bonds::Vector{Tuple{Int, Int, Int, Float64, Float64}}
    )

    e_valence = 0

    for (a,b,c,valence_angle, k) ∈ bonds
        ra = @SVector [rs[1, a], rs[2, a], rs[3, a]]
        rb = @SVector [rs[1, b], rs[2, b], rs[3, b]]
        rc = @SVector [rs[1, c], rs[2, c], rs[3, c]]

        rba = ra - rb
        rbc = rc - rb

        currenct_angle = acos(dot(rba, rbc) / (norm(rba) * norm(rbc)))        
        e_valence += 2 * k * (currenct_angle - valence_angle)^2
    end
    return e_valence
end

function potential_energy(sr::SimulationResult, time::Real)
    e_potential = 0
    coordinates = get_position(sr, time)
    return potential_energy(coordinates, sr.simulation)
end

function total_energy(sr::SimulationResult, time::Real)
    e_kin = kinetic_energy(sr, time)
    e_pot = potential_energy(sr, time)
    e_kin + e_pot
end

function initial_energy(simulation::NBodySimulation)
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation.system)
    return potential_energy(u0, simulation) + kinetic_energy(v0, simulation) 
end

# Instead of treating NBodySimulation as a DiffEq problem and passing it into a solve method
# it is better to use a specific function for n-body simulations.
function run_simulation(s::NBodySimulation, alg_type=Tsit5(), args...; kwargs...)
    solution = solve(ODEProblem(s), alg_type, args...; kwargs...)
    return SimulationResult(solution, s)
end

# this should be a method for integrators designed for the SecondOrderODEProblem (It is worth somehow to sort them from other algorithms)
function run_simulation(s::NBodySimulation, alg_type::Union{VelocityVerlet,DPRKN6,Yoshida6}, args...; kwargs...)
    cb = obtain_callbacks_for_so_ode_problem(s)
    solution = solve(SecondOrderODEProblem(s), alg_type, args...; callback=cb, kwargs...)
    return SimulationResult(solution, s)
end

function obtain_callbacks_for_so_ode_problem(s::NBodySimulation)
    callback_array = Vector{DECallback}()

    if s.thermostat isa AndersenThermostat
        push!(callback_array, get_andersen_thermostating_callback(s))
    end

    return CallbackSet(tuple(callback_array...)...)
end

function get_andersen_thermostating_callback(s::NBodySimulation)
    p = s.thermostat::AndersenThermostat
    n = length(s.system.bodies)
    v_dev = sqrt(p.kb * p.T / s.system.bodies[1].m)

    condition = function (u, t, integrator)
        true
    end
    affect! = function (integrator)
        for i = 1:n
            if randn() < p.ν * (integrator.t - integrator.tprev)
                integrator.u.x[1][:,i] .= v_dev * randn(3)
            end
        end
    end
    cb = DiscreteCallback(condition, affect!)
end

@recipe function generate_data_for_scatter(sr::SimulationResult{<:PotentialNBodySystem}, time::Real=0.0)
    solution = sr.solution
    n = length(sr.simulation.system.bodies)

    if :gravitational ∈ keys(sr.simulation.system.potentials)
    
        xlim --> 1.1 * [minimum(solution[1,1:n,:]), maximum(solution[1,1:n,:])]
        ylim --> 1.1 * [minimum(solution[2,1:n,:]), maximum(solution[2,1:n,:])]        
    
        for i in 1:n
            @series begin
                label --> "Orbit $i"
                vars --> (3 * (i - 1) + 1, 3 * (i - 1) + 2)
                solution
            end
        end
    else
        borders = sr.simulation.boundary_conditions
        if borders isa PeriodicBoundaryConditions
            xlim --> 1.1 * [borders[1], borders[2]]
            ylim --> 1.1 * [borders[3], borders[4]]
            zlim --> 1.1 * [borders[5], borders[6]]
        elseif borders isa CubicPeriodicBoundaryConditions
            xlim --> 1.1 * [0, borders.L]
            ylim --> 1.1 * [0, borders.L]
            zlim --> 1.1 * [0, borders.L]
        end
        positions = get_position(sr, time)
            
        seriestype --> :scatter
        markersize --> 5

        (positions[1,:], positions[2,:], positions[3,:])
        #(positions[1,:], positions[2,:])
    end
end

function distancies(result::SimulationResult, time::Real)
    n = length(result.simulation.system.bodies)
    cc = get_position(result, time)

    d = Float64[]
    for i = 1:n
        for j = 1:n
            if i != j
                push!(d, norm(vec(cc[:,i] - cc[:,j])))
            end
        end
    end
    return d
end

@recipe function initial_distribution(sr::SimulationResult{<:WaterSPCFw}, time::Real=0.0)

    n = length(sr.simulation.system.bodies)

    borders = sr.simulation.boundary_conditions
    
    cc = get_position(sr, time)
 
    if borders isa PeriodicBoundaryConditions
        xlim --> 1.1 * [borders[1], borders[2]]
        ylim --> 1.1 * [borders[3], borders[4]]
        zlim --> 1.1 * [borders[5], borders[6]]
    elseif borders isa CubicPeriodicBoundaryConditions
        xlim --> 1.1 * [0, borders.L]
        ylim --> 1.1 * [0, borders.L]
        zlim --> 1.1 * [0, borders.L]

        
        map!(x ->  x -= borders.L * floor(x / borders.L), cc, cc)
    end
    seriestype --> :scatter

    @series begin
        label --> "O"
        markersize --> 8
        markercolor --> :red
        (cc[1,1:3:3 * n - 2], cc[2,1:3:3 * n - 2], cc[3,1:3:3 * n - 2])
    end

    @series begin
        label --> "H"
        markersize --> 4
        markercolor --> :green
        x = vcat(cc[1,2:3:3 * n - 1], cc[1,3:3:3 * n])
        y = vcat(cc[2,2:3:3 * n - 1], cc[2,3:3:3 * n])
        z = vcat(cc[3,2:3:3 * n - 1], cc[3,3:3:3 * n])
        (x, y, z)
    end
end