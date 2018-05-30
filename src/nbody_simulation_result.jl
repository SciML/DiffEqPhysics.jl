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
    show(stream, sr.simulation.system)
    print(stream, "Time steps: ") 
    show(stream, length(sr.solution.t))
    println(stream)
    print(stream, "t: ") 
    show(stream, sr.solution.t)
end

(sr::SimulationResult)(args...; kwargs...) = return sr.solution(args...; kwargs...)

# Instead of treating NBodySimulation as a DiffEq problem and passing it into a solve method
# it is better to use a specific function for n-body simulations.
function run_simulation(s::NBodySimulation, alg_type=Tsit5(), args...; kwargs...)
    solution = solve(ODEProblem(s), alg_type, args...; kwargs...)
    return SimulationResult(solution, s)
end

# this should be a method for integrators designed for the SecondOrderODEProblem (It is worth somehow to sort them from other algorithms)
function run_simulation(s::NBodySimulation, alg_type::Union{VelocityVerlet,DPRKN6,Yoshida6}, args...; kwargs...)
    solution = solve(SecondOrderODEProblem(s), alg_type, args...; kwargs...)
    return SimulationResult(solution, s)
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

function get_position(result::SimulationResult, time_idx::Integer=1, i::Integer=0)
    get_position(result, result.solution.t[time_idx], i)
end

function get_masses(system::NBodySystem)
    n = length(system.bodies)
    masses = zeros(n)
    for i = 1:n
        masses[i] = system.bodies[i].m
    end
    return masses
end

function temperature(result::SimulationResult, time::Real)
    kb = 1.38e-23
    velocities = get_velocity(result, time)
    masses = get_masses(simulation.system)
    temperature = mean(sum(velocities.^2, 1) .* masses) / (3kb)
    return temperature
end

function pressure(sr::SimulationResult, time::Real)
end

function total_energy(sr::SimulationResult, time::Real)
end

function potential_energy()

end

function kinetic_energy(sr::SimulationResult, time::Real)
    vs = get_velocity(sr, time)
    masses = get_masses(sr.simulation.system)
    ke = sum(norm.(sum(vs.^2, 1)) .* masses / 2)
    return ke
end

# In future it seems to be convenient to load data for particles from a file
function load_system_from_file(path::AbstractString)
    f = open(path)
    data = readlines(f)
    n = size(data, 1)
    data = load_file(path);
    bodies = Vector{<:Body}
    for i = 1:n
        #body = compile_body(data[i])
        #push!(bodies, body)
    end
    
    return NBodySystem(bodies)
end

@recipe function g(sr::SimulationResult{<:PotentialNBodySystem}, time::Real=0)
    solution = sr.solution
    n = div(size(solution[1], 2), 2)

    if :gravitational ∈ keys(sr.simulation.system.potentials)
        borders = sr.simulation.boundary_conditions
    
        xlim --> 1.1 * [minimum(solution[1,1:n,:]), maximum(solution[1,1:n,:])]
        ylim --> 1.1 * [minimum(solution[2,1:n,:]), maximum(solution[2,1:n,:])]  
        #zlim --> 1.1 * [minimum(solution[3,1:n,:]), maximum(solution[3,1:n,:])]  
    
        seriestype --> :scatter
        markersize --> 5

        positions = get_position(sr, time)
        (positions[1,:], positions[2,:], positions[3,:])
        #(positions[1,:], positions[2,:])
    else
    
        xlim --> 1.1 * [minimum(solution[1,1:n,:]), maximum(solution[1,1:n,:])]
        ylim --> 1.1 * [minimum(solution[2,1:n,:]), maximum(solution[2,1:n,:])]        
    
        for i in 1:n
            @series begin
                label --> "Orbit $i"
                vars --> (3 * (i - 1) + 1, 3 * (i - 1) + 2)
                solution
            end
        end
    end
end

# This iterator interface is implemented specifically for making animation.
# Probably, there will be a wrapper for this in the future.
Base.start(::SimulationResult) = 1

Base.done(sr::SimulationResult, state) = state > length(sr.solution.t)

function Base.next(sr::SimulationResult, state) 
    positions = get_position(sr, state)

    if sr.simulation.boundary_conditions isa PeriodicBoundaryConditions
        L = sr.simulation.boundary_conditions[2]
        map!(x -> begin while x > L x -= L end; x end, positions, positions)
        map!(x -> begin while x < 0 x += L end; x end, positions, positions)
    end

    (positions[1,:], positions[2,:], positions[3,:]), state + 1
    #(positions[1,:], positions[2,:]), state + 1
end