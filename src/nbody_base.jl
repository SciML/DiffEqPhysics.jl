#Structures for bodies and systems under symulations

# Fields for position, velocity and mass of a particle should be required for every descendant
abstract type Body

end

struct SphericalBody{mType, qType, cType <: AbstractFloat} <: Body
    r :: SVector{3, cType}  #required 
    v :: SVector{3, cType}  #required
    m :: mType              #required
    R :: cType
    q :: qType 
    mm :: SVector{3, qType}
end

abstract type NBodySystem
end

struct MDSystem{bType<:Body} <: NBodySystem
    bodies :: Array{bType, 1}
    thermostats :: Array{Symbol, 1}
    collisions :: Array{Symbol, 1}
end

struct ChargedParticles{bType<:Body} <: NBodySystem
    bodies :: Array{bType, 1}
end

struct GravitationalInteractions{bType<:Body} <: NBodySystem
    bodies :: Array{bType, 1}
end

struct GenericInteractingBodies{bType<:Body} <: NBodySystem
    bodies :: Array{bType, 1}
    pairwise_potential
end

# It is convenient to load data for particles from a file
function load_system_from_file(path::AbstractString)
    f = open(path)
    data = readlines(f)
    n = size(data, 1)
    data = load_file(path);
    bodies = Body[]
    for i=1:n
        #body = compile_body(data[i])
        #push!(bodies, body)
    end
    
    return NBodySystem(bodies)
end

# This structure defines conditions under wich we test our sysmem of n-bodies
struct NBodySimulation{cType <: AbstractFloat}
    limiting_boundary :: SVector{6, cType}
    tspan :: Tuple{Float64, Float64}
    external_electric_field
    external_magnetic_field
    external_gravitational_field
end

NBodySimulation(limiting_boundary :: SVector{6, cType}, tspan :: Tuple{Float64, Float64}) = NBodySimulation(limiting_boundary, tspan, x->0, x->0, x->0)

# Instead of treating NBodySimulation as a DiffEq problem and passing it into a solve method
# it is better to use a specific function for n-body simulations.
function run_simulation(s::NBodySimulation, args...; kwargs...)
    #create a DEProblem depending on types of interactions between particles
    solution = solve(problem) 
    return SimulationResult(solution)
end

# This function can be used in a force-oriented approach of calculations
function acceleration(
    system::NBodySystem, 
    external_electric_field,
    external_magnetic_field,
    external_gravitational_field)
end

struct SimulationResult

end

# SimulationResult sould provide an interface for working with properties of a separate particle
# and with physical properties of the whole system.
function SimulationResult(solution :: AbstractTimeseriesSolution)
    #build simulation result out a diff eq solution
end

function temperature(result :: SimulationResult, time :: Float64)
end

function pressure(result :: SimulationResult, time :: Float64)
end

function total_energy(result :: SimulationResult, time :: Float64)
end