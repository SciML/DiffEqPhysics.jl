using StaticArrays, DiffEqBase
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
    k :: AbstractFloat
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
struct NBodySimulation{sType, cType <: AbstractFloat}
    system :: sType
    limiting_boundary :: SVector{6, cType}
    tspan :: Tuple{Float64, Float64}
    external_electric_field
    external_magnetic_field
    external_gravitational_field
end

NBodySimulation(system::NBodySystem, limiting_boundary :: SVector{6, Float64}, tspan :: Tuple{Float64, Float64}) = 
    NBodySimulation(system, limiting_boundary, tspan, x->0, x->0, x->0)

# Instead of treating NBodySimulation as a DiffEq problem and passing it into a solve method
# it is better to use a specific function for n-body simulations.
function run_simulation(s::NBodySimulation, args...; kwargs...)
    #create a DEProblem depending on types of interactions between particles
    solution = solve(problem) 
    return SimulationResult(solution)
end

# Function as this one can be used in a force-oriented approach of calculations
function acceleration(
    system :: NBodySystem, 
    external_electric_field,
    external_magnetic_field,
    external_gravitational_field)
end

# Function, scpecifically designed for a system of ChargedParticles
function acceleration(
    rs,
    i,
    n,
    system ::  ChargedParticles, 
    external_electric_field)

    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j=1:n
        if j!=i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            accel += system.k*system.bodies[i].q*system.bodies[j].q/mi*(ri-rj)/norm(ri-rj)^3
        end
    end
    
    accel += system.bodies[i].q*external_electric_field(ri);

    return accel
end


#= 
    Here u will be an array consisting of a consecutive positional coordinates r 
    following by components of velocities v
=#
function DiffEqBase.ODEProblem(simulation :: NBodySimulation{ChargedParticles, <:AbstractFloat})
    bodies = simulation.system.bodies;
    n = length(bodies)
    u0 = zeros(3, 2*n)
    
    for i=1:n
        u0[:, i] = bodies[i].r 
        u0[:, n+i] = bodies[i].v
    end    
    
    function ode_system!(du, u, p, t)
        du[:, 1:n] = @view u[:, n+1:2n];
        @inbounds for i=1:n
            du[:, n+i] .= acceleration(u[:, 1:n], i, n, simulation.sysmem, simulation.external_electric_field);
        end 
    end

    ODEProblem(ode_system!, u0, simulation.tspan)
end

# SimulationResult sould provide an interface for working with properties of a separate particle
# and with physical properties of the whole system.
struct SimulationResult

end

function SimulationResult(solution :: AbstractTimeseriesSolution)
    #build simulation result out a diff eq solution
end

function temperature(result :: SimulationResult, time :: Float64)
end

function pressure(result :: SimulationResult, time :: Float64)
end

function total_energy(result :: SimulationResult, time :: Float64)
end


# An example of the API usage
p1 = SphericalBody(SVector(-10.0, 0.0, 0.0), SVector(0.0, 10.0, 0.0), 1, 5e-6, 1e-3, SVector(1, 10.0, -2))
p2 = SphericalBody(SVector(10.0, 0.0, 0.0), SVector(0.0, -10.0, 0.0), 1, 15e-7, -4e-3, SVector(0.0, 0.0, 0.0))
system = ChargedParticles([p1, p2], 9e9)
simulation = NBodySimulation(system, SVector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), (0.0, 1.0))
result = run_simulation(simulation)