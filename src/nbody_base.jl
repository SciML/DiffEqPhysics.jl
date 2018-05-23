using StaticArrays, DiffEqBase, OrdinaryDiffEq
#Structures for bodies and systems under symulations

# Fields for position, velocity and mass of a particle should be required for every descendant
abstract type Body
end

abstract type NBodySystem
end

# This structure defines conditions under wich we test our system of n-bodies
struct NBodySimulation{sType <: NBodySystem, cType <: AbstractFloat}
    system :: sType
    limiting_boundary :: SVector{6, cType}
    boundary_conditions :: Symbol
    tspan :: Tuple{Float64, Float64}
    external_electric_field
    external_magnetic_field
    external_gravitational_field
end

NBodySimulation(system::NBodySystem, limiting_boundary :: SVector{6, Float64}, tspan :: Tuple{Float64, Float64}; boundary_conditions=:none) = 
    NBodySimulation(system, limiting_boundary, boundary_conditions, tspan, x->0, x->0, x->0);

# SimulationResult sould provide an interface for working with properties of a separate particle
# and with physical properties of the whole system.
struct SimulationResult
    solution :: AbstractTimeseriesSolution
end

struct SphericalBody{mType, qType, cType <: AbstractFloat} <: Body
    r :: SVector{3, cType}  #required 
    v :: SVector{3, cType}  #required
    m :: mType              #required
    R :: cType
    q :: qType 
    mm :: SVector{3, qType}
end

struct MassBody{mType <: AbstractFloat, cType <: AbstractFloat} <:Body
    r :: SVector{3, cType}
    v :: SVector{3, cType}
    m :: mType     
end

struct ChargedParticle{fType <: AbstractFloat} <:Body
    r :: SVector{3, fType}
    v :: SVector{3, fType}
    m :: fType 
    q :: fType     
end

#=
function SimulationResult(solution :: AbstractTimeseriesSolution)
    #build simulation result out a diff eq solution
end
=#

# Instead of treating NBodySimulation as a DiffEq problem and passing it into a solve method
# it is better to use a specific function for n-body simulations.
function run_simulation(s::NBodySimulation, args...; kwargs...)
    solution = solve(ODEProblem(s), args...; kwargs...)
    return SimulationResult(solution)
end

struct ChargedParticles{bType<:ChargedParticle} <: NBodySystem
    bodies :: Array{bType, 1}
    k :: AbstractFloat
end

struct LennardJonesParameters
    ϵ :: AbstractFloat
    σ :: AbstractFloat
    R :: AbstractFloat
end

struct PotentialNBodySystem <: NBodySystem
    bodies :: Vector{<:Body}
    potentials :: Vector{Symbol}
    lj_parameters :: LennardJonesParameters
end

function PotentialNBodySystem(
    bodies :: Vector{<:Body}; 
    potentials :: Vector{Symbol} = [], 
    lj_epsilon = nothing,
    lj_sigma = nothing,
    lj_range = nothing)
    
    if :lennard_jones ∈ potentials
        if lj_epsilon == nothing || lj_sigma == nothing
            error("The Lennard-Jones potential is switched on for this system. Further correct calculations require `lj_epsilon` and `lj_sigma` to be set explicitly.")
        end
        if lj_range == nothing
            warn("The interaction limit `lj_range` for the Lennard-Jones potential wasn't set. The corresponding forces will be calculated between all the particles of the system.")
        end
    end
    lj_parameters = LennardJonesParameters(lj_epsilon, lj_sigma, lj_range)
    PotentialNBodySystem(bodies, potentials, lj_parameters)
end

function gather_bodies_initial_coordinates(system :: NBodySystem)
    bodies = simulation.system.bodies;
    n = length(bodies)
    u0 = zeros(3, 2*n)
    m = zeros(n)

    for i=1:n
        u0[:, i] = bodies[i].r 
        u0[:, n+i] = bodies[i].v
        m[i] = bodies[i].m
    end 

    (u0, n)
end

#= 
    Here u will be an array consisting of a consecutive positional coordinates r 
    following by components of velocities v
=#
function DiffEqBase.ODEProblem(simulation :: NBodySimulation{<:ChargedParticles, <:AbstractFloat})
    (u0, n) = gather_bodies_initial_coordinates(simulation.system)   
    
    function ode_system!(du, u, p, t)
        du[:, 1:n] = @view u[:, n+1:2n];
        @inbounds for i=1:n
            du[:, n+i] .= pairwise_electrostatic_acceleration(u[:, 1:n], i, simulation.system);
        end 
    end

    return ODEProblem(ode_system!, u0, simulation.tspan)
end

# Function, scpecifically designed for a system of ChargedParticles
function pairwise_electrostatic_acceleration(
    rs,
    i :: Int,
    system :: NBodySystem)

    n = length(system.bodies)
    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j=1:n
        if j!=i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            accel += system.k*system.bodies[i].q*system.bodies[j].q/system.bodies[i].m*(ri-rj)/norm(ri-rj)^3
        end
    end    

    return accel
end

function DiffEqBase.ODEProblem(simulation :: NBodySimulation{PotentialNBodySystem, <:AbstractFloat})
    (u0, n) = gather_bodies_initial_coordinates(simulation.system)
    L = simulation.limiting_boundary[2]

    acceleration_functions = []
    for potential in simulation.system.potentials
        if potential == :lennard_jones
            push!(acceleration_functions, (u,i,system) -> pairwise_lennard_jones_acceleration(u,i,system, simulation.boundary_conditions, simulation.limiting_boundary[2]))
        end
    end

    function ode_system!(du, u, p, t)
        du[:, 1:n] = @view u[:, n+1:2n];

        if simulation.boundary_conditions == :periodic
            for i = 1:n
                for j = 1:3
                    while (u[j,i] > L) u[j,i] -= L end
                    while (u[j,i] < 0) u[j,i] += L end
                end
            end
        end

        @inbounds for i=1:n
            for acceleration in acceleration_functions
                du[:, n+i] .= acceleration(u[:, 1:n], i, simulation.system);                
            end
        end 
    end

    return ODEProblem(ode_system!, u0, simulation.tspan)
end

function pairwise_lennard_jones_acceleration(rs,
    i :: Int,
    system :: NBodySystem,
    boundary_conditions :: Symbol,
    L :: AbstractFloat)

    n = length(system.bodies)
    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j=1:n
        if j!=i
            for add1 = -1:1, add2 = -1:1, add3 = -1:1
                rj = @SVector [rs[1, j]+add1*L, rs[2, j]+add2*L, rs[3, j]+add3*L]
                rij = norm(ri-rj)
                if rij<system.lj_parameters.R
                    accel += (2*(system.lj_parameters.σ/rij)^13-(system.lj_parameters.σ/rij)^7 )*(ri-rj)/rij
                end
            end
        end
    end    

    return -24*system.lj_parameters.ϵ/system.lj_parameters.σ*accel
end

function temperature(result :: SimulationResult, time :: Float64)
end

function pressure(result :: SimulationResult, time :: Float64)
end

function total_energy(result :: SimulationResult, time :: Float64)
end

# In future it seems to be convenient to load data for particles from a file
function load_system_from_file(path::AbstractString)
    f = open(path)
    data = readlines(f)
    n = size(data, 1)
    data = load_file(path);
    bodies = Vector{<:Body}
    for i=1:n
        #body = compile_body(data[i])
        #push!(bodies, body)
    end
    
    return NBodySystem(bodies)
end