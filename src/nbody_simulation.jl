using StaticArrays, DiffEqBase, OrdinaryDiffEq, RecipesBase

include("./nbody_bodies.jl")
include("./nbody_basic_potentials.jl")
include("./nbody_system.jl")
include("./nbody_boundary_conditions.jl")


# This structure defines conditions under wich we test our system of n-bodies
# With this wrapping we can make such fields as `boundary_conditions` necessary for every simulation
# while allowing one to describe a particular systme of N interacting particles
struct NBodySimulation{sType <: NBodySystem,bcType <: BoundaryConditions}
    system::sType
    tspan::Tuple{Float64,Float64}
    boundary_conditions::bcType
    external_electric_field
    external_magnetic_field
    external_gravitational_field
end

function NBodySimulation(system::BasicPotentialSystem,
    tspan::Tuple{Float64,Float64},
    boundary_conditions::BoundaryConditions,
    external_electric_field,
    external_magnetic_field,
    external_gravitational_field)
    
    potential_system = PotentialNBodySystem(system)
    NBodySimulation(potential_system, tspan, boundary_conditions, external_electric_field, external_magnetic_field, external_gravitational_field)
end

NBodySimulation(system::NBodySystem, tspan::Tuple{Float64,Float64}, boundary_conditions::BoundaryConditions) = 
    NBodySimulation(system, tspan, boundary_conditions, x -> 0, x -> 0, x -> 0)

NBodySimulation(system::NBodySystem, tspan::Tuple{Float64,Float64}) = 
    NBodySimulation(system, tspan, InfiniteBox(), x -> 0, x -> 0, x -> 0)

function Base.show(stream::IO, s::NBodySimulation)
    print(stream, "Timespan: ")
    show(stream, s.tspan)
    println(stream)
    print(stream, "Boundary conditions: ")
    show(stream, s.boundary_conditions)
    println(stream)
    show(stream, s.system)
end

include("./nbody_to_ode.jl")
include("./nbody_simulation_result.jl")