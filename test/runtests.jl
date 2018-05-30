using DiffEqPhysics, ForwardDiff
using Base.Test

test_solve(prob...) = mapreduce(p->solve(p, VelocityVerlet(), dt=1//2).u, ==, prob)

include("hamiltonian_test.jl")
include("nbody_test.jl")
include("./../src/nbody_simulation.jl")
include("nbody_lennard_jones_test.jl")
include("nbody_electrostatics_test.jl")
include("nbody_gravitational_test.jl")
include("nbody_custom_potential_test.jl")
include("nbody_magnetostaic_test.jl")
