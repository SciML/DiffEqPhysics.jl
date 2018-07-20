using DiffEqPhysics, ForwardDiff, StaticArrays, LinearAlgebra
using Test

test_solve(prob...) = mapreduce(p->solve(p, VelocityVerlet(), dt=1//2).u, ==, prob)
include("hamiltonian_test.jl")
include("nbody_test.jl")
include("nbody_lennard_jones_test.jl")
include("nbody_electrostatics_test.jl")
include("nbody_gravitational_test.jl")
include("nbody_custom_potential_test.jl")
include("nbody_magnetostaic_test.jl")
include("nbody_thermostat_test.jl")
include("nbody_water_test.jl")
