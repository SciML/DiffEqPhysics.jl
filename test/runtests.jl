using DiffEqPhysics, ForwardDiff, OrdinaryDiffEq
using StaticArrays, LinearAlgebra
using Test

test_solve(prob...) = mapreduce(p->solve(p, VelocityVerlet(), dt=1//2).u, ==, prob)
include("hamiltonian_test.jl")
include("nbody_test.jl")
