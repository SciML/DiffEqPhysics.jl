using DiffEqPhysics, ForwardDiff, OrdinaryDiffEq
using StaticArrays, LinearAlgebra
using SafeTestsets, Test

test_solve(prob...) = mapreduce(p->solve(p, VelocityVerlet(), dt=1//2).u, ==, prob)

@safetestset "Hamiltonian Test" begin include("hamiltonian_test.jl") end
@safetestset "N-Body Test" begin include("nbody_test.jl") end
