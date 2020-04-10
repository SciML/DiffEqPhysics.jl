using DiffEqPhysics, ForwardDiff, OrdinaryDiffEq
using StaticArrays, LinearAlgebra
using SafeTestsets, Test

@safetestset "Hamiltonian Test" begin include("hamiltonian_test.jl") end
#@safetestset "N-Body Test" begin include("nbody_test.jl") end
