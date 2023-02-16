module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, RecursiveArrayTools
using ForwardDiff, StaticArraysCore, RecipesBase
using Random, Printf, LinearAlgebra

using DiffEqBase: NullParameters

include("hamiltonian.jl")
include("plot.jl")

export HamiltonianProblem, plot_orbits

end # module
