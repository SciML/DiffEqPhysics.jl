__precompile__()

module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq, RecursiveArrayTools
using ForwardDiff, StaticArrays, RecipesBase

include("nbody.jl")
include("hamiltonian.jl")
include("lagrangian.jl")
include("plot.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem, plot_orbits

end # module
