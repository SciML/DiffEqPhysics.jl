__precompile__()

module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq, RecursiveArrayTools
using ForwardDiff, StaticArrays, RecipesBase

include("nbody_problem.jl")
include("hamiltonian.jl")
include("plot.jl")
include("nbody_simulation.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem, plot_orbits

export NBodyGravProblem, MassBody, GravImagingData, solve


end # module
