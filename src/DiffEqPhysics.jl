__precompile__()

module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq
using ForwardDiff, StaticArrays, Plots

include("nbody.jl")
include("hamiltonian.jl")
include("plot.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem, plot_orbits

end # module
