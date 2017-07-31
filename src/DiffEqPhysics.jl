__precompile__()

module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq, ForwardDiff

include("problem.jl")
include("hamiltonian.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem

end # module
