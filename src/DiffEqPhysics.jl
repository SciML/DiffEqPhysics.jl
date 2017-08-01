__precompile__()

module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq, ForwardDiff

include("nbody.jl")
include("hamiltonian.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem

end # module
