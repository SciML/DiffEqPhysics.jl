__precompile__()

module DiffEqPhysics

using OrdinaryDiffEq, ForwardDiff

include("problems.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem

end # module
