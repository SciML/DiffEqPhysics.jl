__precompile__()

module DiffEqPhysics

using OrdinaryDiffEq

include("problems.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem

end # module
