__precompile__()

module DiffEqPhysics

using DiffEqBase

include("problems.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem

end # module
