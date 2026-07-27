module DiffEqPhysics

using RecursiveArrayTools: ArrayPartition
using DifferentiationInterface: AutoForwardDiff, derivative, gradient, gradient!
import ForwardDiff
using StaticArraysCore: StaticArraysCore
using RecipesBase: @recipe, @series
using SciMLBase: SciMLBase, NullParameters, ODEProblem, DynamicalODEFunction,
    isinplace, AbstractDynamicalODEProblem, numargs, AbstractSciMLSolution,
    TooManyArgumentsError, TooFewArgumentsError, FunctionArgumentsError

include("hamiltonian.jl")
include("plot.jl")

export HamiltonianProblem, orbitplot, plot_orbits

end # module
