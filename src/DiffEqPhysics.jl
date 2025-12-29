module DiffEqPhysics

using Reexport: @reexport
using DiffEqBase: DiffEqBase
using RecursiveArrayTools: RecursiveArrayTools, ArrayPartition
@reexport using DiffEqBase, RecursiveArrayTools
using ForwardDiff: ForwardDiff
using StaticArraysCore: StaticArraysCore
using RecipesBase: RecipesBase, @recipe, @series
using SciMLBase: SciMLBase, NullParameters, ODEProblem, DynamicalODEFunction,
                 isinplace, AbstractDynamicalODEProblem, numargs, DESolution,
                 TooManyArgumentsError, TooFewArgumentsError, FunctionArgumentsError

include("hamiltonian.jl")
include("plot.jl")

export HamiltonianProblem, plot_orbits

end # module
