module DiffEqPhysics

using PrecompileTools: @compile_workload, @setup_workload
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

@setup_workload begin
    @compile_workload begin
        H(p, q, _, _) = p^2 / 2 + q^2 / 2
        HamiltonianProblem(H, 1.0, 0.0, (0.0, 1.0))

        dp(p, q, _, _) = -q
        dq(p, q, _, _) = p
        HamiltonianProblem((dp, dq), 1.0, 0.0, (0.0, 1.0))
    end
end

end # module
