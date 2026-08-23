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

# The SciML common interface that DiffEqPhysics reexports (see the second `export`
# below), so that `using DiffEqPhysics` on its own is enough to build a
# `HamiltonianProblem`, hand it to a solver, drive the integrator from a callback, and
# take apart the result. `HamiltonianProblem` returns an `ODEProblem` wrapping a
# `DynamicalODEFunction` over an `ArrayPartition` state, so the dynamical problem types
# and `ArrayPartition` are what the documented examples index through -- `sol.u[i].x[1]`
# in the orbit recipe, `sol2[2, :]` in `examples/pendulum.jl`. Every name stays owned and
# documented upstream, by SciMLBase and RecursiveArrayTools respectively.
using SciMLBase: CallbackSet, ContinuousCallback, DiscreteCallback,
    DynamicalODEProblem, EnsembleAnalysis, EnsembleDistributed, EnsembleProblem,
    EnsembleSerial, EnsembleSolution, EnsembleSplitThreads, EnsembleSummary,
    EnsembleThreads, ODEFunction, ODESolution, ReturnCode, SecondOrderODEProblem,
    VectorContinuousCallback, init, remake, solve, solve!, step!, successful_retcode,
    terminate!, u_modified!
using RecursiveArrayTools: DiffEqArray, VectorOfArray

include("hamiltonian.jl")
include("plot.jl")

export HamiltonianProblem, orbitplot, plot_orbits

# Reexported SciML common interface; approved via `reexports_allow` in test/qa/qa.jl.
export ArrayPartition, CallbackSet, ContinuousCallback, DiffEqArray, DiscreteCallback,
    DynamicalODEFunction, DynamicalODEProblem, EnsembleAnalysis, EnsembleDistributed,
    EnsembleProblem, EnsembleSerial, EnsembleSolution, EnsembleSplitThreads,
    EnsembleSummary, EnsembleThreads, NullParameters, ODEFunction, ODEProblem,
    ODESolution, ReturnCode, SecondOrderODEProblem, VectorContinuousCallback,
    VectorOfArray, init, remake, solve, solve!, step!, successful_retcode, terminate!,
    u_modified!

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
