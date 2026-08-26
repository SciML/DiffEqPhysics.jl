using SciMLTesting, DiffEqPhysics, Test

# The SciML common interface DiffEqPhysics deliberately reexports so that `using
# DiffEqPhysics` is enough to build a `HamiltonianProblem`, solve it, and take the
# resulting `ArrayPartition` solution apart, as the README and `examples/` do. Owned and
# documented upstream (SciMLBase, RecursiveArrayTools); kept in sync with the reexport
# `export` block in src/DiffEqPhysics.jl.
const REEXPORTS = (
    :ArrayPartition, :CallbackSet, :ContinuousCallback, :DiffEqArray, :DiscreteCallback,
    :DynamicalODEFunction, :DynamicalODEProblem, :EnsembleAnalysis, :EnsembleDistributed,
    :EnsembleProblem, :EnsembleSerial, :EnsembleSolution, :EnsembleSplitThreads,
    :EnsembleSummary, :EnsembleThreads, :NullParameters, :ODEFunction, :ODEProblem,
    :ODESolution, :ReturnCode, :SecondOrderODEProblem, :VectorContinuousCallback,
    :VectorOfArray, :init, :remake, :solve, :solve!, :step!, :successful_retcode,
    :terminate!, :u_modified!,
)

run_qa(DiffEqPhysics; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using DiffEqPhysics`, so
    # the allow-list cannot drift into approving names the package no longer provides.
    @testset "$name" for name in REEXPORTS
        @test name in names(DiffEqPhysics)
        @test isdefined(@__MODULE__, name)
    end
end
