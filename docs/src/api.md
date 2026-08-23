# [API](@id API)

The following functions and types make up the public API of DiffEqPhysics.jl.

## Problem Types

```@docs
HamiltonianProblem
```

## Plotting

```@docs
plot_orbits
orbitplot
```

## Reexported SciML common interface

`using DiffEqPhysics` also brings in the parts of the SciML common interface needed to
hand a [`HamiltonianProblem`](@ref) to a solver and take the result apart, so they do not
have to be imported separately. There are two owners, and each name stays documented by
the package that owns it.

### From [SciMLBase](https://docs.sciml.ai/SciMLBase/stable/)

`HamiltonianProblem` builds an `ODEProblem` wrapping a `DynamicalODEFunction`, so these
are the problem, solving and solution names that go with it:

  - Problems: `ODEProblem`, `DynamicalODEProblem`, `SecondOrderODEProblem`,
    `EnsembleProblem`
  - Functions: `ODEFunction`, `DynamicalODEFunction`
  - Solutions: `ODESolution`, `EnsembleSolution`, `EnsembleSummary`
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`, and the `EnsembleAnalysis` module
  - Solving: `solve`, `solve!`, `init`, `step!`, `remake`
  - Integrator interface: `u_modified!`, `terminate!`
  - Return status: `ReturnCode`, `successful_retcode`
  - Callbacks: `ContinuousCallback`, `DiscreteCallback`, `VectorContinuousCallback`,
    `CallbackSet`
  - `NullParameters`

### From [RecursiveArrayTools](https://docs.sciml.ai/RecursiveArrayTools/stable/)

A `HamiltonianProblem` state is an `ArrayPartition` of `(p, q)`, and the saved solution
is a `VectorOfArray` of those, which is how the documented examples index a solution
(`sol.u[i].x[1]` in the orbit recipe, `sol[2, :]` in `examples/pendulum.jl`):

  - `ArrayPartition`, `VectorOfArray`, `DiffEqArray`

DiffEqPhysics defines no solvers of its own — the algorithm passed to `solve` comes from
a solver package such as
[OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/), which must be loaded
separately. Anything else from SciMLBase or RecursiveArrayTools must be imported from
that package directly.
