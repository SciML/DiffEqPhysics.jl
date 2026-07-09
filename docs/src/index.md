# DiffEqPhysics.jl

DiffEqPhysics.jl provides physics-based problem types for the [SciML](https://sciml.ai)
ecosystem. The main feature is [`HamiltonianProblem`](@ref), which allows you to define and
solve Hamiltonian systems using automatic differentiation.

For N-body gravitational simulations, please use
[NBodySimulator.jl](https://github.com/SciML/NBodySimulator.jl).

## Installation

```julia
using Pkg
Pkg.add("DiffEqPhysics")
```

## Quick Start: Simple Pendulum

Define a Hamiltonian system by specifying the Hamiltonian function `H(p, q, params, t)`:

```julia
using DiffEqPhysics, OrdinaryDiffEq

# Hamiltonian for a simple pendulum: H = p²/(2ml²) + mgl(1 - cos(q))
function H(p, q, params, t)
    g, m, l = params
    return p^2 / (2 * m * l^2) + m * g * l * (1 - cos(q))
end

# Parameters: gravitational acceleration, mass, length
params = (9.81, 1.0, 1.0)

# Initial conditions: momentum p₀ and angle q₀
p₀ = 0.0
q₀ = π / 4  # 45 degrees

# Create and solve the Hamiltonian problem
prob = HamiltonianProblem(H, p₀, q₀, (0.0, 10.0), params)
sol = solve(prob, SofSpa10(), dt = 0.01)
```

The equations of motion `q̇ = ∂H/∂p` and `ṗ = -∂H/∂q` are automatically computed using
ForwardDiff.jl.

## Features

  - **Automatic differentiation**: Derivatives are computed automatically from the
    Hamiltonian function.
  - **Multiple input types**: Works with scalars, `StaticArrays`, or regular
    `AbstractArray`s.
  - **Symplectic solvers**: Use symplectic integrators from OrdinaryDiffEq.jl (like
    `SofSpa10`) to preserve the Hamiltonian structure.
  - **Manual derivatives**: Optionally provide your own derivative functions `(dp, dq)` for
    better performance.

## Providing Manual Derivatives

For better performance, you can provide the derivative functions directly:

```julia
# dp = -∂H/∂q, dq = ∂H/∂p
dp(p, q, params, t) = -params[1] * params[2] * params[3] * sin(q)  # -mgl*sin(q)
dq(p, q, params, t) = p / (params[2] * params[3]^2)                # p/(ml²)

prob = HamiltonianProblem((dp, dq), p₀, q₀, (0.0, 10.0), params)
```

## Plotting Orbits

For solutions of N-body-style problems, [`plot_orbits`](@ref) and [`orbitplot`](@ref) draw
the trajectory of each body.

## Public API

See the [API](@ref API) page for the full list of documented public functions and types.
