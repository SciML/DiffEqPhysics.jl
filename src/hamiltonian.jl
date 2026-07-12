const TOO_MANY_ARGUMENTS_ERROR_MESSAGE = """
All methods for the model function `H` had too many arguments. A
Hamiltonian `H` must define either `H(p, q, param)` or `H(p, q, param, t)`. This error
can be thrown if you define an Hamiltonian for example as `H(p, q, param1, param2, t)`.
For more information on the required number of arguments for the function
you were defining, consult the documentation for the `HamiltonianProblem`.
"""

struct HamiltonianTooManyArgumentsError <: Exception
    fname::String
    f::Any
end

function Base.showerror(io::IO, e::HamiltonianTooManyArgumentsError)
    println(io, TOO_MANY_ARGUMENTS_ERROR_MESSAGE)
    print(io, "Offending function: ")
    printstyled(io, e.fname; bold = true, color = :red)
    println(io, "\nMethods:")
    return println(io, methods(e.f))
end

const TOO_FEW_ARGUMENTS_ERROR_MESSAGE = """
All methods for the Hamiltonian `H` had too few arguments. A
Hamiltonian `H` must define either `H(p, q, param)` or `H(p, q, param, t)`. This error
can be thrown if you define an Hamiltonian for example as `H(p, q)`.
Note that `param` must be in the arguments list even if it's not used.
For more information on the required number of arguments for the function
you were defining, consult the documentation for the `HamiltonianProblem`.
"""

struct HamiltonianTooFewArgumentsError <: Exception
    fname::String
    f::Any
end

function Base.showerror(io::IO, e::HamiltonianTooFewArgumentsError)
    println(io, TOO_FEW_ARGUMENTS_ERROR_MESSAGE)
    print(io, "Offending function: ")
    printstyled(io, e.fname; bold = true, color = :red)
    println(io, "\nMethods:")
    return println(io, methods(e.f))
end

const ARGUMENTS_ERROR_MESSAGE = """
Methods dispatches for the Hamiltonian `H` do not match the required number.
A Hamiltonian `H` must define either `H(p, q, param)` or `H(p, q, param, t)`. This error
can be thrown if you define an Hamiltonian for example as `H(p)`.
Note that arguments must be in the arguments list even if it's not used.
For more information on the required number of arguments for the function
you were defining, consult the documentation for the `HamiltonianProblem`.
"""

struct HamiltonianFunctionArgumentsError <: Exception
    fname::String
    f::Any
end

function Base.showerror(io::IO, e::HamiltonianFunctionArgumentsError)
    println(io, ARGUMENTS_ERROR_MESSAGE)
    print(io, "Offending function: ")
    printstyled(io, e.fname; bold = true, color = :red)
    println(io, "\nMethods:")
    return println(io, methods(e.f))
end

struct HamiltonianProblem{iip} <: AbstractDynamicalODEProblem end

"""
    HamiltonianProblem(H, p0, q0, tspan, param=nothing; kwargs...)
    HamiltonianProblem((dp, dq), p0, q0, tspan, param=nothing; kwargs...)

Define a Hamiltonian system as a SciML dynamical ODE problem.

The Hamiltonian form evolves canonical momenta `p` and coordinates `q` according to

```math
\\dot{q} = \\partial H / \\partial p, \\qquad
\\dot{p} = -\\partial H / \\partial q.
```

# Arguments

- `H`: Hamiltonian function. It should support `H(p, q, param, t)`. The deprecated
  three-argument form `H(p, q, param)` is still accepted.
- `(dp, dq)`: Pair of derivative functions, where `dp` computes `-∂H/∂q` and `dq`
  computes `∂H/∂p`. This form skips automatic differentiation.
- `p0`: Initial canonical momentum. Scalars, static arrays, and other
  `AbstractArray`s are supported.
- `q0`: Initial canonical coordinate with a shape compatible with `p0`.
- `tspan`: Time span passed to the generated `ODEProblem`.
- `param`: Parameter object passed through to `H`, `dp`, and `dq`. If omitted,
  `NullParameters()` is used.

# Keyword Arguments

- `kwargs...`: Additional keyword arguments forwarded to the generated `ODEProblem`.

# Returns

- An `ODEProblem` containing a `DynamicalODEFunction`.

# Interface

- Scalar and static-array states use out-of-place derivative functions with
  signatures `dp(p, q, param, t)` and `dq(p, q, param, t)`.
- Mutable `AbstractArray` states use in-place derivative functions with signatures
  `dp!(dp, p, q, param, t)` and `dq!(dq, p, q, param, t)`.
- The automatic-differentiation constructor differentiates `H` with ForwardDiff.jl.

!!! note
    Prefer `H(p, q, param, t)` for new code. `H(p, q, param)` is deprecated and
    supported only for backward compatibility.

# Examples

```julia
using DiffEqPhysics

H(p, q, params, t) = p^2 / 2 + params.k * q^2 / 2
prob = HamiltonianProblem(H, 1.0, 0.0, (0.0, 10.0), (; k = 2.0))
```
"""
function HamiltonianProblem(
        H, p0::S, q0::T, tspan, param = NullParameters(); kwargs...
    ) where {S, T}
    iip = T <: AbstractArray && !(T <: StaticArraysCore.SArray) && S <: AbstractArray &&
        !(S <: StaticArraysCore.SArray)
    return HamiltonianProblem{iip}(H, p0, q0, tspan, param; kwargs...)
end

struct PhysicsTag end

function generic_derivative(q0, hami, x)
    return ForwardDiff.gradient(hami, x)
end

function generic_derivative(q0::Number, hami, x)
    return ForwardDiff.derivative(hami, x)
end

function HamiltonianProblem{false}(
        (dp, dq)::Tuple{Any, Any}, p0, q0, tspan, param = NullParameters(); kwargs...
    )
    return ODEProblem(
        DynamicalODEFunction{false}(dp, dq),
        ArrayPartition(p0, q0), tspan, param; kwargs...
    )
end
function HamiltonianProblem{true}(
        (dp, dq)::Tuple{Any, Any}, p0, q0, tspan, param = NullParameters(); kwargs...
    )
    return ODEProblem(
        DynamicalODEFunction{true}(dp, dq), ArrayPartition(p0, q0), tspan, param; kwargs...
    )
end

function HamiltonianProblem{false}(H, p0, q0, tspan, param = NullParameters(); kwargs...)
    try
        isinplace(H, 4)
    catch e
        if e isa TooManyArgumentsError
            throw(HamiltonianTooManyArgumentsError(e.fname, e.f))
        elseif e isa TooFewArgumentsError
            throw(HamiltonianTooFewArgumentsError(e.fname, e.f))
        elseif e isa FunctionArgumentsError
            throw(HamiltonianFunctionArgumentsError(e.fname, e.f))
        end
    end

    if 4 in numargs(H)
        dp = (p, q, param, t) -> generic_derivative(q0, q -> -H(p, q, param, t), q)
        dq = (p, q, param, t) -> generic_derivative(q0, p -> H(p, q, param, t), p)
    else
        issue_depwarn()
        dp = (p, q, param, t) -> generic_derivative(q0, q -> -H(p, q, param), q)
        dq = (p, q, param, t) -> generic_derivative(q0, p -> H(p, q, param), p)
    end
    return HamiltonianProblem{false}((dp, dq), p0, q0, tspan, param; kwargs...)
end

function HamiltonianProblem{true}(H, p0, q0, tspan, param = NullParameters(); kwargs...)
    try
        isinplace(H, 4)
    catch e
        if e isa TooManyArgumentsError
            throw(HamiltonianTooManyArgumentsError(e.fname, e.f))
        elseif e isa TooFewArgumentsError
            throw(HamiltonianTooFewArgumentsError(e.fname, e.f))
        elseif e isa FunctionArgumentsError
            throw(HamiltonianFunctionArgumentsError(e.fname, e.f))
        end
    end
    let cp = ForwardDiff.GradientConfig(PhysicsTag(), p0),
            cq = ForwardDiff.GradientConfig(PhysicsTag(), q0), vfalse = Val(false)

        if 4 in numargs(H)
            dp = (
                Δp, p, q, param,
                t,
            ) -> ForwardDiff.gradient!(Δp, q -> -H(p, q, param, t), q, cq, vfalse)
            dq = (
                Δq, p, q, param,
                t,
            ) -> ForwardDiff.gradient!(Δq, p -> H(p, q, param, t), p, cp, vfalse)
        else
            issue_depwarn()
            dp = (
                Δp, p, q, param,
                t,
            ) -> ForwardDiff.gradient!(Δp, q -> -H(p, q, param), q, cq, vfalse)
            dq = (
                Δq, p, q, param,
                t,
            ) -> ForwardDiff.gradient!(Δq, p -> H(p, q, param), p, cp, vfalse)
        end
        return HamiltonianProblem{true}((dp, dq), p0, q0, tspan, param; kwargs...)
    end
end

function issue_depwarn()
    return Base.depwarn(
        "Hamiltonians with 3 arguments are deprecated; please use `H(p, q, params, t)`",
        :HamiltonianProblem
    )
end
