struct HamiltonianProblem{iip} <: DiffEqBase.AbstractDynamicalODEProblem end

"""
    HamiltonianProblem(H, p0, q0, tspan, param=nothing; kwargs...)
    HamiltonianProblem((dp, dq), p0, q0, tspan, param=nothing; kwargs...)

Define a physical system by its Hamiltonian function `H(p, q, param)` or the function pair
`dp = -∂H/∂q` and `dq = ∂H/∂p`.

The equations of motion are then given by `q̇ = ∂H/∂p = dq` and `ṗ  = -∂H/∂q = dp`.

The initial values for canonical impulses `p0` and coordinates `q0` may be scalars, `SArray`s or
other `AbstractArray`s. Their type determines the type of functions `dp` and `dq`.

In the first two cases the derivative functions require the signatures
`dp(p, q, param, t)` and `dq(p, q, param, t)` while
in the latter case the partial derivatives use mutating signatures
`dp!(Δp, p, q, param, t)` and `dq!(Δq, p, q, param, t)` with predefined arrays `Δp` and `Δq`.

If the Hamiltonian function is given, `dp` and `dq` are calculated automatically using
AD (`ForwardDiff`).


!!! note
It is assumed, that `H` does not depend on `t`, while `dp` and `dq` do.
It is possible to convert a time-dependent problem given by `H(p, q, param, t)` into a
time-independent one
by adding coordinates `q₀ = t` and `p₀ = E` with the new function
`H([t,p], [E,q], params) := H(p, q, params, t) + E`. 
"""
function HamiltonianProblem(H, p0::T, q0::T, tspan, param=nothing; kwargs...) where T
  iip = T <: AbstractArray && !(T <: SArray)
  HamiltonianProblem{iip}(H, p0, q0, tspan, param; kwargs...)
end

struct PhysicsTag end

function generic_derivative(q0, hami, x)
    ForwardDiff.gradient(hami, x)
end

function generic_derivative(q0::Number, hami, x)
    ForwardDiff.derivative(hami, x)
end

function HamiltonianProblem{false}((dp, dq)::Tuple{Any,Any}, p0, q0, tspan, param=nothing; kwargs...)
    return ODEProblem(DynamicalODEFunction{false}(dp, dq), ArrayPartition(p0, q0), tspan, param; kwargs...)
end
function HamiltonianProblem{true}((dp, dq)::Tuple{Any,Any}, p0, q0, tspan, param=nothing; kwargs...)
    return ODEProblem(DynamicalODEFunction{true}(dp, dq), ArrayPartition(p0, q0), tspan, param; kwargs...)
end

function HamiltonianProblem{false}(H, p0, q0, tspan, param=nothing; kwargs...)
    dp = function (p, q, param, t)
        generic_derivative(q0, q -> -H(p, q, param), q)
    end
    dq = function (p, q, param, t)
        generic_derivative(q0, p -> H(p, q, param), p)
    end
    return HamiltonianProblem{false}((dp, dq), p0, q0, tspan, param; kwargs...)
end

function HamiltonianProblem{true}(H, p0, q0, tspan, param=nothing; kwargs...)
    let cfg = ForwardDiff.GradientConfig(PhysicsTag(), p0), cfg2 = ForwardDiff.GradientConfig(PhysicsTag(), q0)
        dp = function (Δp, p, q, param, t)
            fun2 = q -> -H(p, q, param)
            ForwardDiff.gradient!(Δp, fun2, q, cfg2, Val{false}())
        end
        dq = function (Δq, p, q, param, t)
            fun1 = p-> H(p, q, param)
            ForwardDiff.gradient!(Δq, fun1, p, cfg, Val{false}())
        end
        return HamiltonianProblem{true}((dp, dq), p0, q0, tspan, param; kwargs...)
    end
end
