struct HamiltonianProblem{iip} <: DiffEqBase.AbstractDynamicalODEProblem end

function HamiltonianProblem(H, p0, q0, tspan, param=nothing; kwargs...)
  iip = (typeof(q0) <: AbstractArray) && !(typeof(q0) <: SArray)
  HamiltonianProblem{iip}(H, p0, q0, tspan, param; kwargs...)
end

struct PhysicsTag end

function generic_derivative(q0, hami, x)
    ForwardDiff.gradient(hami, x)
end

function generic_derivative(q0::Number, hami, x)
    ForwardDiff.derivative(hami, x)
end

function HamiltonianProblem{false}(H, p0, q0, tspan, param=nothing; kwargs...)
    dp = function (p, q, param, t)
        generic_derivative(q0, q -> -H(p, q, param), q)
    end
    dq = function (p, q, param, t)
        generic_derivative(q0, p -> H(p, q, param), p)
    end
    return ODEProblem(DynamicalODEFunction{false}(dp, dq), ArrayPartition(p0, q0), tspan, param; kwargs...)
end

function HamiltonianProblem{true}(H, p0, q0, tspan, param=nothing; kwargs...)
    let cfg = ForwardDiff.GradientConfig(PhysicsTag(), p0), cfg2 = ForwardDiff.GradientConfig(PhysicsTag(), q0)
        dp = function (dv, p, q, param, t)
            fun2 = q -> -H(p, q, param)
            ForwardDiff.gradient!(dv, fun2, q, cfg2, Val{false}())
        end
        dq = function (dx, p, q, param, t)
            fun1 = p-> H(p, q, param)
            ForwardDiff.gradient!(dx, fun1, p, cfg, Val{false}())
        end
        return ODEProblem(DynamicalODEFunction{true}(dp, dq), ArrayPartition(p0, q0), tspan, param; kwargs...)
    end
end
