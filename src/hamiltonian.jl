struct HamiltonianProblem{iip} <: DiffEqBase.AbstractDynamicalODEProblem end

function HamiltonianProblem(H,p0,q0,tspan,p=nothing;kwargs...)
  iip = (typeof(q0) <: AbstractArray) && !(typeof(q0) <: SArray)
  HamiltonianProblem{iip}(H,p0,q0,tspan,p;kwargs...)
end

struct PhysicsTag end

function generic_derivative(q0, hami, x)
    ForwardDiff.gradient(hami, x)
end

function generic_derivative(q0::Number, hami, x)
    ForwardDiff.derivative(hami, x)
end

function HamiltonianProblem{false}(H,p0,q0,tspan,p=nothing;kwargs...)
    dp = function (v,x,p,t)
        generic_derivative(q0, x->-H(v, x, p), x)
    end
    dq = function (v,x,p,t)
        generic_derivative(q0, v->H(v, x, p), v)
    end
    return ODEProblem(DynamicalODEFunction{false}(dp,dq), ArrayPartition(p0,q0), tspan, p; kwargs...)
end

function HamiltonianProblem{true}(H,p0,q0,tspan,p=nothing;kwargs...)
    cfg2 = ForwardDiff.GradientConfig(PhysicsTag(), q0)
    dp = function (dv,v,x,p,t)
        fun2 = x->-H(v, x, p)
        ForwardDiff.gradient!(dv, fun2, x, cfg2, Val{false}())
    end
    cfg = ForwardDiff.GradientConfig(PhysicsTag(), p0)
    dq = function (dx,v,x,p,t)
        fun1 = v-> H(v, x, p)
        ForwardDiff.gradient!(dx, fun1, v, cfg, Val{false}())
    end
    return ODEProblem(DynamicalODEFunction{true}(dp,dq), ArrayPartition(p0,q0), tspan, p; kwargs...)
end
