struct HamiltonianProblem{iip} <: DiffEqBase.AbstractDynamicalODEProblem end

function HamiltonianProblem(H,q0,p0,tspan,p=nothing;kwargs...)
  iip = (typeof(q0) <: AbstractArray) && !(typeof(q0) <: SArray)
  HamiltonianProblem{iip}(H,q0,p0,tspan,p;kwargs...)
end

struct PhysicsTag end

function HamiltonianProblem{T}(H,p0,q0,tspan,p=nothing;kwargs...) where T
    if T == false
        if typeof(q0) <: Number
            dp = function (v,x,p,t)
                ForwardDiff.derivative(x->-H(v, x, p), x)
            end
            dq = function (v,x,p,t)
                ForwardDiff.derivative(v->H(v, x, p), v)
            end

            return ODEProblem(DynamicalODEFunction{T}(dp,dq), ArrayPartition(p0,q0), tspan, p; kwargs...)
        else
            dp = function (v,x,p,t)
                ForwardDiff.gradient(x->-H(v, x, p), x)
            end
            dq = function (v,x,p,t)
                ForwardDiff.gradient(v->H(v, x, p), v)
            end

            return ODEProblem(DynamicalODEFunction{T}(dp,dq), ArrayPartition(p0,q0), tspan, p; kwargs...)
        end
    else
        let cfg = ForwardDiff.GradientConfig(PhysicsTag(), p0), cfg2 = ForwardDiff.GradientConfig(PhysicsTag(), q0)
            dp = function (dv,v,x,p,t)
                fun2 = x->-H(v, x, p)
                ForwardDiff.gradient!(dv, fun2, x, cfg2, Val{false}())
            end
            dq = function (dx,v,x,p,t)
                fun1 = v-> H(v, x, p)
                ForwardDiff.gradient!(dx, fun1, v, cfg, Val{false}())
            end
            return ODEProblem(DynamicalODEFunction{T}(dp,dq), ArrayPartition(p0,q0), tspan, p; kwargs...)
        end
    end
end
