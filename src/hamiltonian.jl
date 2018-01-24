using OrdinaryDiffEq, ForwardDiff

struct HamiltonianProblem{iip} <: AbstractDynamicalODEProblem end

function HamiltonianProblem(H,q0,p0,tspan,p=nothing;kwargs...)
  iip = (typeof(q0) <: AbstractArray) && !(typeof(q0) <: SArray)
  HamiltonianProblem{iip}(H,q0,p0,tspan,p;kwargs...)
end

struct PhysicsTag end

function HamiltonianProblem{T}(H,q0,p0,tspan,p=nothing;kwargs...) where T
    if T == false
        if typeof(q0) <: Number
            dq = function (x,v,p,t)
                ForwardDiff.derivative(v->H(x, v, p), v)
            end
            dp = function (x,v,p,t)
                ForwardDiff.derivative(x->-H(x, v, p), x)
            end

            return ODEProblem{T}(DynamicalODEFunction{T}(dq,dp), (q0,p0), tspan,
                               p,HamiltonianProblem{false}(); kwargs...)
        else
            dq = function (x,v,p,t)
                ForwardDiff.gradient(v->H(x, v, p), v)
            end
            dp = function (x,v,p,t)
                ForwardDiff.gradient(x->-H(x, v, p), x)
            end

            return ODEProblem{T}(DynamicalODEFunction{T}(dq,dp), (q0,p0), tspan,
                               p,HamiltonianProblem{false}(); kwargs...)
        end
    else
        let cfg = ForwardDiff.GradientConfig(PhysicsTag(), p0), cfg2 = ForwardDiff.GradientConfig(PhysicsTag(), q0)
            dq = function (dx,x,v,p,t)
                fun1 = v-> H(x, v, p)
                ForwardDiff.gradient!(dx, fun1, v, cfg, Val{false}())
            end
            dp = function (dv,x,v,p,t)
                fun2 = x->-H(x, v, p)
                ForwardDiff.gradient!(dv, fun2, x, cfg2, Val{false}())
            end
            return ODEProblem{T}(DynamicalODEFunction{T}(dq,dp), (q0,p0),
                                 tspan,p,HamiltonianProblem{true}(); kwargs...)
        end
    end
end
