using OrdinaryDiffEq, ForwardDiff

struct HamiltonianProblem{iip} <: AbstractDynamicalODEProblem end

function HamiltonianProblem(H,q0,p0,tspan;kwargs...)
  iip = (typeof(q0) <: AbstractArray) && !(typeof(q0) <: SArray)
  HamiltonianProblem{iip}(H,q0,p0,tspan;kwargs...)
end

struct PhysicsTag end

function HamiltonianProblem{T}(H,q0,p0,tspan;kwargs...) where T
    if T == false
        dq = function (t,x,v)
            ForwardDiff.derivative(v->H(x, v), v)
        end
        dp = function (t,x,v)
            ForwardDiff.derivative(x->-H(x, v), x)
        end

        return ODEProblem{T}(DynamicalODEFunction{T}(dq,dp), (q0,p0), tspan,
                           HamiltonianProblem{false}(); kwargs...)
    else
        cfg = ForwardDiff.GradientConfig(PhysicsTag(), p0)
        dq = function (t,x,v,dx)
            fun1 = v-> H(x, v)
            ForwardDiff.gradient!(dx, fun1, v, cfg, Val{false}())
        end

        cfg2 = ForwardDiff.GradientConfig(PhysicsTag(), q0)
        dp = function (t,x,v,dv)
            fun2 = x->-H(x, v)
            ForwardDiff.gradient!(dv, fun2, x, cfg2, Val{false}())
        end

        return ODEProblem{T}(DynamicalODEFunction{T}(dq,dp), (q0,p0),
                             tspan,HamiltonianProblem{true}(); kwargs...)
    end
end
