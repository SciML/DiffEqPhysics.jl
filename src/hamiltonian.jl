using OrdinaryDiffEq, ForwardDiff

struct HamiltonianProblem{iip} <: AbstractDynamicalODEProblem end

function HamiltonianProblem(H,q0,p0,tspan;kwargs...)
  iip = (typeof(q0) <: AbstractArray) && !(typeof(q0) <: SArray)
  HamiltonianProblem{iip}(H,q0,p0,tspan;kwargs...)
end

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
        dq = function (t,x,v,dx)
            ForwardDiff.gradient!(dx, v-> H(x, v), v)
        end
        dp = function (t,x,v,dv)
            ForwardDiff.gradient!(dv, x->-H(x, v), x)
        end

        return ODEProblem{T}(DynamicalODEFunction{T}(dq,dp), (q0,p0),
                             tspan,HamiltonianProblem{true}(); kwargs...)
    end
end
