using OrdinaryDiffEq, ForwardDiff

function HamiltonianProblem(H,q0::Real,p0::Real,tspan::Tuple{T,T};kwargs...) where T<:Real
    @inline function dq(t,x,v)
        ForwardDiff.derivative(v->H(x, v), v)
    end
    @inline function dp(t,x,v)
        ForwardDiff.derivative(x->-H(x, v), x)
    end

    return ODEProblem((dq,dp), (q0,p0), tspan; kwargs...)
end

function HamiltonianProblem(H,q0::AbstractVector,p0::AbstractVector,tspan::Tuple{T,T};kwargs...) where T<:Real
    @assert length(q0) == length(p0)
    @inline function dq(t,x,v,dx)
        ForwardDiff.gradient!(dx, v-> H(x, v), v)
    end
    @inline function dp(t,x,v,dv)
        ForwardDiff.gradient!(dv, x->-H(x, v), x)
    end

    return ODEProblem((dq,dp), (q0,p0), tspan; kwargs...)
end

