using OrdinaryDiffEq, ForwardDiff, RecursiveArrayTools

function NBodyProblem(potential,M,u0,v0,tspan; kwargs...)
    # check number of particles
    @assert length(u0) == length(v0)
    ind = i -> div(i-1, dim) + 1
    N   = length(u0)
    dim = div(N, length(M))
    @assert dim*length(M) == N
    @assert dim ∈ (2, 3)

    let fun = FWrapper{typeof(potential),typeof(tspan[1]),typeof(M)}(potential,tspan[1],M), cfg = ForwardDiff.GradientConfig(fun, u0)
        function acceleration!(t, x, v, dv)
            fun.t = t
            ForwardDiff.gradient!(dv, fun, x, cfg)
            for x in dv.x
                x ./= -M
            end
        end
        return SecondOrderODEProblem{true}(acceleration!, u0, v0, tspan; kwargs...)
    end
end

mutable struct FWrapper{F,T,MType}
    f::F
    t::T
    M::MType
end
(f::FWrapper)(u) = f.f(f.t,u.x...,f.M)
