function NBodyProblem(potential,M,v0,u0,tspan,p=nothing; kwargs...)
    # check number of particles
    @assert length(u0) == length(v0)
    ind = i -> div(i-1, dim) + 1
    N   = length(u0)
    dim = div(N, length(M))
    @assert dim*length(M) == N
    @assert dim ∈ (2, 3)

    let fun = FWrapper{typeof(potential),typeof(tspan[1]),typeof(M),typeof(p)}(potential,tspan[1],M,p), cfg = ForwardDiff.GradientConfig(fun, u0)
        function acceleration!(dv,v,x,p,t)
            fun.t = t
            ForwardDiff.gradient!(dv, fun, x, cfg)
            for x in dv.x
                x ./= -M
            end
        end
        return SecondOrderODEProblem{true}(acceleration!, v0, u0, tspan, p; kwargs...)
    end
end

mutable struct FWrapper{F,T,MType,P}
    f::F
    t::T
    M::MType
    p::P
end
(f::FWrapper)(u) = f.f(f.p,f.t,u.x...,f.M)
