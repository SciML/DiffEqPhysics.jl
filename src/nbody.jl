using OrdinaryDiffEq, ForwardDiff, RecursiveArrayTools

function NBodyProblem(potential,mass,u0,v0,tspan; kwargs...)
    # check number of particles
    @assert length(u0) == length(v0)
    ind = i -> div(i-1, dim) + 1
    N   = length(u0)
    dim = div(N, length(mass))
    @assert dim*length(mass) == N
    @assert dim ∈ (2, 3)

    if dim == 3
        divd = N÷3
        u0 = ArrayPartition(u0[1:3:end], u0[2:3:end],        u0[3:3:end])
        v0 = ArrayPartition(v0[1:divd],  v0[(divd+1):2divd], v0[(2divd+1):end])
    else
        divd = N÷2
        u0 = ArrayPartition(u0[1:2:end], u0[2:2:end])
        v0 = ArrayPartition(v0[1:divd],  v0[(divd+1):end])
    end

    f_wrapper(t,u) = potential(t,u.x...,mass)

    function acceleration!(t, x, v, dv)
        N = length(v)
        config_length = N > 10 ? 10 : N
        fun = q -> f_wrapper(t, q)
        cfg = ForwardDiff.GradientConfig(fun, u0, ForwardDiff.Chunk{config_length}())
        ForwardDiff.gradient!(dv, fun, x, cfg)
        for i in eachindex(dv)
            dv[i] /= -mass[ind(i)]
        end
    end
    SecondOrderODEProblem{true}(acceleration!, u0, v0, tspan; kwargs...)
end

