using OrdinaryDiffEq, ForwardDiff

function NBodyProblem(potential,mass_matrix,u0,v0,tspan; kwargs...)
    # check number of particles
    @assert length(u0) == length(v0)
    @assert isdiag(mass_matrix)
    ind = i -> div(i-1, dim) + 1
    N   = length(u0)
    dim = length(u0.x)
    @assert dim*size(mass_matrix, 1) == N
    @assert dim ∈ (2, 3)

    f_wrapper(t,u) = potential(t,u.x...,mass_matrix)

    function acceleration!(t, x, v, dv)
        N = length(v)
        config_length = N > 10 ? 10 : N
        fun = q -> f_wrapper(t, q)
        cfg = ForwardDiff.GradientConfig(fun, u0, ForwardDiff.Chunk{config_length}())
        ForwardDiff.gradient!(dv, fun, x, cfg)
        for i in eachindex(dv)
            dv[i] /= -mass_matrix.diag[ind(i)]
        end
    end
    SecondOrderODEProblem{true}(acceleration!, u0, v0, tspan; kwargs...)
end

