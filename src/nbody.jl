using OrdinaryDiffEq, ForwardDiff

function NBodyProblem(potential,mass_matrix,u0,v0,tspan; kwargs...)
    # check number of particles
    @assert length(u0) == length(v0)
    @assert isdiag(mass_matrix)
    ind = i -> div(i-1, dim) + 1
    N   = length(u0)
    dim = div(N, size(mass_matrix, 1))
    @assert dim*size(mass_matrix, 1) == N
    @assert dim ∈ (2, 3)

    function f_wrapper2(t,u)
        x = @view u[1:2:end]
        y = @view u[2:2:end]
        potential(t,x,y,mass_matrix)
    end

    function f_wrapper3(t,u)
        x = @view u[1:3:end]
        y = @view u[2:3:end]
        z = @view u[3:3:end]
        potential(t,x,y,z,mass_matrix)
    end

    function acceleration!(t, x, v, dv)
        N = length(v)
        config_length = N > 10 ? 10 : N
        fun = dim == 2 ? q -> f_wrapper2(t, q) : q -> f_wrapper3(t, q)
        cfg = ForwardDiff.GradientConfig(fun, u0, ForwardDiff.Chunk{config_length}())
        ForwardDiff.gradient!(dv, fun, x, cfg)
        for i in eachindex(dv)
            dv[i] /= -mass_matrix.diag[ind(i)]
        end
    end
    SecondOrderODEProblem{true}(acceleration!, u0, v0, tspan; kwargs...)
end
