struct OrbitPlot
    sol
    body_names
    dim
end
@recipe function f(p::OrbitPlot)
    sol = p.sol
    body_names = p.body_names
    dim = p.dim
    @assert dim ∈ (2, 3)
    N = length(sol.u[1].x[1].x[1])
    body_names = body_names==nothing ? ["orbit $i" for i in 1:N] : body_names
    ind = i -> i:N:length(sol.u[1].x[1])
    for i in 1:N
        @series begin
            vars=(ind(i)...,)
            label --> body_names[i]
            vars --> vars
            sol
        end
    end
end

function orbitplot(sol::DESolution;body_names=nothing,dim=3,kwargs...)
    RecipesBase.plot(OrbitPlot(sol,body_names,dim);kwargs...)
end

export orbitplot

function plot_orbits(sol;body_names=nothing,dim=3,kwargs...)
    @assert dim ∈ (2, 3)
    N = length(sol.u[1].x[1].x[1])
    body_names = body_names==nothing ? ["orbit $i" for i in 1:N] : body_names
    ind = i -> i:N:length(sol.u[1].x[1])
    p = plot(sol, vars=(ind(1)...,), lab=body_names[1], kwargs...)
    for i in 2:N
        plot!(p, sol, vars=(ind(i)...,), lab=body_names[i], kwargs...)
    end
    p
end
