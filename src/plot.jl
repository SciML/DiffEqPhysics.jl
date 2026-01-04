"""
Internal type for orbit plotting recipe.
"""
struct OrbitPlot
    sol::Any
    body_names::Any
    dim::Any
end
@recipe function f(p::OrbitPlot)
    sol = p.sol
    body_names = p.body_names
    dim = p.dim
    @assert dim ∈ (2, 3)
    N = length(sol.u[1].x[1].x[1])
    body_names = body_names === nothing ? ["orbit $i" for i in 1:N] : body_names
    ind = i -> i:N:length(sol.u[1].x[1])
    for i in 1:N
        @series begin
            vars = (ind(i)...,)
            label --> body_names[i]
            vars --> vars
            sol
        end
    end
end

"""
    orbitplot(sol; body_names=nothing, dim=3, kwargs...)

Plot the orbits from an N-body simulation solution using Plots.jl recipes.

# Arguments
- `sol`: A `DESolution` from solving an N-body problem
- `body_names`: Optional vector of names for each body (defaults to "orbit 1", "orbit 2", etc.)
- `dim`: Dimension of the plot, either `2` or `3` (default: `3`)
- `kwargs...`: Additional keyword arguments passed to the plot recipe

# Returns
A plot object showing the trajectories of all bodies.
"""
function orbitplot(sol::DESolution; body_names = nothing, dim = 3, kwargs...)
    return RecipesBase.plot(OrbitPlot(sol, body_names, dim); kwargs...)
end

export orbitplot

"""
    plot_orbits(sol; body_names=nothing, dim=3, kwargs...)

Plot the orbits from an N-body simulation solution.

# Arguments
- `sol`: A solution from solving an N-body problem
- `body_names`: Optional vector of names for each body (defaults to "orbit 1", "orbit 2", etc.)
- `dim`: Dimension of the plot, either `2` or `3` (default: `3`)
- `kwargs...`: Additional keyword arguments passed to `plot`

# Returns
A plot object showing the trajectories of all bodies.

See also: [`orbitplot`](@ref)
"""
function plot_orbits(sol; body_names = nothing, dim = 3, kwargs...)
    @assert dim ∈ (2, 3)
    N = length(sol.u[1].x[1].x[1])
    body_names = body_names === nothing ? ["orbit $i" for i in 1:N] : body_names
    ind = i -> i:N:length(sol.u[1].x[1])
    p = plot(sol, vars = (ind(1)...,), lab = body_names[1], kwargs...)
    for i in 2:N
        plot!(p, sol, vars = (ind(i)...,), lab = body_names[i], kwargs...)
    end
    return p
end
