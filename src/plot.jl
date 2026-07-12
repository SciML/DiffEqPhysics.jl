# Internal type for the orbit plotting recipe.
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

Return a recipe-based orbit plot for an N-body-style solution.

# Arguments

- `sol`: SciML solution whose saved states store position variables in the
  N-body layout used by DiffEqPhysics-style orbit examples.

# Keyword Arguments

- `body_names`: Optional names for the plotted bodies. When `nothing`, labels are
  generated as `"orbit 1"`, `"orbit 2"`, and so on.
- `dim`: Plot dimension. Must be `2` or `3`.
- `kwargs...`: Additional keyword arguments forwarded to `RecipesBase.plot`.

# Returns

- A plot object returned by the active plotting backend.

# Examples

```julia
orbitplot(sol; body_names = ["Sun", "Planet"], dim = 2)
```
"""
function orbitplot(sol::AbstractSciMLSolution; body_names = nothing, dim = 3, kwargs...)
    return RecipesBase.plot(OrbitPlot(sol, body_names, dim); kwargs...)
end

export orbitplot

"""
    plot_orbits(sol; body_names=nothing, dim=3, kwargs...)

Plot the orbit of each body in an N-body-style solution.

# Arguments

- `sol`: SciML solution whose saved states store body coordinates in the layout
  expected by the DiffEqPhysics orbit plotting recipe.

# Keyword Arguments

- `body_names`: Optional names for the plotted bodies. When `nothing`, labels are
  generated as `"orbit 1"`, `"orbit 2"`, and so on.
- `dim`: Plot dimension. Must be `2` or `3`.
- `kwargs...`: Additional keyword arguments forwarded to `plot`.

# Returns

- A plot object with one trajectory series per body.

# Examples

```julia
plot_orbits(sol; body_names = ["Sun", "Planet"], dim = 2)
```

See also: [`orbitplot`](@ref).
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
