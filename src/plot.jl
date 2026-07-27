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
    orbitplot(sol; body_names=nothing, dim=3)

Create a recipe-based orbit-plot input for an N-body-style solution.

# Arguments

- `sol`: SciML solution whose saved states store position variables in the
  N-body layout used by DiffEqPhysics-style orbit examples.

# Keyword Arguments

- `body_names`: Optional names for the plotted bodies. When `nothing`, labels are
  generated as `"orbit 1"`, `"orbit 2"`, and so on.
- `dim`: Plot dimension. Must be `2` or `3`.

# Returns

- An object that a plotting backend supporting RecipesBase recipes can render.

# Examples

```julia
using Plots

plot(orbitplot(sol; body_names = ["Sun", "Planet"], dim = 2))
```
"""
function orbitplot(sol::AbstractSciMLSolution; body_names = nothing, dim = 3)
    @assert dim ∈ (2, 3)
    return OrbitPlot(sol, body_names, dim)
end

export orbitplot

"""
    plot_orbits(sol; body_names=nothing, dim=3)

Create a recipe-based orbit-plot input for an N-body-style solution.

# Arguments

- `sol`: SciML solution whose saved states store body coordinates in the layout
  expected by the DiffEqPhysics orbit plotting recipe.

# Keyword Arguments

- `body_names`: Optional names for the plotted bodies. When `nothing`, labels are
  generated as `"orbit 1"`, `"orbit 2"`, and so on.
- `dim`: Plot dimension. Must be `2` or `3`.

# Returns

- An object that a plotting backend supporting RecipesBase recipes can render.

# Examples

```julia
using Plots

plot(plot_orbits(sol; body_names = ["Sun", "Planet"], dim = 2))
```

See also: [`orbitplot`](@ref).
"""
function plot_orbits(sol; body_names = nothing, dim = 3)
    return orbitplot(sol; body_names, dim)
end
