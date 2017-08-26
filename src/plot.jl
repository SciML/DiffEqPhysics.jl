using Plots

function plot_orbits(sol;body_names=nothing,dim=3,kwargs...)
    @assert dim ∈ (2, 3)
    N = length(sol.u[1].x[1]) ÷ dim
    body_names = body_names==nothing ? ["orbit $i" for i in 1:N] : body_names
    body_names = reshape(body_names, N, 1)
    ind = i -> (dim*i+1):(dim*(i+1))
    p = plot(sol, vars=(ind(0)...), lab=body_names[1], kwargs...)
    for i in 1:N-1
        plot!(p, sol, vars=(ind(i)...), lab=body_names[i+1])
    end
    p
end

