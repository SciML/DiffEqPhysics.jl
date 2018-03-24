module NBodyGravitational

using DiffEqBase, OrdinaryDiffEq, StaticArrays, Plots

#=
    Represents a body/particle in an N-body gravitational problem
    m - mass
    r - initial position
    v - initial velocity
=#
struct MassBody{mType, cType <: AbstractFloat}
    m :: mType
    r :: SVector{3, cType}
    v :: SVector{3, cType}
end

#=
    Contains the necessary information for solving an N-body gravitational problem
    bodies - interacting bodies/particles
    tspan - time, for which we want to evaluate the solution
    G - gravitational constant
=#
struct NBodyGravProblem
    bodies :: Vector{MassBody}
    G :: Float64
    tspan :: Tuple{Float64, Float64}
end


#should be optimized using the 3d Newton's law
function gravitational_acceleration(G::Float64, rs, ms::Vector{mType}, i::Integer, n::Integer) where {mType<:AbstractFloat}
    accel = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[i, 1], rs[i, 2], rs[i, 3]]
    for j=1:n
        if j!=i
            rj = @SVector [rs[j, 1], rs[j, 2], rs[j, 3]]
            accel -= G*ms[j]*(ri-rj)/norm(ri-rj)^3
        end
    end
    
    return accel
end


#= 
    As a quick approach to solving an N-body problem we will just 
    transform NBodyGravProblem to ODEProblem
    Here u will be an array consisting of a consecutive positional coordinates r 
    following by components of velocities v
=#
function DiffEqBase.ODEProblem(x::NBodyGravProblem)
    n = length(x.bodies)
    u0 = zeros(n, 6);
    m = zeros(n)
    
    for i=1:n
        u0[i, 1:3] = x.bodies[i].r;
        u0[i, 4:6] = x.bodies[i].v 
        m[i] = x.bodies[i].m
    end    
    
    function ode_system!(du, u, p, t)
        du[:, 1:3] = @view u[:, 4:6];
        for i=1:n
            du[i, 4:6] .= gravitational_acceleration(x.G, u[:, 1:3], m, i, n);
        end 
    end

    ODEProblem(ode_system!, u0, x.tspan)
end

#= 
    Since the second Newton's law is an ODE of the 2nd order,
    we can also transform NBodyGravProblem into SecondOrderODEProblem.
=#
function DiffEqBase.SecondOrderODEProblem(x::NBodyGravProblem)
    n = length(x.bodies)
    u0 = zeros(n, 3);
    v0 = zeros(n, 3);
    m = zeros(n)
    
    for i=1:n
        u0[i, 1:3] = x.bodies[i].r;
        v0[i, 1:3] = x.bodies[i].v 
        m[i] = x.bodies[i].m
    end    
    
    function gravitation!(dv,v,u,p,t)
        for i=1:n            
            dv[i, 1:3] .= gravitational_acceleration(x.G, u[:, 1:3], m, i, n);
        end 
    end

    SecondOrderODEProblem(gravitation!, v0, u0, x.tspan)
end

import DiffEqBase: solve
function solve(x::NBodyGravProblem, args...; transform_order=1, kwargs...)
    if transform_order == 1
        solve(ODEProblem(x), args...; kwargs...)
    elseif transform_order == 2
        solve(SecondOrderODEProblem(x), args...; kwargs...)        
    else
        throw(ArgumentError("unsupported transformation order of the problem (choose 1 or 2)."))
    end
end



# Create a gif file representing moving bodies and their orbits
function plot_xy_scattering(solution::ODESolution, path::AbstractString, duration::AbstractFloat = 3.0; plot_kwargs...)
    fps = 15
    n = size(solution[1], 1)
    
    xlim = [minimum(solution[:,1,:]), maximum(solution[:,1,:])]
    ylim = [minimum(solution[:,2,:]), maximum(solution[:,2,:])]
    
    
    pl = plot()
    for i=1:n
        plot!(pl, solution, vars=(i,i+n); linecolor = :black, label = "Orbit $i", xlim = 1.1*xlim, ylim = 1.1*ylim, plot_kwargs...)
    end
    
    tmax = maximum(solution.t);
    tstep = tmax/(fps*duration);
    animation_data = solution(0.0:tstep:tmax)
    pos0 = animation_data[1];
    
    scatter!([pos0[i, 1] for i=1:n]', [pos0[i, 2] for i=1:n]'; markersize = 10, label = (["$i" for i=1:n]), plot_kwargs...)
    
    anim = @animate for i=1:length(animation_data)
        for j=1:n
           pl[j+n] = animation_data[i][j,1], animation_data[i][j,2]
        end
    end
    gif(anim, path, fps = fps) 
end

#= 
    Creates a gif file representing moving bodies with a trace/tail of their movement
    consisting of `ntrail` points of their orbits
=#
function plot_xy_trailing(solution::ODESolution, path::AbstractString; ntrail::Int = 3, duration::AbstractFloat = 3.0)
    # ntrail - number of points for displaying a trailing path; the trail will correspond to the velocity of a body
    fps = 15
    n = size(solution[1], 1)
    tmax = maximum(solution.t);
    tstep = tmax/(fps*duration);
    nshift = Int(floor(ntrail/2));
    xy_data = solution(tstep : tstep : tmax-tstep)
    
    xlim = [minimum(xy_data[:,1,:]), maximum(xy_data[:,1,:])]
    ylim = [minimum(xy_data[:,2,:]), maximum(xy_data[:,2,:])]
    
    anim = @animate for i=1+nshift:size(xy_data, 3)-nshift
        plot(xy_data[1:n,1,i-nshift:i+nshift]',xy_data[1:n,2,i-nshift:i+nshift]', xlim = xlim, ylim = ylim, lw=5)
    end
    gif(anim, path, fps = 15)
end

# A rather dull function for automatic plotting orbits of many bodies
function plot_xy(solution::ODESolution; kwargs...)
    n = size(solution[1],1)
    
    xlim = [minimum(solution[:,1,:]), maximum(solution[:,1,:])]
    ylim = [minimum(solution[:,2,:]), maximum(solution[:,2,:])]
    
    pl = plot()
    for i=1:n
        plot!(pl, solution, vars=(i, i+n); label = "Orbit $i", kwargs...)
    end
    plot!(pl, xlim = 1.1*xlim, ylim = 1.1*ylim)
    return pl
end

export NBodyGravProblem, MassBody, plot_xy, plot_xy_trailing, plot_xy_scattering, solve

end