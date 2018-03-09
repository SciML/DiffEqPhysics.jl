module NBodyGravitational

using DiffEqBase, DifferentialEquations, StaticArrays, Plots, StaticArrays


grav_const = 6.673e-11

#=
    Represents a body/particle in an N-body gravitational problem
    m - mass
    r - initial position
    v - initial velocity
=#
struct MassBody{mType, cType}
    m :: mType
    r :: Vector{cType}
    v :: Vector{cType}
    
    function MassBody{mType,cType}(m::mType, r::Vector{cType}, v::Vector{cType})
        @assert length(r)==length(v) "The length of the position (r) and velocity (v) components should match"
        new(m,r,v)
    end
end

MassBody(m::mType, r :: Vector{cType}, v :: Vector{cType}) where {mType <: AbstractFloat, cType  <: AbstractFloat} = MassBody{mType, cType}(m,r,v)

#=
    Contains the necessary information for solving an N-body gravitational problem
    bodies - interacting bodies/particles
    tSpan - time, for which we want to evaluate the solution
    G - gravitational constant
=#
struct NBodyGravProblem
    bodies :: Array{MassBody, 1}
    tSpan :: Tuple{Float64, Float64}
    G :: Float64
end

NBodyGravProblem(bodies, tSpan::Tuple{Float64, Float64}) = NBodyGravProblem(bodies, tSpan, grav_const)

#= 
    As a quick approach to solving an N-body problem we will just 
    transform NBodyGravProblem to an ODEProblem

    Here u will be an array consisting of a consecutive positional coordinates r 
    following by components of velocities v
=#
import Base.convert
convert(::Type{ODEProblem}, x:: NBodyGravProblem) =
begin
    n = length(x.bodies)
    u0 = zeros(6*n);
    m = zeros(n)
    
    for i=1:n
        ith = 3(i-1)+1;        
        u0[ith:ith+2] = x.bodies[i].r;
        u0[3n+ith:3n+ith+2] = x.bodies[i].v 
        m[i] = x.bodies[i].m
    end    
    
    function ode_system!(du, u, p, t)
        du[1:3n] = @view u[3n+1:6n];
        for i=1:n
            ith = 3(i-1)+1;
            ri = @SVector [u[ith] , u[ith+1], u[ith+2]]
            du[3n+ith:3n+ith+2] = 
            begin
                accel = @SVector [0.0, 0.0, 0.0]; 
                for j=1:n
                    if j!=i
                        jth = 3(j-1)+1;
                        rj = @SVector [u[jth], u[jth+1], u[jth+2]]
                        accel -= x.G*m[j]*(ri-rj)/norm(ri-rj)^3
                    end
                end
                accel
            end
        end 
    end

    ODEProblem(ode_system!, u0, x.tSpan)
end

import DiffEqBase.solve
solve(x::NBodyGravProblem; kwargs...) = 
begin
    return DiffEqBase.solve(convert(ODEProblem,x); alg_hints=[:stiff], kwargs...)
end

# Create a gif file representing moving bodies and their orbits
function plot_xy_scattering(solution::ODESolution, path::AbstractString, duration::AbstractFloat = 3.0; plot_kwargs...)
    fps = 15
    n = Int(length(solution[1])/6)
    
    pl = plot()
    for i=1:n
        plot!(pl, solution, vars = (3(i-1)+1,3(i-1)+2); linecolor = :black, label = "Orbit $i", plot_kwargs...)
    end
    
    tmax = maximum(solution.t);
    tstep = tmax/(fps*duration);
    animation_data = solution(0.0:tstep:tmax)
    pos0 = animation_data[1];
    
    scatter!([pos0[3(i-1)+1] for i=1:n]', [pos0[3(i-1)+2] for i=1:n]'; markersize = 10, label = (["$i" for i=1:n]), plot_kwargs...)
    
    anim = @animate for i=1:length(animation_data)
        for j=1:n
           pl[j+n] = animation_data[i][3(j-1)+1], animation_data[i][3(j-1)+2]
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
    n = Int(length(solution[1])/6)
    tmax = maximum(solution.t);
    tstep = tmax/(fps*duration);
    nshift = Int(floor(ntrail/2));
    animation_data = solution(tstep : tstep : tmax-tstep)
    xy_data = zeros(length(animation_data), 2n)
    
    for i=1:n
        xy_data[:,i] = [d[3(i-1)+1] for d in  animation_data]
        xy_data[:,i+n] = [d[3(i-1)+2] for d in  animation_data]
    end
    
    xlim = 1.1*[minimum(minimum(xy_data[:, 1:n])), maximum(maximum(xy_data[:, 1:n]))]
    ylim = 1.1*[minimum(minimum(xy_data[:, n+1:2n])), maximum(maximum(xy_data[:, n+1:2n]))]
    anim = @animate for i=1+nshift:length(animation_data)-nshift
        plot(xy_data[i-nshift:i+nshift,1:n],xy_data[i-nshift:i+nshift,n+1:2n], xlim = xlim, ylim = ylim, lw=5)
    end
    gif(anim, path, fps = 15)
end

# A rather dull function for automatic plotting orbits of many bodies
function plot_xy(solution::ODESolution; kwargs...)
    n = Int(length(solution[1])/6)
    
    xlim = [Inf, -Inf]
    ylim = [Inf, -Inf]
    
    pl = plot()
    for i=1:n
        xindex = 3(i-1)+1;
        yindex = 3(i-1)+2;
        
        plot!(pl, solution, vars = (xindex, yindex); label = "Orbit $i", kwargs...)
        
        xn_min = minimum(solution[xindex, :]);
        xn_max = maximum(solution[xindex, :]);
        yn_min = minimum(solution[yindex, :]);
        yn_max = maximum(solution[yindex, :]);
        
        if xn_min < xlim[1] xlim[1] = xn_min end
        if xn_max > xlim[2] xlim[2] = xn_max end
        if yn_min < ylim[1] ylim[1] = yn_min end
        if yn_max > ylim[2] ylim[2] = yn_max end
    end
    plot!(pl, xlim = 1.1*xlim, ylim = 1.1*ylim)
    return pl
end

export NBodyGravProblem, MassBody, plot_xy, plot_xy_trailing, plot_xy_scattering

end