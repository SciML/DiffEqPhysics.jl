using Plots, DiffEqBase

# Create a gif file representing moving bodies and their orbits
function plot_xy_scattering(solution::ODESolution, path::AbstractString, duration::AbstractFloat = 3.0; plot_kwargs...)
    fps = 15
    n = Int(size(solution[1], 2)/2)
    
    xlim = [minimum(solution[1,1:n,:]), maximum(solution[1,1:n,:])]  
    ylim = [minimum(solution[2,1:n,:]), maximum(solution[2,1:n,:])]  
    
    
    pl = plot()
    for i=1:n
        plot!(pl, solution, vars=(3*(i-1)+1, 3*(i-1)+2); linecolor = :black, label = "Orbit $i", xlim = 1.1*xlim, ylim = 1.1*ylim, plot_kwargs...)
    end
    
    tmax = maximum(solution.t);
    tstep = tmax/(fps*duration);
    animation_data = solution(0.0:tstep:tmax)
    pos0 = animation_data[1];
    
    scatter!([pos0[1, i] for i=1:n]', [pos0[2, i] for i=1:n]'; markersize = 10, label = (["$i" for i=1:n]), plot_kwargs...)
    
    anim = @animate for i=1:length(animation_data)
        for j=1:n
           pl[j+n] = animation_data[i][1,j], animation_data[i][2,j]
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
    
    xlim = [minimum(xy_data[1,1:n,:]), maximum(xy_data[1,1:n,:])]
    ylim = [minimum(xy_data[2,1:n,:]), maximum(xy_data[2,1:n,:])]
    
    anim = @animate for i=1+nshift:size(xy_data, 3)-nshift
        plot(xy_data[1,1:n,i-nshift:i+nshift]',xy_data[2,1:n,i-nshift:i+nshift]', xlim = xlim, ylim = ylim, lw=5)
    end
    gif(anim, path, fps = 15)
end

# A rather dull function for automatic plotting orbits of many bodies
function plot_xy(solution::ODESolution; kwargs...)
    n = size(solution[1],1)
    
    xlim = [minimum(solution[1,1:n,:]), maximum(solution[1,1:n,:])]  
    ylim = [minimum(solution[2,1:n,:]), maximum(solution[2,1:n,:])] 
    
    pl = plot()
    for i=1:n
        plot!(pl, solution, vars=(3*(i-1)+1, 3*(i-1)+2); label = "Orbit $i", kwargs...)
    end
    plot!(pl, xlim = 1.1*xlim, ylim = 1.1*ylim)
    return pl
end