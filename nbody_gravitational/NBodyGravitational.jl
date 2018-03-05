module NBodyGravitational

using DiffEqBase, DifferentialEquations, StaticArrays, Plots, StaticArrays


grav_const = 6.673e-11

struct MassBody{mType, cType}
    m :: mType
    r :: SVector{3, cType}
    v :: SVector{3, cType}
    
    function MassBody{mType,cType}(m::mType,r::AbstractArray{cType,1}, v::AbstractArray{cType,1})
        @assert length(r)==length(v) "The length of the position (r) and velocity (v) components should match"
        new(m,r,v)
    end
end

MassBody(m::mType, r :: SVector{3, cType}, v :: SVector{3, cType}) where {mType <: AbstractFloat, cType  <: AbstractFloat} = MassBody{mType, cType}(m,r,v)

struct NBodyGravProblem
    bodies :: AbstractArray{MassBody, 1}
    tSpan :: Tuple{Float64, Float64}
    G :: Float64
end

NBodyGravProblem(bodies, tSpan::Tuple{Float64, Float64}) = NBodyGravProblem(bodies, tSpan, grav_const)

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
    
    function grav_acceleration(m2, r1, r2)
        return -x.G*m2*(r1-r2)/norm(r1-r2)^3;
    end
    
    
    function ode_system!(du, u, p, t)
        for i=1:n
            ith = 3(i-1)+1;        
            du[ith:ith+2] = u[3n+ith:3n+ith+2];
            du[3n+ith:3n+ith+2] = 
            begin
                accel=zeros(3); 
                j=1; 
                while j<i 
                    jth = 3(j-1)+1;
                    accel += grav_acceleration(m[j],u[ith:ith+2],u[jth:jth+2])
                    j+=1
                end
                j+=1;
                while j<=n 
                    jth = 3(j-1)+1;
                    accel += grav_acceleration(m[j],u[ith:ith+2],u[jth:jth+2])
                    j+=1
                end
                accel
            end
        end 
    end

    ODEProblem(ode_system!, u0, x.tSpan)
end

import DiffEqBase.solve
solve(x::NBodyGravProblem) = 
begin
    return DiffEqBase.solve(convert(ODEProblem,x), alg_hints=[:stiff])
end

function plot_xy_scattering(solution::ODESolution, path::AbstractString, duration::AbstractFloat = 3.0)
    fps = 15
    n = Int(length(solution[1])/6)
    
    pl = plot()
    for i=1:n
        plot!(pl, solution, vars = (3(i-1)+1,3(i-1)+2), linecolor = :black, label = "Orbit $i")
    end
    
    tmax = maximum(solution.t);
    tstep = tmax/(fps*duration);
    animation_data = solution(0.0:tstep:tmax)
    pos0 = animation_data[1];
    
     scatter!([pos0[3(i-1)+1] for i=1:n]', [pos0[3(i-1)+2] for i=1:n]', markersize = 10, label = (["$i" for i=1:n]))
    
    anim = @animate for i=1:length(animation_data)
        for j=1:n
           pl[j+n] = animation_data[i][3(j-1)+1], animation_data[i][3(j-1)+2]
        end
    end
    gif(anim, path, fps = fps) 
end


function plot_xy_trailing(solution::ODESolution, path::AbstractString; ntrail::Int = 3, duration::AbstractFloat = 3.0)
    # ntrail - number of points for displaing a trailing path; the trail will correspond to the velocity of a body
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

export NBodyGravProblem, MassBody, plot_xy_trailing, plot_xy_scattering

end