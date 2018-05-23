using DiffEqBase, OrdinaryDiffEq, StaticArrays, RecipesBase

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
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
    for j=1:n
        if j!=i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
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
    u0 = zeros(3, 2*n)
    m = zeros(n)
    
    for i=1:n
        u0[:, i] = x.bodies[i].r 
        u0[:, n+i] = x.bodies[i].v
        m[i] = x.bodies[i].m
    end    
    
    function ode_system!(du, u, p, t)
        du[:, 1:n] = @view u[:, n+1:2n];
        @inbounds for i=1:n
            du[:, n+i] .= gravitational_acceleration(x.G, u[:, 1:n], m, i, n);
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
    u0 = zeros(3, n);
    v0 = zeros(3, n);
    m = zeros(n)
    
    for i=1:n
        u0[:, i] = x.bodies[i].r;
        v0[:, i] = x.bodies[i].v 
        m[i] = x.bodies[i].m
    end    
    
    function gravitation!(dv,v,u,p,t)
        @inbounds for i=1:n            
            dv[:, i] .= gravitational_acceleration(x.G, u, m, i, n);
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

struct GravImagingData
    solution::ODESolution
end

@recipe function g(data::GravImagingData)
    solution = data.solution
    
    n = Int(size(solution[1],2)/2)
    
    xlim --> 1.1*[minimum(solution[1,1:n,:]), maximum(solution[1,1:n,:])]
    ylim --> 1.1*[minimum(solution[2,1:n,:]), maximum(solution[2,1:n,:])]        
    
    for i in 1:n
        @series begin
            label --> "Orbit $i"
            vars --> (3*(i-1)+1, 3*(i-1)+2)
            solution
        end
    end
end
