include("../src/nbody_base.jl")

function generate_bodies(n::Int, m::AbstractFloat, T::AbstractFloat, L::AbstractFloat, kb::AbstractFloat)
    velocity_directions = generate_random_directions(n)
    mean_velocity  = sqrt(3kb*T/m )
    bodies = MassBody[]
    for i=1:n
        r = @SVector rand(3);
        v = mean_velocity*velocity_directions[i]
        body = MassBody(L*r, v, m)
        push!(bodies, body)
    end
    return bodies
end

function generate_random_directions(n::Int)
    theta = acos.(1-2*rand(n));
    phi = 2*pi*rand(n);
    directions = [@SVector [sin(theta[i]).*cos(phi[i]), sin(theta[i]).*sin(phi[i]), cos(theta[i])] for i=1:n]
end

T = 120.0 # °K
σ = 3.4e-10 # m
ρ = 1374 # kg/m^3
m = 39.95*1.6747*1e-24 # kg
L = 20σ #102.29σ 
N = floor(Int, ρ*L^3/m)
R = 2.25σ   
τ = 1e-14
kb= 1.38e-23
ϵ = T*kb
bodies = generate_bodies(N, m, T, L, kb)
ljSystem = PotentialNBodySystem(bodies, potentials = [:lennard_jones], lj_epsilon = ϵ, lj_sigma = σ, lj_range = R);
simulation = NBodySimulation(ljSystem, SVector(0.0, L, 0.0, L, 0.0, L), (0.0, 100τ));
result = run_simulation(simulation, Tsit5(), dt = τ)
