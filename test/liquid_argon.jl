#include("../src/nbody_simulation.jl")

function generate_bodies_randomly(n::Int, m::AbstractFloat, mean_velocity::AbstractFloat, L::AbstractFloat)
    velocity_directions = generate_random_directions(n)
    velocities =  sqrt(mean_velocity) * velocity_directions
    bodies = MassBody[]
    for i = 1:n
        r = @SVector rand(3);
        v = velocities[i]
        body = MassBody(L * r, v, m)
        push!(bodies, body)
    end
    return bodies
end

function generate_bodies_in_cell_nodes(n::Int, m::AbstractFloat, v_dev::AbstractFloat, L::AbstractFloat)

    dL = L / n^(1 / 3)
    rng = MersenneTwister(n);
    velocities = v_dev * randn(rng, Float64, (3, n))
    bodies = MassBody[]

    count = 1
    for x = 0:dL:L, y = 0:dL:L, z = 0:dL:L        
        if count > n
            break
        end
        r = SVector(x, y, z)
        v = SVector{3}(velocities[:,count])
        body = MassBody(r, v, m)
        push!(bodies, body)
        count += 1           
    end
    return bodies
end

function generate_random_directions(n::Int)
    theta = acos.(1 - 2 * rand(n));
    phi = 2 * pi * rand(n);
    directions = [@SVector [sin(theta[i]) .* cos(phi[i]), sin(theta[i]) .* sin(phi[i]), cos(theta[i])] for i = 1:n]
end

T = 90.0 # °K
kb = 1.38e-23 # J/K
ϵ = 120 * kb
σ = 3.4e-10 # m
ρ = 1374 # kg/m^3
m = 39.95 * 1.6747 * 1e-27 # kg
L = 10.229σ 
N = floor(Int, ρ * L^3 / m)
R = 2.25σ   
v_dev = sqrt(kb * T / m)
τ = 1e-14 # σ/v
t1 = 0.0
t2 = 300τ
#bodies = generate_bodies_randomly(N, m, v_dev, L)
bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)
parameters = LennardJonesParameters(ϵ, σ, R)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
simulation = NBodySimulation(lj_system, (t1, t2), PeriodicBoundaryConditions(L));
#result = run_simulation(simulation, Tsit5())
result = run_simulation(simulation, VelocityVerlet(), dt=τ)