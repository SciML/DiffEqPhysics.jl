include("../src/nbody_simulation.jl")

function generate_bodies(n::Int, m::AbstractFloat, mean_velocity::AbstractFloat, L::AbstractFloat)
    velocity_directions = generate_random_directions(n)
    bodies = MassBody[]
    for i = 1:n
        r = @SVector rand(3);
        v = mean_velocity * velocity_directions[i]
        body = MassBody(L * r, v, m)
        push!(bodies, body)
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
v = sqrt(3 * kb * T / m)
τ = 1e-14 # σ/v
t1 = 0.0
t2 = 800τ
bodies = generate_bodies(N, m, v, L)
parameters = LennardJonesParameters(ϵ, σ, R)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
simulation = NBodySimulation(lj_system, (t1, t2), PeriodicBoundaryConditions(L));
#result = run_simulation(simulation, Tsit5())
result = run_simulation(simulation, VelocityVerlet(), dt=τ)


using Plots

vs = Float64[]
for i = 1:N
    v = get_velocity(result, t2,i)
    push!(vs,norm(v))
end
histogram(vs, bins=0:11, xlim=[0, maximum(vs)], ylim=[0,1])
