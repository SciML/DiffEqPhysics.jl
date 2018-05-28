include("./../src/nbody_simulation.jl")
println("====     Two particles interacting via the Lennard-Jones potential    ====")

m1 = 1.0
m2 = 1.0
r1 = 1.3

t1 = 0.0
t2 = 1.0
τ = (t2 - t1) / 1000

p1 = MassBody(SVector(-r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m1)
p2 = MassBody(SVector(r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m2)

σ = 1.0
ϵ = 1.0
parameters = LennardJonesParameters(ϵ, σ, Inf)
system = PotentialNBodySystem([p1, p2], Dict(:lennard_jones => parameters))
simulation = NBodySimulation(system, (t1, t2))
simResult = run_simulation(simulation, VelocityVerlet(), dt=τ)


r2 = get_position(simResult, t2, 2) - get_position(simResult, t2, 1)
vExpected = sqrt(4ϵ / m1 * ( ((σ / norm(r1))^12 - (σ / norm(r2))^12) - ((σ / norm(r1))^6 - (σ / norm(r2))^6 ) ))
vActual = norm(get_velocity(simResult, t2, 2))

ε = 0.001 * vExpected
@test vExpected ≈ vActual atol = ε