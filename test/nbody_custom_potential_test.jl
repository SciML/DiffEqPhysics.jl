include("./../src/nbody_simulation.jl")

m1 = 1.0
m2 = 1.0

t1 = 0.0
t2 = 1.0
τ = (t2 - t1) / 100

p1 = MassBody(SVector(1, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m1)
p2 = MassBody(SVector(-1, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m2)


struct CustomPotentialParameters <: PotentialParameters
    a::AbstractFloat
end

function get_accelerating_function(p::CustomPotentialParameters, simulation::NBodySimulation)
    (dv, u, v, t, i) -> begin custom_accel = SVector(p.a, 0.0, 0.0); dv .= custom_accel end 
end

parameters = CustomPotentialParameters(1.5)
system = PotentialNBodySystem([p1, p2], Dict(:custom_potential_params => parameters))
simulation = NBodySimulation(system, (t1, t2))
simResult = run_simulation(simulation, VelocityVerlet(), dt=τ)

v2 = get_velocity(simResult, t2, 1)
r2 = get_position(simResult, t2, 1)

ε = 1e-6
@test 1.5 ≈ v2[1] atol = ε
@test 1.75 ≈ r2[1] atol = ε