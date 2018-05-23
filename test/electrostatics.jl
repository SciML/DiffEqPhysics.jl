include("../src/nbody_base.jl")


p1 = ChargedParticle(SVector(-10.0, 0.0, 0.0), SVector(0.0, 10.0, 0.0), 1.0, 1e-3)
p2 = ChargedParticle(SVector(10.0, 0.0, 0.0), SVector(0.0, -10.0, 0.0), 1.0, -4e-3)
system = ChargedParticles([p1, p2], 9e9)
simulation = NBodySimulation(system, SVector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), (0.0, 1.0))
result = run_simulation(simulation, Tsit5())