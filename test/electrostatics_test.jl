
include("./../src/nbody_base.jl")

println("====     Electrostatics Test     ====")

# small mass with negative charge rotating around more massive object with positive charge
r = 100.0
q1 = 1e-3
q2 = -1e-3
m1 = 100.0
m2 = 0.1
k = 9e9
v2 = sqrt(abs(k*q1*q2/m2/r))
t = 2*pi*r/v2
p1 = ChargedParticle(SVector(0.0, 0.0, 0.0), SVector(0.0, 0, 0.0), m1, q1)
p2 = ChargedParticle(SVector(r, 0.0, 0.0), SVector(0.0, v2, 0.0), m2, q2)
system = ChargedParticles([p1, p2], k)
simulation = NBodySimulation(system, SVector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0), (0.0, t))
result = run_simulation(simulation, Tsit5())

#temporary workaround
solution = result.solution;
ε = 0.1*r
for j=1:2, i=1:3
    @test solution[1][i,j] ≈ solution[end][i,j] atol=ε
end