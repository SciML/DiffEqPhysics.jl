println("====  Electrostatics Functional Tests  ====")
println("====  One particle rotates around another  ====")

# small mass with negative charge rotating around more massive object with positive charge
r = 100.0
q1 = 1e-3
q2 = -1e-3
m1 = 100.0
m2 = 0.1
k = 9e9
v2 = sqrt(abs(k * q1 * q2 / m2 / r))
t = 2 * pi * r / v2
p1 = ChargedParticle(SVector(0.0, 0.0, 0.0), SVector(0.0, 0, 0.0), m1, q1)
p2 = ChargedParticle(SVector(r, 0.0, 0.0), SVector(0.0, v2, 0.0), m2, q2)
system = ChargedParticles([p1, p2], k)
simulation = NBodySimulation(system, (0.0, t))
sim_result = run_simulation(simulation)


solution = sim_result.solution;
ε = 0.1 * r
for j = 1:2, i = 1:3
    @test solution[1][i,j] ≈ solution[end][i,j] atol = ε
end

println("====  Two positive charges repelling from each other  ====")
#   ("====               <---⊕-----r1-----⊕--->                ====")
q1 = 1e-3 # C
q2 = 1e-3 # C
m1 = 1.0 # kg
m2 = 1.0 # kg
r1 = 1 # m

E0 = k * q1 * q2 / r1 # initial energy
t1 = 0.0
t2 = 1.0
τ = (t2 - t1) / 1000

p1 = ChargedParticle(SVector(-r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m1, q1)
p2 = ChargedParticle(SVector(r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m2, q2)
system = ChargedParticles([p1, p2], k)
simulation = NBodySimulation(system, (t1, t2))
sim_result = run_simulation(simulation, VelocityVerlet(), dt=τ)

r2 = get_position(sim_result, t2, 2) - get_position(sim_result, t2, 1)
v_expected = sqrt(k * q1 * q2 / m1 * (1 / norm(r1) - 1 / norm(r2)))
v_actual = norm(get_velocity(sim_result, t2, 2))

ε = 0.001 * v_expected
@test v_expected ≈ v_actual atol = ε

for coordinates in sim_result
    @test length(coordinates) == 3
    for i=1:3
        @test length(coordinates[1]) == 2
    end
end


let default_potential = ElectrostaticParameters()
    @test 9e9 == default_potential.k
end

let
    io = IOBuffer()
    potential1 = ElectrostaticParameters()
    potential2 = ElectrostaticParameters(1.0)
    @test sprint(io -> show(io, potential1)) == "Electrostatic:\n\tk:9.0e9\n"
    @test sprint(io -> show(io, potential2)) == "Electrostatic:\n\tk:1.0\n"
end