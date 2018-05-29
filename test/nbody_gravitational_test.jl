using DiffEqPhysics, Base.Test, StaticArrays, OrdinaryDiffEq 

println("====  Gravitational Functional Test  ====")
println("====  Choreography cycle test case   ====")

G = 1

# Testing the well-known figure Eight
m1 = MassBody(SVector(-0.995492, 0.0, 0.0), SVector(-0.347902, -0.53393, 0.0), 1.0)
m2 = MassBody(SVector(0.995492, 0.0, 0.0), SVector(-0.347902, -0.53393, 0.0), 1.0)
m3 = MassBody(SVector(0.0, 0.0, 0.0), SVector(0.695804, 1.067860, 0.0), 1.0)
tspan = (0.0, 2pi);
system = GravitationalSystem([m1, m2, m3], G)
simulation = NBodySimulation(system, tspan)
sim_result = run_simulation(simulation)
solution_simo_3 = sim_result.solution;
ε = 0.1
for j = 1:3, i = 1:3
    @test solution_simo_3[1][i,j] ≈ solution_simo_3[end][i,j] atol = ε
end


#the same NBodyGravProblem converted into SecondOrderODEProblem
sim_result = run_simulation(simulation, DPRKN6())
solution_simo_3_2nd = sim_result.solution;
ε = 0.001
for i = 1:3, j = 1:3
    @test solution_simo_3_2nd[1][9 + 3(i - 1) + j] ≈ solution_simo_3_2nd[end][9 + 3(i - 1) + j] atol = ε
end

#and using symplectic integrator VelocityVerlet
sim_result = run_simulation(simulation, VelocityVerlet(), dt=pi / 130)
solution_simo_3_2nd = sim_result.solution;
ε = 0.001
for i = 1:3, j = 1:3
    @test solution_simo_3_2nd[1][9 + 3(i - 1) + j] ≈ solution_simo_3_2nd[end][9 + 3(i - 1) + j] atol = ε
end

#using Yoshida6
sim_result = run_simulation(simulation, Yoshida6(), dt=pi / 12)
solution_simo_3_2nd = sim_result.solution;
ε = 0.001
for i = 1:3, j = 1:3
    @test solution_simo_3_2nd[1][9 + 3(i - 1) + j] ≈ solution_simo_3_2nd[end][9 + 3(i - 1) + j] atol = ε
end


# We need greater precision for a bigger choreography
m1 = MassBody(SVector(1.657666, 0.0, 0.0), SVector(0.0, -0.593786, 0.0), 1.0)
m2 = MassBody(SVector(0.439775, -0.169717, 0.0), SVector(1.822785, 0.128248, 0.0), 1.0)
m3 = MassBody(SVector(-1.268608, -0.267651, 0.0), SVector(1.271564, 0.168645, 0.0), 1.0)
m4 = MassBody(SVector(-1.268608, 0.267651, 0.0), SVector(-1.271564, 0.168645, 0.0), 1.0)
m5 = MassBody(SVector(0.439775, 0.169717, 0.0), SVector(-1.822785, 0.128248, 0.0), 1.0)
tspan = (0.0, 2pi);
system = GravitationalSystem([m1, m2, m3, m4, m5], G)
simulation = NBodySimulation(system, tspan)
sim_result = run_simulation(simulation, Tsit5(), abstol=1e-10, reltol=1e-10)
solution_simo_5 = sim_result.solution;

ε = 0.01
for j = 1:5, i = 1:3
    @test solution_simo_5[1][i,j] ≈ solution_simo_5[end][i,j] atol = ε
end