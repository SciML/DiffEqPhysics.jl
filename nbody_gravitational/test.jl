println("====     NBodyGravProblem Test     ====")
println("====  Choreography cycle test case   ====")

include("NBodyGravitational.jl")

using Base.Test, NBodyGravitational, StaticArrays, OrdinaryDiffEq 
G = 1

# Testing the well-known figure Eight
m1 = MassBody(1.0, SVector(-0.995492, 0.0, 0.0), SVector(-0.347902, -0.53393, 0.0))
m2 = MassBody(1.0, SVector(0.995492, 0.0, 0.0), SVector(-0.347902, -0.53393, 0.0))
m3 = MassBody(1.0, SVector(0.0, 0.0, 0.0), SVector(0.695804, 1.067860, 0.0))
tspan = (0.0, 2pi);
problem = NBodyGravProblem([m1, m2, m3], G, tspan)
solution_simo_3 = solve(problem, Tsit5());
ε = 0.1
for j=1:3, i=1:3
    @test solution_simo_3[1][i,j] ≈ solution_simo_3[end][i,j] atol=ε
end


#the same NBodyGravProblem converted into SecondOrderODEProblem
solution_simo_3_2nd = solve(problem, DPRKN6(), transform_order=2);
ε = 0.001
for i=1:3, j=1:3
    @test solution_simo_3_2nd[1][9+3(i-1)+j] ≈ solution_simo_3_2nd[end][9+3(i-1)+j] atol=ε
end

#and using symplectic integrator VelocityVerlet
solution_simo_3_2nd = solve(problem, VelocityVerlet(), dt=pi/130, transform_order=2);
ε = 0.001
for i=1:3, j=1:3
    @test solution_simo_3_2nd[1][9+3(i-1)+j] ≈ solution_simo_3_2nd[end][9+3(i-1)+j] atol=ε
end

#using Yoshida6
solution_simo_3_2nd = solve(problem, Yoshida6(), dt=pi/12, transform_order=2);
ε = 0.001
for i=1:3, j=1:3
    @test solution_simo_3_2nd[1][9+3(i-1)+j] ≈ solution_simo_3_2nd[end][9+3(i-1)+j] atol=ε
end


# We need greater precision for a bigger choreography
m1 = MassBody(1.0, SVector(1.657666, 0.0, 0.0), SVector(0.0, -0.593786, 0.0))
m2 = MassBody(1.0, SVector(0.439775, -0.169717, 0.0), SVector(1.822785, 0.128248, 0.0))
m3 = MassBody(1.0, SVector(-1.268608, -0.267651, 0.0), SVector(1.271564, 0.168645, 0.0))
m4 = MassBody(1.0, SVector(-1.268608, 0.267651, 0.0), SVector(-1.271564, 0.168645, 0.0))
m5 = MassBody(1.0, SVector(0.439775, 0.169717, 0.0), SVector(-1.822785, 0.128248, 0.0))
tspan = (0.0, 2pi);
problem = NBodyGravProblem([m1, m2, m3, m4, m5], G, tspan)
solution_simo_5 = solve(problem, Tsit5(), abstol=1e-10, reltol=1e-10);

ε = 0.01
for j=1:5, i=1:3
    @test solution_simo_5[1][i,j] ≈ solution_simo_5[end][i,j] atol=ε
end