using DifferentialEquations, Plots

#Solving the simple pendulum with a traditional ODE method
#==========================================================#
function pendulum(du, u, params, t)
    #reference: http://www.pgccphy.net/ref/advmech.pdf page 6

    g = params[1] #gravitational acceleration
    m = params[2] #mass
    l = params[3] #length

    θ = u[1]
    ℒ = u[2]

    d_θ = ℒ / (m * l^2)
    d_ℒ = -m * g * l * sin(θ)

    du .= [d_θ, d_ℒ]
    return nothing
end

g = 9.81
m = 2.0
l = 1.0
params = [g, m, l]
u0 = [1.0, 1.0]

prob1 = ODEProblem(pendulum, u0, (0.0, 100.0), params)
sol1 = solve(prob1, AutoVern7(Rodas5()), dt = 0.05)

#==========================================================#

#Solving the simple pendulum with the DiffEqPhysics.jl HamiltonianProblem()
#==========================================================#
function H(ℒ, θ, params, t)
    g = params[1] #gravitational acceleration
    m = params[2] #mass
    l = params[3] #length

    return ℒ^2 / (2 * m * l^2) + m * g * l * (1 - cos(θ))
end

g = 9.81
m = 2.0
l = 1.0
params = [g, m, l]
θ₀ = 1.0
ℒ₀ = 1.0

prob2 = HamiltonianProblem(H, ℒ₀, θ₀, (0.0, 100.0), params)
sol2 = solve(prob2, SofSpa10(), dt = 0.05);
#==========================================================#

#Plotting
#==========================================================#
plot(sol1, vars = 1, tspan = (0, 20), label = "ODE Method")
plot!(sol2.t, sol2[2, :], xlim = (0, 20), label = "Hamiltonian Method")
#They produce the same solution!
