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

    d_θ = ℒ/(m*l^2)
    d_ℒ = -m*g*l*sin(θ)

    du .= [d_θ, d_ℒ]
    return nothing
end

g = 9.81
m = 2.
l = 1.
params = [g,m,l]
u0 = [1.0,1.0]

prob = ODEProblem(pendulum, u0, (0., 100.), params)
sol1 = solve(prob, AutoVern7(Rodas5()), dt = .05)

#==========================================================#



#Solving the simple pendulum with the DiffEqPhysics.jl HamiltonianProblem()
#==========================================================#
function H(θ, ℒ, params)
    g = params[1] #gravitational acceleration
    m = params[2] #mass
    l = params[3] #length

    return ℒ^2/(2*m*l^2) + m*g*l*(1-cos(θ))
end

g = 9.81
m = 2.
l = 1.
params = [g,m,l]
θ₀ = 1.
ℒ₀ = 1.

prob = HamiltonianProblem(H, θ₀, ℒ₀, (0., 100.), p=params)
sol2 = solve(prob, SofSpa10(), dt = .05);
#==========================================================#


#Plotting
#==========================================================#
plot(sol1, vars=1, tspan=(0,20), label="ODE Method")
plot!(sol2.t, sol2[1,:], xlim=(0,20), label="Hamiltonian Method")
#They produce the same solution!
