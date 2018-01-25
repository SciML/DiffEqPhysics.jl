using DiffEqPhysics, ForwardDiff, StaticArrays
using Base.Test

println("====    HamiltonianProblem Test    ====")
println("====         Scalar case test      ====")
H(dθ ,θ, p) = dθ/2 - 9.8*cos(θ)
p0, q0   = rand(2)
acc      = (v, x, p, t) -> ForwardDiff.derivative(x->-H(v[1],x,p), x[1])
vel      = (v, x, p, t) -> ForwardDiff.derivative(v-> H(v,x[1],p), v[1])

prob1    = HamiltonianProblem(H, p0, q0, (0.,10.))
prob_1   = DynamicalODEProblem(acc,vel,p0,q0,(0.,10.))

@test test_solve(prob1, prob_1)

println("====   Static Vector case test      ====")
p0 = @SVector rand(2)
q0 = @SVector rand(2)
H(dθ, θ, p) = dθ[1]/2 + dθ[2]/2 - 9.8*cos(θ[1]) - 9.8*cos(θ[2])
prob1    = HamiltonianProblem(H, p0, q0, (0.,10.))
acc      = (v, x, p, t) -> ForwardDiff.gradient(x->-H(v,x,p), x)
vel      = (v, x, p, t) -> ForwardDiff.gradient(v-> H(v,x,p), v)
prob_1   = DynamicalODEProblem(acc,vel,p0,q0,(0.,10.))
@test test_solve(prob1, prob_1)

println("====         Vector case test      ====")
H(dθ, θ, p) = (dθ/2 - 9.8*cos.(θ))[1]
p0, q0 = [rand(5) for i in 1:2]
acc = (dv, v, x, p, t) -> ForwardDiff.gradient!(dv, x->-H(v,x,p), x)
vel = (dx, v, x, p, t) -> ForwardDiff.gradient!(dx, v-> H(v,x,p), v)
prob2  = HamiltonianProblem(H, p0, q0, (0.,10.))
prob_2 = DynamicalODEProblem(acc,vel,p0,q0, (0.,10.))
@test test_solve(prob2, prob_2)