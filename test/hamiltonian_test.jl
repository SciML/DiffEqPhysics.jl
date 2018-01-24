using DiffEqPhysics, ForwardDiff, StaticArrays
using Base.Test

println("====    HamiltonianProblem Test    ====")
println("====         Scalar case test      ====")
H(θ, dθ, p) = dθ/2 - 9.8*cos(θ)
q0, p0   = rand(2)
acc      = (x, v, p, t) -> ForwardDiff.derivative(x->-H(x,v[1],p), x[1])
vel      = (x, v, p, t) -> ForwardDiff.derivative(v-> H(x[1],v,p), v[1])

prob1    = HamiltonianProblem(H, q0, p0, (0.,10.))
prob_1   = DynamicalODEProblem(vel,acc,q0,p0, (0.,10.))
@test test_solve(prob1, prob_1)

println("====   Static Vector case test      ====")
q0 = @SVector rand(2)
p0 = @SVector rand(2)
H(θ, dθ, p) = dθ[1]/2 + dθ[2]/2 - 9.8*cos(θ[1]) - 9.8*cos(θ[2])
prob1    = HamiltonianProblem(H, q0, p0, (0.,10.))
acc      = (x, v, p, t) -> ForwardDiff.gradient(x->-H(x,v,p), x)
vel      = (x, v, p, t) -> ForwardDiff.gradient(v-> H(x,v,p), v)
prob_1   = DynamicalODEProblem(vel,acc,q0,p0, (0.,10.))
@test test_solve(prob1, prob_1)

println("====         Vector case test      ====")
H(θ, dθ, p) = (dθ/2 - 9.8*cos.(θ))[1]
q0, p0 = [rand(5) for i in 1:2]
acc = (dv, x, v, p, t) -> ForwardDiff.gradient!(dv, x->-H(x,v,p), x)
vel = (dx, x, v, p, t) -> ForwardDiff.gradient!(dx, v-> H(x,v,p), v)
prob2  = HamiltonianProblem(H, q0, p0, (0.,10.))
prob_2 = DynamicalODEProblem(vel,acc,q0,p0, (0.,10.))
@test test_solve(prob2, prob_2)
