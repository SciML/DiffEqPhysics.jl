using DiffEqPhysics
using Base.Test

println("====    HamiltonianProblem Test    ====")
println("====         Scalar case test      ====")
H(θ, dθ) = dθ/2 - 9.8*cos(θ)
q0, p0   = rand(2)
acc      = (t, x, v, dv) -> dv[1] = ForwardDiff.derivative(x->-H(x,v[1]), x[1])
vel      = (t, x, v, dx) -> dx[1] = ForwardDiff.derivative(v-> H(x[1],v), v[1])

prob1    = HamiltonianProblem(H, q0, p0, (0.,10.))
prob_1   = ODEProblem((vel,acc), ([q0],[p0]), (0.,10.))
test_solve(prob...) = mapreduce(p->solve(p, VelocityVerlet(), dt=1//2).u, ==, prob)
@test test_solve(prob1, prob_1)

println("====         Vector case test      ====")
H(θ, dθ) = (dθ/2 - 9.8*cos.(θ))[1]
q0, p0 = [rand(5) for i in 1:2]
acc = (t, x, v, dv) -> ForwardDiff.gradient!(dv, x->-H(x,v), x)
vel = (t, x, v, dx) -> ForwardDiff.gradient!(dx, v-> H(x,v), v)
prob2  = HamiltonianProblem(H, q0, p0, (0.,10.))
prob_2 = ODEProblem((vel,acc), (q0,p0), (0.,10.))
@test test_solve(prob2, prob_2)
