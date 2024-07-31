using Test
using DiffEqPhysics, ForwardDiff, OrdinaryDiffEq
using StaticArrays, LinearAlgebra, Random

test_solve(prob...) = mapreduce(p->solve(p, Tsit5(), dt=1//2).u, ==, prob)

p0, q0   = rand(2)
H(dθ, θ, p) = dθ / 2 - 9.8 * cos(θ)
dp(dθ, θ, p, t) = - 9.8 * sin(θ)
dq(dθ, θ, p, t) = 0.5
acc      = (v, x, p, t) -> ForwardDiff.derivative(x -> -H(v[1], x, p), x[1])
vel      = (v, x, p, t) -> ForwardDiff.derivative(v -> H(v, x[1], p), v[1])
prob_1   = DynamicalODEProblem(acc, vel, p0, q0, (0., 10.))

@testset "prob1($h)" for h in (H, (dp, dq))
    prob1    = HamiltonianProblem(h, p0, q0, (0., 10.))
    @test test_solve(prob1, prob_1)
end

p0 = @SVector rand(2)
q0 = @SVector rand(2)
H4(dθ, θ, p, t) = dθ[1] / 2 + dθ[2] / 2 - 9.8 * cos(θ[1]) - 9.8 * cos(θ[2]) + t
dq(dθ, θ, p, t) = @SVector [0.5, 0.5]
dp(dθ, θ, p, t) = @SVector [-9.8 * sin(θ[1]), -9.8 * sin(θ[2])]
acc      = (v, x, p, t) -> ForwardDiff.gradient(x -> -H4(v, x, p, t), x)
vel      = (v, x, p, t) -> ForwardDiff.gradient(v -> H4(v, x, p, t), v)
prob_2   = DynamicalODEProblem(acc, vel, p0, q0, (0., 10.))

@testset "prob2($h)" for h in (H4, (dp, dq))
    prob2    = HamiltonianProblem(h, p0, q0, (0., 10.))
    @test test_solve(prob2, prob_2)
end

H(dθ, θ, p, t) = (dθ / 2 - 9.8 * cos.(θ))[1]
p0, q0 = [rand(5) for i in 1:2]
dq(dx, dΘ, θ, p, t) = (dx .= [0.5, 0, 0, 0, 0])
dp(dv, dθ, θ, p, t) = (dv .= [(-9.8 * sin.(θ))[1], 0, 0, 0, 0])
acc = (dv, v, x, p, t) -> ForwardDiff.gradient!(dv, x -> -H(v, x, p), x)
vel = (dx, v, x, p, t) -> ForwardDiff.gradient!(dx, v -> H(v, x, p), v)
prob_3 = DynamicalODEProblem(acc, vel, p0, q0, (0., 10.))

@testset "prob3($h)" for h in (H, (dp, dq))
    prob3  = HamiltonianProblem(h, p0, q0, (0., 10.))
    @test test_solve(prob3, prob_3)
end

H4(dθ, θ, p, t) = (dθ / 2 - 9.8 * cos.(θ))[1] + t
p0, q0 = [rand(5) for i in 1:2]
dq(dx, dΘ, θ, p, t) = (dx .= [0.5, 0, 0, 0, 0])
dp(dv, dθ, θ, p, t) = (dv .= [(-9.8 * sin.(θ))[1], 0, 0, 0, 0])
acc = (dv, v, x, p, t) -> ForwardDiff.gradient!(dv, x -> -H4(v, x, p, t), x)
vel = (dx, v, x, p, t) -> ForwardDiff.gradient!(dx, v -> H4(v, x, p, t), v)
prob_4 = DynamicalODEProblem(acc, vel, p0, q0, (0., 10.))

@testset "prob4($h)" for h in (H4, (dp, dq))
    prob4  = HamiltonianProblem(h, p0, q0, (0., 10.))
    @test test_solve(prob4, prob_4)
end
