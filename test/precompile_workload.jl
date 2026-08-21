using DiffEqPhysics, Test

@testset "Precompile workload" begin
    H(p, q, _, _) = p^2 / 2 + q^2 / 2
    @test HamiltonianProblem(H, 1.0, 0.0, (0.0, 1.0)).tspan == (0.0, 1.0)

    dp(p, q, _, _) = -q
    dq(p, q, _, _) = p
    @test HamiltonianProblem((dp, dq), 1.0, 0.0, (0.0, 1.0)).tspan == (0.0, 1.0)
end
