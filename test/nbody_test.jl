using RecursiveArrayTools

@testset "NBodyProblem Test" begin
    G = 2.95912208286e-4
    M = [
        1.00000597682,
        0.000954786104043,
        0.000285583733151,
        0.0000437273164546,
        0.0000517759138449,
        1 / 1.3e8,
    ]
    invM = inv.(M)
    planets = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

    pos_x = [0.0, -3.5023653, 9.0755314, 8.3101420, 11.4707666, -15.5387357]
    pos_y = [0.0, -3.8169847, -3.0458353, -16.2901086, -25.7294829, -25.2225594]
    pos_z = [0.0, -1.5507963, -1.6483708, -7.2521278, -10.8169456, -3.1902382]
    pos = ArrayPartition(pos_x, pos_y, pos_z)

    vel_x = [0.0, 0.00565429, 0.00168318, 0.00354178, 0.00288930, 0.00276725]
    vel_y = [0.0, -0.00412490, 0.00483525, 0.00137102, 0.00114527, -0.00170702]
    vel_z = [0.0, -0.00190589, 0.00192462, 0.00055029, 0.00039677, -0.00136504]
    vel = ArrayPartition(vel_x, vel_y, vel_z)

    tspan = (0.0, 200_000)

    ∑ = sum
    N = 6
    function potential(p, t, x, y, z, M)
        -G * ∑(i -> ∑(j -> (M[i] * M[j]) /
                      sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2), 1:(i - 1)),
          2:N)
    end
    nprob = NBodyProblem(potential, M, vel, pos, tspan)
    sol = solve(nprob, Yoshida6(), dt = 100)

    # Make sure the distances from the sun stay small enough
    f = (x, y, z, j) -> sqrt((x[1] - x[j])^2 + (y[1] - y[j])^2 + (z[1] - z[j])^2)
    for i in 1:length(sol)
        x, y, z = sol[i].x[1].x
        @test f(x, y, z, 2) < 6
        @test f(x, y, z, 3) < 10.5
        @test f(x, y, z, 4) < 21
        @test f(x, y, z, 5) < 32
        @test f(x, y, z, 6) < 50
    end
end
