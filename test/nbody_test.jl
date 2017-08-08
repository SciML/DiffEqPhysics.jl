using DiffEqPhysics, ForwardDiff, Base.Test

G = 2.95912208286e-4
M = Diagonal(repeat([1.00000597682, 0.000954786104043, 0.000285583733151,
        0.0000437273164546, 0.0000517759138449, 1/1.3e8], inner=3))       # repeat each mass three time
                                                                          # becasue there are 3 degrees
                                                                          # of freedom for each planet/sun

invM = inv(M).diag    # the diagonal of the inverse of the mass matrix
                      # will explain it later

planets = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

pos = [
    0.,0.,0.,                              -3.5023653,-3.8169847,-1.5507963,
    9.0755314,-3.0458353,-1.6483708,       8.3101420,-16.2901086,-7.2521278,
    11.4707666,-25.7294829,-10.8169456,    -15.5387357,-25.2225594,-3.1902382
]
vel = [
    0,0,0,                               0.00565429,-0.00412490,-0.00190589,
    0.00168318,0.00483525,0.00192462,    0.00354178,0.00137102,0.00055029,
    0.00288930,0.00114527,0.00039677,    0.00276725, -0.00170702, -0.00136504
]

tspan = (0.,20)

const ∑ = sum
const N = 5

ø = i -> (3i+1):(3(i+1))    # Generate indices like 1:3, 4:6 ...
V(q) = -G*∑(i->∑(j->(M[3i+1,3i+1]*M[3j+1,3j+1])/norm(q[ø(i)]-q[ø(j)]), 0:i-1), 1:N)

function acceleration(t, x, v, dv)
    ForwardDiff.gradient!(dv, V, x)
    for i in eachindex(dv)
        dv[i] *= -invM[i]
    end
end

prob = SecondOrderODEProblem(acceleration, pos, vel, tspan)

mass = [1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 1/1.3e8]
potential(t, x, y, z, M) = -G*∑(i->∑(j->(M[i+1]*M[j+1])/sqrt((x[i+1]-x[j+1])^2 + (y[i+1]-y[j+1])^2 + (z[i+1]-z[j+1])^2), 0:i-1), 1:N)
nprob = NBodyProblem(potential, mass, pos, vel, tspan)

@test_broken test_solve(prob, nprob)
