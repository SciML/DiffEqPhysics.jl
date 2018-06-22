include("../src/nbody_simulation.jl")

function generate_bodies_in_cell_nodes(n::Int, m::Real, v_dev::Real, L::Real)
   
    rng = MersenneTwister(n);
    velocities = v_dev * randn(rng, Float64, (3, n))
    bodies = MassBody[]

    count = 1
    dL = L / (ceil(n^(1 / 3)))
    for x = dL/2:dL:L, y = dL/2:dL:L, z = dL/2:dL:L  
    #for x = 1.5*dL:2dL:5.6*dL, y = 1.5*dL:2dL:5.6*dL, z = 1.5*dL:2dL:5.6*dL        
        if count > n
            break
        end
        r = SVector(x, y, z)
        v = SVector{3}(velocities[:,count])
        body = MassBody(r, v, m)
        push!(bodies, body)
        count += 1           
    end
    return bodies
end

function generate_bodies_in_line(n::Int, m::Real, v_dev::Real, L::Real)
    dL = L / (ceil(n^(1 / 3)))
    n_line = floor(Int, L / dL)
    rng = MersenneTwister(n);
    velocities = v_dev * randn(rng, Float64, (3, n_line))
    bodies = MassBody[]
    x = y = L / 2
    for i = [3] #1:n_line-1
        r = SVector(x, y, i * dL)
        v = SVector{3}(velocities[:,i])
        body = MassBody(r, v, m)
        push!(bodies, body)  
    end
    return bodies
end

const qe = 1.6e-19
const Na = 6.022e23
const T = 298.16 # °K
const kb = 1.38e-23 # J/K
const ϵOO = 0.1554253*4184/Na
const σOO = 3.165492e-10 # m
const ρ = 1000 # kg/m^3
const mO = 15.999 * 1.6747 * 1e-27 # kg
const mH = 1.00794 * 1.6747 * 1e-27#
const mH2O = mO+2*mH
const N = 216#floor(Int, ρ * L^3 / m)
const L = (mH2O*N/ρ)^(1/3)#10.229σ
const R = 9e-10 # 3*σOO  
const Rel = 0.45*L
const v_dev = 0#sqrt(kb * T /mH2O)
const τ = 1e-15 # σ/v
const t1 = 0τ
const t2 = 15τ 
const k_bond = 1059.162*4184*1e20/Na # J/m^2
const k_angle = 75.90*4184/Na # J/rad^2
const rOH = 1.012e-10 # m
const ∠HOH = 113.24*pi/180 # rad
const qH = 0.41*qe
const qO = -0.84*qe
const k = 9e9 #
#bodies = generate_bodies_randomly(N, m, v_dev, L)
bodies = generate_bodies_in_cell_nodes(N, mH2O, v_dev, L)
#bodies = generate_bodies_in_line(N, mH2O, v_dev, L)
jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
e_parameters = ElectrostaticParameters(k, Rel)
spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
pbc = CubicPeriodicBoundaryConditions(L)
water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters);
simulation = NBodySimulation(water, (t1, t2), pbc);
#result = run_simulation(simulation, Tsit5())
result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)

using Plots
import GR
time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))
#@time animate(result, "D:/$Nactual H20 particles with $timesteps timesteps $time_now.gif")

#plot(simulation)