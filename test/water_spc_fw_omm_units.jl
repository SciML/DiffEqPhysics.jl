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

const Na = 6.022e23
const T = 298.16 # °K
const kb = 8.3144598e-3 # kJ/(K*mol)
const ϵOO = 0.1554253*4.184 # kJ 
const σOO = 0.3165492 # nm
const ρ = 997/1.6747# Da/nm^3
const mO = 15.999 # Da
const mH = 1.00794 # Da
const mH2O = mO+2*mH
const N = 216#floor(Int, ρ * L^3 / m)
const L = (mH2O*N/ρ)^(1/3)#10.229σ
const R = 0.9 # ~3*σOO  
const Rel = 0.49*L
const v_dev = sqrt(3*kb * T /mH2O)
const τ = 0.5e-3 # ps
const t1 = 0τ
const t2 = 3000τ # ps
const k_bond = 1059.162*4.184*1e2 # kJ/(mol*nm^2)
const k_angle = 75.90*4.184 # kJ/(mol*rad^2)
const rOH = 0.1012 # nm
const ∠HOH = 113.24*pi/180 # rad
const qH = 0.41
const qO = -0.82
const k = 138.935458 #
bodies = generate_bodies_in_cell_nodes(N, mH2O, v_dev, L)
jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
e_parameters = ElectrostaticParameters(k, Rel)
spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
pbc = CubicPeriodicBoundaryConditions(L)
water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters);
simulation = NBodySimulation(water, (t1, t2), pbc, kb);
#result = run_simulation(simulation, Tsit5())
result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)


time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

(rs, grf) = @time rdf(result)
(ts, dr2) = @time msd(result)

using JLD
#save("D:/water $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2)
#save("D:/water $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2, "e_tot", e_tot, "e_kin", e_kin, "e_pot", e_pot)

using Plots
import GR
plot(rs, grf, xlim=[0, 0.4999L], label=["Radial distribution function"],ylabel="g(r)", xlabel="r, nm")

#@time animate(result, "D:/$Nactual H20 particles with $timesteps timesteps $time_now.gif")
