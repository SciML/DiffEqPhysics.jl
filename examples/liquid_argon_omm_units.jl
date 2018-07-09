using DiffEqPhysics

T = 120.0 # °K
kb = 8.3144598e-3 # kJ/(K*mol)
ϵ = T * kb
σ = 0.34 # nm
ρ = 1374/1.6747# Da/nm^3
m = 39.95# Da
N = 1000
L = (m*N/ρ)^(1/3)#10.229σ
R = 0.5*L   
v_dev = sqrt(kb * T / m)
bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

τ = 0.5e-3 # ps or 1e-12 s
t1 = 0.0
t2 = 10000τ

parameters = LennardJonesParameters(ϵ, σ, R)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
thermostat = BerendsenThermostat(90, 10τ)
thermostat = NoseHooverThermostat(90, 1)
pbc = CubicPeriodicBoundaryConditions(L)
simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)


using Plots
t = t1:τ:result.solution.t[end-1]
temper = temperature.(result, t)
plot(t, temper, ylim=[0,200], xlabel="t, ps", ylabel = "T, °K", label="Temperature, °K" )