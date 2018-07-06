println("====  Andersen thermostat testing    ====")

let 
    T = 120.0 # °K
    kb = 1.38e-23 # J/K
    ϵ = T * kb
    σ = 3.4e-10 # m
    ρ = 1374 # kg/m^3
    m = 39.95 * 1.6747 * 1e-27 # kg
    L = 5σ # 10.229σ
    N = 3 # floor(Int, ρ * L^3 / m)
    R = 2.25σ   
    v_dev = sqrt(3*kb * T / m)
    r1 = SVector(L/3, L/3, 2*L/3)
    r2 = SVector(L/3, 2*L/3, L/3)
    r3 = SVector(2*L/3, L/3, L/3)
    v1 = SVector(0, 0, -v_dev)
    v2 = SVector(0, -v_dev, 0)
    v3 = SVector(-v_dev, 0, 0)
    p1 = MassBody(r1, v1, m)
    p2 = MassBody(r2, v2, m)
    p3 = MassBody(r3, v3, m)

    τ = 1e-14
    t1 = 0.0
    t2 = 200τ

    parameters = LennardJonesParameters(ϵ, σ, R)
    lj_system = PotentialNBodySystem([p1, p2, p3], Dict(:lennard_jones => parameters));
    thermostat = AndersenThermostat(0.2, T, kb)
    simulation = NBodySimulation(lj_system, (t1, t2), PeriodicBoundaryConditions(L), thermostat);
    result = run_simulation(simulation, VelocityVerlet(), dt=τ)
    

    T1 = temperature(result, t1) 
    T2 = temperature(result, t2)
    ε = 2
    @test abs(T2-T1)/T1 ≈ 1.0 atol = ε
end