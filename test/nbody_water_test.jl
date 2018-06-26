let
    const qe = 1.6e-19
    const Na = 6.022e23
    const T = 298.16 # °K
    const kb = 1.38e-23 # J/K
    const ϵOO = 0.1554253*4184/Na
    const σOO = 3.165492e-10 # m
    const ρ = 900 # kg/m^3
    const mO = 15.999 * 1.6747 * 1e-27 # kg
    const mH = 1.00794 * 1.6747 * 1e-27#
    const mH2O = mO+2*mH
    const N = 216#floor(Int, ρ * L^3 / m)
    const L = (mH2O*N/ρ)^(1/3)#10.229σ
    const R = 9e-10 # 3*σOO  
    const Rel = 0.45*L
    const v_dev = sqrt(kb * T /mH2O)
    const τ = 1e-15 # σ/v
    const t1 = 0τ
    const t2 = 10τ 
    const k_bond = 1059.162*4184*1e20/Na # J/m^2
    const k_angle = 75.90*4184/Na # J/rad^2
    const rOH = 1.012e-10 # m
    const ∠HOH = 113.24*pi/180 # rad
    const qH = 0.41*qe
    const qO = -0.84*qe
    const k = 9e9 #

    r1 = SVector(L/3, L/3, 2*L/3)
    r2 = SVector(L/3, 2*L/3, L/3)
    r3 = SVector(2*L/3, L/3, L/3)
    v1 = SVector(0, 0, -v_dev)
    v2 = SVector(0, -v_dev, 0)
    v3 = SVector(-v_dev, 0, 0)
    p1 = MassBody(r1, v1, mH2O)
    p2 = MassBody(r2, v2, mH2O)
    p3 = MassBody(r3, v3, mH2O)

    bodies = [p1, p2, p3]
    jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
    e_parameters = ElectrostaticParameters(k, Rel)
    spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
    pbc = CubicPeriodicBoundaryConditions(L)
    water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters);
    simulation = NBodySimulation(water, (t1, t2), pbc);

    result = run_simulation(simulation, VelocityVerlet(), dt=τ)

    cc = get_position(result, t2)
    ε = 0.01
    for i=1:3
        indO = 3*(i-1)+1
        @test (rOH - norm(cc[:,indO]-cc[:,indO+1]))/rOH ≈ 0.0 atol = ε
        @test (rOH - norm(cc[:,indO]-cc[:,indO+2]))/rOH ≈ 0.0 atol = ε
        ang = acos(dot(cc[:,indO]-cc[:,indO+1],cc[:,indO]-cc[:,indO+2])/(norm(cc[:,indO]-cc[:,indO+1])*norm(cc[:,indO]-cc[:,indO+2])))
        @test (ang - ∠HOH)/∠HOH ≈ 0.0 atol = ε
    end

    e_tot_1 = total_energy.(result, t1)
    e_tot_2 = total_energy.(result, t2)
    #@test (e_tot_1 - e_tot_2)/e_tot_1 ≈ 0.0 atol = ε
end