# SimulationResult sould provide an interface for working with properties of a separate particle
# and with physical properties of the whole system.
struct SimulationResult{sType <: NBodySystem}
    solution::AbstractTimeseriesSolution
    simulation::NBodySimulation{sType}
end

function Base.show(stream::IO, sr::SimulationResult)
    print(stream, "N: ") 
    show(stream, length(sr.simulation.system.bodies))
    println(stream)
    show(stream, sr.simulation)
    print(stream, "Time steps: ") 
    show(stream, length(sr.solution.t))
    println(stream)
    print(stream, "t: ", minimum(sr.solution.t), ", ", maximum(sr.solution.t)) 
    println(stream)
end

(sr::SimulationResult)(args...; kwargs...) = return sr.solution(args...; kwargs...)


# This iterator interface is implemented specifically for making animation.
# Probably, there will be a wrapper for this in the future.
Base.start(::SimulationResult) = 1

Base.done(sr::SimulationResult, state) = state > length(sr.solution.t)

function Base.next(sr::SimulationResult, state) 
    (sr, sr.solution.t[state]), state + 1
end

function get_velocity(sr::SimulationResult, time::Real, i::Integer=0)
    if typeof(sr.solution[1]) <: RecursiveArrayTools.ArrayPartition
        velocities = sr(time).x[1]
        n = size(velocities, 2)
        if i <= 0
            return velocities[:, 1:end]
        else
            return velocities[:, i]
        end
    else
        velocities = sr(time)
        n = div(size(velocities, 2), 2)
        if i <= 0
            return @view velocities[:, n + 1:end]
        else
            return @view velocities[:, n + i]
        end
    end
end

function get_position(sr::SimulationResult, time::Real, i::Integer=0)
    if typeof(sr.solution[1]) <: RecursiveArrayTools.ArrayPartition
        positions = sr(time).x[2]
        n = size(positions, 2)
    else
        positions = sr(time)
        n = div(size(positions, 2), 2)
    end

    if i <= 0
        return @view positions[:, 1:n]
    else
        return @view positions[:, i]
    end
end

function get_masses(system::NBodySystem)
    n = length(system.bodies)
    masses = zeros(n)
    for i = 1:n
        masses[i] = system.bodies[i].m
    end
    return masses
end

function get_masses(system::WaterSPCFw)
    n = length(system.bodies)
    ms = zeros(Real, 3 * n)
    for i = 1:n
        ms[3 * (i - 1) + 1] = system.mO
        ms[3 * (i - 1) + 2] = system.mH
        ms[3 * (i - 1) + 3] = system.mH
    end 
    return ms
end

function temperature(result::SimulationResult, time::Real)
    kb = result.simulation.kb
    velocities = get_velocity(result, time)
    masses = get_masses(result.simulation.system)
    temperature = mean(sum(velocities.^2, 1) .* masses) / (3kb)
    return temperature
end

function temperature(result::SimulationResult{<:WaterSPCFw}, time::Real)
    kb = result.simulation.kb
    system = result.simulation.system
    n = length(system.bodies)
    vs = get_velocity(result, time)
    mH2O = system.mO + 2 * system.mH
    v2 = zeros(n)
    for i = 1:n
        indO, indH1, indH2 = 3 * (i - 1) + 1, 3 * (i - 1) + 2, 3 * (i - 1) + 3
        v_c = @. (vs[:,indO] * system.mO + vs[:,indH1] * system.mH + vs[:,indH2] * system.mH) / mH2O
        v2[i] = dot(v_c, v_c)
    end
    temperature = mean(v2) * mH2O / (3kb)
    return temperature
end

function kinetic_energy(velocities, masses)
    ke = sum(dot(vec(sum(velocities.^2, 1)), masses / 2))
    return ke
end

function kinetic_energy(sr::SimulationResult, time::Real)
    vs = get_velocity(sr, time)
    ms = get_masses(sr.simulation.system)
    return kinetic_energy(vs, ms)
end

function potential_energy(coordinates, simulation::NBodySimulation)
    e_potential = 0
    system = simulation.system
    n = length(system.bodies)
    if :lennard_jones ∈ keys(system.potentials)
        p = system.potentials[:lennard_jones]
        (ms, indxs) = obtain_data_for_lennard_jones_interaction(system)
        e_potential += lennard_jones_potential(p, indxs, coordinates, simulation.boundary_conditions)
    end

    if :electrostatic ∈ keys(system.potentials)
        p = system.potentials[:electrostatic]
        (qs, ms, indxs, exclude) = obtain_data_for_electrostatic_interaction(simulation.system)
        e_potential += electrostatic_potential(p, indxs, exclude, qs, coordinates, simulation.boundary_conditions)
    end
    e_potential
end

function potential_energy(coordinates, simulation::NBodySimulation{<:WaterSPCFw})
    e_potential = 0
    system = simulation.system
    p = system.lj_parameters
    (ms, indxs) = obtain_data_for_lennard_jones_interaction(system)
    e_potential += lennard_jones_potential(p, indxs, coordinates, simulation.boundary_conditions)

    p = system.e_parameters
    (qs, ms, indxs, exclude) = obtain_data_for_electrostatic_interaction(simulation.system)
    e_potential += electrostatic_potential(p, indxs, exclude, qs, coordinates, simulation.boundary_conditions)

    p = system.scpfw_parameters
    (ms, neighbouhoods) = obtain_data_for_harmonic_bond_interaction(simulation.system, p)
    e_potential += harmonic_bonds_potential(p, coordinates, ms, neighbouhoods)

    p = system.scpfw_parameters
    (ms, bonds) = obtain_data_for_valence_angle_harmonic_interaction(simulation.system)
    e_potential += valence_angle_harmonic_potential(coordinates, bonds)
    e_potential
end

function lennard_jones_potential(p::LennardJonesParameters, indxs::Vector{<:Integer}, coordinates, pbc::BoundaryConditions)
    e_lj = 0
    n = length(indxs)
    @inbounds for ind_i = 1:n
        i = indxs[ind_i]
        ri = @SVector [coordinates[1, i], coordinates[2, i], coordinates[3, i]]
        for ind_j = ind_i + 1:n    
            j = indxs[ind_j]          
            rj = @SVector [coordinates[1, j], coordinates[2, j], coordinates[3, j]]
            
            (rij, rij_2, success) = apply_boundary_conditions!(ri, rj, pbc, p.R2)
            
            if !success
                rij_2 = p.R2
            end
            
            σ_rij_6 = (p.σ2 / rij_2)^3
            σ_rij_12 = σ_rij_6^2
            e_lj += (σ_rij_12 - σ_rij_6 )
        end
    end 

    return 4 * p.ϵ * e_lj
end

function electrostatic_potential(p::ElectrostaticParameters, indxs::Vector{<:Integer}, exclude::Dict{Int,Vector{Int}}, qs, rs, pbc::BoundaryConditions)
    e_el = 0

    n = length(indxs)
    for ind_i = 1:n
        i = indxs[ind_i]
        ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
        e_el_i = 0
        for ind_j = ind_i + 1:n    
            j = indxs[ind_j] 
            if !in(j, exclude[i])
                rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]

                (rij, rij_2, success) = apply_boundary_conditions!(ri, rj, pbc, p.R2)
                if !success
                    rij = p.R
                end
                
                e_el_i += qs[j] / norm(rij)
            end
        end    
        e_el += e_el_i * qs[i]
    end

    return e_el * p.k
end

function harmonic_bonds_potential(p::SPCFwParameters,
    rs,
    ms::Vector{<:Real},
    neighborhoods::Dict{Int,Vector{Tuple{Int,Float64}}})

    e_harmonic = 0
    
    @inbounds for (i, neighborhood) ∈ neighborhoods
        ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]
        for (j, k) in neighborhood
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            rij = ri - rj
            r = norm(rij)
            d = r - p.rOH
            e_harmonic += d^2 * k
        end
    end
    return e_harmonic / 4
end

function valence_angle_harmonic_potential(rs,
    bonds::Vector{Tuple{Int,Int,Int,Float64,Float64}})

    e_valence = 0

    for (a, b, c, valence_angle, k) ∈ bonds
        ra = @SVector [rs[1, a], rs[2, a], rs[3, a]]
        rb = @SVector [rs[1, b], rs[2, b], rs[3, b]]
        rc = @SVector [rs[1, c], rs[2, c], rs[3, c]]

        rba = ra - rb
        rbc = rc - rb

        currenct_angle = acos(dot(rba, rbc) / (norm(rba) * norm(rbc)))        
        e_valence += k * (currenct_angle - valence_angle)^2
    end
    return e_valence / 2
end

function potential_energy(sr::SimulationResult, time::Real)
    e_potential = 0
    coordinates = get_position(sr, time)
    return potential_energy(coordinates, sr.simulation)
end

function total_energy(sr::SimulationResult, time::Real)
    e_kin = kinetic_energy(sr, time)
    e_pot = potential_energy(sr, time)
    e_kin + e_pot
end

function initial_energy(simulation::NBodySimulation)
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation.system)
    ms = get_masses(simulation.system)
    return potential_energy(u0, simulation) + kinetic_energy(v0, ms) 
end

# Instead of treating NBodySimulation as a DiffEq problem and passing it into a solve method
# it is better to use a specific function for n-body simulations.
function run_simulation(s::NBodySimulation, alg_type=Tsit5(), args...; kwargs...)
    solution = solve(ODEProblem(s), alg_type, args...; kwargs...)
    return SimulationResult(solution, s)
end

# this should be a method for integrators designed for the SecondOrderODEProblem (It is worth somehow to sort them from other algorithms)
function run_simulation(s::NBodySimulation, alg_type::Union{VelocityVerlet,DPRKN6,Yoshida6}, args...; kwargs...)
    cb = obtain_callbacks_for_so_ode_problem(s)
    solution = solve(SecondOrderODEProblem(s), alg_type, args...; callback=cb, kwargs...)
    return SimulationResult(solution, s)
end

function obtain_callbacks_for_so_ode_problem(s::NBodySimulation)
    callback_array = Vector{DECallback}()

    if s.thermostat isa AndersenThermostat
        push!(callback_array, get_andersen_thermostating_callback(s))
    end

    return CallbackSet(tuple(callback_array...)...)
end

function get_andersen_thermostating_callback(s::NBodySimulation)
    p = s.thermostat::AndersenThermostat
    n = length(s.system.bodies)
    v_dev = sqrt(p.kb * p.T / s.system.bodies[1].m)

    condition = function (u, t, integrator)
        true
    end
    affect! = function (integrator)
        for i = 1:n
            if randn() < p.ν * (integrator.t - integrator.tprev)
                @. integrator.u.x[1][:,i] = v_dev * randn()
            end
        end
    end
    cb = DiscreteCallback(condition, affect!)
end

@recipe function generate_data_for_scatter(sr::SimulationResult{<:PotentialNBodySystem}, time::Real=0.0)
    solution = sr.solution
    n = length(sr.simulation.system.bodies)

    if :gravitational ∈ keys(sr.simulation.system.potentials)
    
        xlim --> 1.1 * [minimum(solution[1,1:n,:]), maximum(solution[1,1:n,:])]
        ylim --> 1.1 * [minimum(solution[2,1:n,:]), maximum(solution[2,1:n,:])]        
    
        for i in 1:n
            @series begin
                label --> "Orbit $i"
                vars --> (3 * (i - 1) + 1, 3 * (i - 1) + 2)
                solution
            end
        end
    else
        positions = get_position(sr, time)
        borders = sr.simulation.boundary_conditions
        if borders isa PeriodicBoundaryConditions
            xlim --> 1.1 * [borders[1], borders[2]]
            ylim --> 1.1 * [borders[3], borders[4]]
            zlim --> 1.1 * [borders[5], borders[6]]
        elseif borders isa CubicPeriodicBoundaryConditions
            xlim --> 1.1 * [0, borders.L]
            ylim --> 1.1 * [0, borders.L]
            zlim --> 1.1 * [0, borders.L]
            map!(x ->  x -= borders.L * floor(x / borders.L), positions, positions)    
        end

        
        seriestype --> :scatter
        markersize --> 5

        (positions[1,:], positions[2,:], positions[3,:])
        #(positions[1,:], positions[2,:])
    end
end

function distancies(result::SimulationResult, time::Real)
    n = length(result.simulation.system.bodies)
    cc = get_position(result, time)

    d = Float64[]
    for i = 1:n
        for j = 1:n
            if i != j
                push!(d, norm(cc[:,i] - cc[:,j]))
            end
        end
    end
    return d
end

@recipe function initial_distribution(sr::SimulationResult{<:WaterSPCFw}, time::Real=0.0)

    n = length(sr.simulation.system.bodies)

    borders = sr.simulation.boundary_conditions
    
    cc = get_position(sr, time)
 
    if borders isa PeriodicBoundaryConditions
        xlim --> 1.1 * [borders[1], borders[2]]
        ylim --> 1.1 * [borders[3], borders[4]]
        zlim --> 1.1 * [borders[5], borders[6]]
    elseif borders isa CubicPeriodicBoundaryConditions
        xlim --> 1.1 * [0, borders.L]
        ylim --> 1.1 * [0, borders.L]
        zlim --> 1.1 * [0, borders.L]

        
        map!(x ->  x -= borders.L * floor(x / borders.L), cc, cc)
    end
    seriestype --> :scatter

    @series begin
        label --> "O"
        markersize --> 8
        markercolor --> :red
        (cc[1,1:3:3 * n - 2], cc[2,1:3:3 * n - 2], cc[3,1:3:3 * n - 2])
    end

    @series begin
        label --> "H"
        markersize --> 4
        markercolor --> :green
        x = vcat(cc[1,2:3:3 * n - 1], cc[1,3:3:3 * n])
        y = vcat(cc[2,2:3:3 * n - 1], cc[2,3:3:3 * n])
        z = vcat(cc[3,2:3:3 * n - 1], cc[3,3:3:3 * n])
        (x, y, z)
    end
end

function rdf(sr::SimulationResult)
    n = length(sr.simulation.system.bodies)
    pbc = sr.simulation.boundary_conditions
    
    (ms, indxs) = obtain_data_for_lennard_jones_interaction(sr.simulation.system)
    indlen = length(indxs)

    maxbin = 1000
    dr = pbc.L / maxbin
    hist = zeros(maxbin)
    for t ∈ sr.solution.t
        cc = get_position(sr, t)
        for ind_i = 1:indlen
            i = indxs[ind_i]
            ri = @SVector [cc[1, i], cc[2, i], cc[3, i]]
            for ind_j = ind_i + 1:indlen    
                j = indxs[ind_j]    
                rj = @SVector [cc[1, j], cc[2, j], cc[3, j]]

                (rij, rij_2, success) = apply_boundary_conditions!(ri, rj, pbc, (0.5 * pbc.L)^2)

                if success
                    rij_1 = sqrt(rij_2)
                    bin = ceil(Int, rij_1 / dr)
                    if bin > 1 && bin <= maxbin
                        hist[bin] += 2 
                    end
                end
            end
        end
    end

    c = 4 / 3 * π * indlen / pbc.L^3

    gr = zeros(maxbin)
    rs = zeros(maxbin)
    tlen = length(sr.solution.t)
    for bin = 1:maxbin
        rlower = (bin - 1) * dr
        rupper = rlower + dr
        nideal = c * (rupper^3 - rlower^3)
        gr[bin] = (hist[bin] / (tlen * indlen)) / nideal
        rs[bin] = rlower + dr / 2
    end

    return (rs, gr)
end

function msd(sr::SimulationResult{<:PotentialNBodySystem})
    n = length(sr.simulation.system.bodies)
    
    (ms, indxs) = obtain_data_for_lennard_jones_interaction(sr.simulation.system)
    indlen = length(indxs)

    ts = sr.solution.t
    tlen = length(ts)
    dr2 = zeros(length(ts))

    cc0 = get_position(sr, ts[1])

    for t = 1:tlen
        cc = get_position(sr, ts[t])
        for ind_i = 1:indlen
            i = indxs[ind_i]
            dr = @SVector [cc[1, i] - cc0[1,i], cc[2, i] - cc0[2,i], cc[3, i] - cc0[3,i]]
            
            dr2[t] += dot(dr, dr)
        end
        dr2[t] /= indlen
    end

    (ts, dr2)
end

function msd(sr::SimulationResult{<:WaterSPCFw})
    n = length(sr.simulation.system.bodies)
    
    mO = sr.simulation.system.mO
    mH = sr.simulation.system.mH

    ts = sr.solution.t
    tlen = length(ts)
    dr2 = zeros(tlen)

    cc0 = get_position(sr, ts[1])
            
    for t = 1:tlen
        cc = get_position(sr, ts[t])
        for i = 1:n
            indO, indH1, indH2 = 3 * (i - 1) + 1, 3 * (i - 1) + 2, 3 * (i - 1) + 3
            dr = ((cc[:,indO] - cc0[:,indO]) * mO + (cc[:,indH1] - cc0[:,indH1]) * mH + (cc[:,indH2] - cc0[:,indH2]) * mH) / (2 * mH + mO)
            
            dr2[t] += dot(dr, dr)
        end
        dr2[t] /= n
    end

    (ts, dr2)
end

# coordinates should be in angstroms
function save_to_pdb(sr::SimulationResult{<:WaterSPCFw}, path )
    f = open(path,"w+")
    L = 10*sr.simulation.boundary_conditions.L
    strL = @sprintf("%9.3f",L)
    strA = @sprintf("%7.2f",90)
    #println(f,"CRYST1",lpad(strL,9),lpad(strL,9),lpad(strL,9),lpad(strA,7),lpad(strA,7),lpad(strA,7))
    n = length(sr.simulation.system.bodies)
    count = 0
    for t in sr.solution.t
        cc = 10*get_position(sr, t)
        map!(x ->  x -= L * floor(x / L), cc, cc)
        count+=1      
        println(f,rpad("MODEL",10), count)
        println(f,"REMARK 250 time=$t picoseconds")
        for i ∈ 1:n
            indO, indH1, indH2 = 3 * (i - 1) + 1, 3 * (i - 1) + 2, 3 * (i - 1) + 3
            
            println(f,"HETATM",lpad(indO,5),"  ",rpad("O",4),"HOH",lpad(i,6),"    ", lpad(@sprintf("%8.3f",cc[1,indO]),8), lpad(@sprintf("%8.3f",cc[2,indO]),8), lpad(@sprintf("%8.3f",cc[3,indO]),8),lpad("1.00",6),lpad("0.00",6), lpad("",10), lpad("O",2))
            println(f,"HETATM",lpad(indH1,5),"  ",rpad("H1",4),"HOH",lpad(i,6),"    ", lpad(@sprintf("%8.3f",cc[1,indH1]),8), lpad(@sprintf("%8.3f",cc[2,indH1]),8), lpad(@sprintf("%8.3f",cc[3,indH1]),8),lpad("1.00",6),lpad("0.00",6), lpad("",10), lpad("H",2))
            println(f,"HETATM",lpad(indH2,5),"  ",rpad("H2",4),"HOH",lpad(i,6),"    ", lpad(@sprintf("%8.3f",cc[1,indH2]),8), lpad(@sprintf("%8.3f",cc[2,indH2]),8), lpad(@sprintf("%8.3f",cc[3,indH2]),8),lpad("1.00",6),lpad("0.00",6), lpad("",10), lpad("H",2))
            end
        println(f,"ENDMDL")
    end
    close(f)
end

function save_to_pdb(sr::SimulationResult, path )
    f = open(path,"w+")
    n = length(sr.simulation.system.bodies)
    L = 10*sr.simulation.boundary_conditions.L
    count = 0
    for t in sr.solution.t
        cc = 10*get_position(sr, t)
        map!(x ->  x -= L * floor(x / L), cc, cc)
        count+=1      
        println(f,rpad("MODEL",10), count)
        println(f,"REMARK 250 time=$t steps")
        for i ∈ 1:n
            
            println(f,"HETATM",lpad(i,5),"  ",rpad("Ar",4),"Ar",lpad(i,6),"    ", lpad(@sprintf("%8.3f",cc[1,i]),8), lpad(@sprintf("%8.3f",cc[2,i]),8), lpad(@sprintf("%8.3f",cc[3,i]),8),lpad("1.00",6),lpad("0.00",6), lpad("",10), lpad("Ar",2))
        end
        println(f,"ENDMDL")
    end
    close(f)
end