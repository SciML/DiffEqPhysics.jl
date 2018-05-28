abstract type NBodySystem
end

abstract type BasicPotentialSystem <: NBodySystem
end

struct ChargedParticles{bType <: ChargedParticle} <: BasicPotentialSystem
    bodies::Vector{bType}
    k::AbstractFloat
end

struct GravitationalSystem{bType <: MassBody} <: BasicPotentialSystem
    bodies::Vector{bType}
    G::Real
end

struct CustomAccelerationSystem{bType <: Body} <: NBodySystem
    bodies::Vector{bType}
    acceleration # f(u, v, i, system)
    parameters::Vector{<:Number}
end

struct PotentialNBodySystem{bType <: Body} <: NBodySystem
    bodies::Vector{bType}
    potentials::Dict{Symbol,<:PotentialParameters}
end

function PotentialNBodySystem(bodies::Vector{<:Body}; potentials::Vector{Symbol}=[])

    paramteres = Dict{Symbol,PotentialParameters}()
    
    if :lennard_jones ∈ potentials
        paramteres[:lennard_jones] = LennardJonesParameters()
    end

    if :electrostatic ∈ potentials 
        paramteres[:electrostatic] = ElectrostaticParameters()
    end

    if :gravitational ∈ potentials   
        paramteres[:gravitational] = GravitationalParameters()
    end

    PotentialNBodySystem(bodies, paramteres)
end

function Base.show(stream::IO, s::PotentialNBodySystem)
    print(stream, "Potentials: ")
    print(stream, "\n")
    for (key, potential) in s.potentials
        show(stream, potential)
    end
end

PotentialNBodySystem(system::PotentialNBodySystem) = system

function PotentialNBodySystem(system::ChargedParticles) 
    pp = ElectrostaticParameters(system.k)
    potential = Dict{Symbol,PotentialParameters}(:electrostatic => pp)
    PotentialNBodySystem(system.bodies, potential)
end

function PotentialNBodySystem(system::GravitationalSystem)
    pp = GravitationalParameters(system.G)
    potential =  Dict{Symbol,PotentialParameters}(:gravitational => pp)
    PotentialNBodySystem(system.bodies, potential)
end