abstract type PotentialParameters
end

struct LennardJonesParameters <: PotentialParameters
    ϵ::AbstractFloat
    σ::AbstractFloat
    R::AbstractFloat
    σ2::AbstractFloat
    R2::AbstractFloat
end

function LennardJonesParameters(ϵ::AbstractFloat, σ::AbstractFloat, R::AbstractFloat)
    LennardJonesParameters(ϵ, σ, R, σ^2, R^2)
end

function LennardJonesParameters()
    LennardJonesParameters(1.0, 1.0, 2.5)
end

function Base.show(stream::IO, pp::LennardJonesParameters)
    print(stream, "Lennard-Jones:\n")
    print(stream, "\tϵ:"); show(stream, pp.ϵ); print(stream, "\n")
    print(stream, "\tσ:"); show(stream, pp.σ); print(stream, "\n")
    print(stream, "\tR:"); show(stream, pp.R); print(stream, "\n")
end

struct GravitationalParameters <: PotentialParameters
    G::AbstractFloat
end

function GravitationalParameters()
    GravitationalParameters(6.67408e-11)
end

function Base.show(stream::IO, pp::GravitationalParameters)
    print(stream, "Gravitational:\n")
    print(stream, "\tG:"); show(stream, pp.G); print(stream, "\n")
end

struct ElectrostaticParameters <: PotentialParameters
    k::AbstractFloat
end

function ElectrostaticParameters()
    ElectrostaticParameters(9e9)
end

function Base.show(stream::IO, pp::ElectrostaticParameters)
    print(stream, "Electrostatic:\n")
    print(stream, "\tk:"); show(stream, pp.k); print(stream, "\n")
end