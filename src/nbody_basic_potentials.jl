abstract type PotentialParameters
end

struct LennardJonesParameters{pType <: Real} <: PotentialParameters
    ϵ::pType
    σ::pType
    R::pType
    σ2::pType
    R2::pType
end

function LennardJonesParameters(ϵ::Real, σ::Real, R::Real)
    LennardJonesParameters(ϵ, σ, R, σ^2, R^2)
end

function LennardJonesParameters()
    LennardJonesParameters(1.0, 1.0, 2.5)
end

function Base.show(stream::IO, pp::LennardJonesParameters)
    println(stream, "Lennard-Jones:")
    print(stream, "\tϵ:"); show(stream, pp.ϵ); println(stream)
    print(stream, "\tσ:"); show(stream, pp.σ); println(stream)
    print(stream, "\tR:"); show(stream, pp.R); println(stream)
end

struct GravitationalParameters{gType <: Real} <: PotentialParameters
    G::gType
end

function GravitationalParameters()
    GravitationalParameters(6.67408e-11)
end

function Base.show(stream::IO, pp::GravitationalParameters)
    println(stream, "Gravitational:")
    print(stream, "\tG:"); show(stream, pp.G); 
    println(stream)
end

struct ElectrostaticParameters{kType <: Real} <: PotentialParameters
    k::kType
end

function ElectrostaticParameters()
    ElectrostaticParameters(9e9)
end

function Base.show(stream::IO, pp::ElectrostaticParameters)
    println(stream, "Electrostatic:")
    print(stream, "\tk:"); show(stream, pp.k); 
    println(stream)
end

struct MagnetostaticParameters{mType <: Real} <: PotentialParameters
    μ_4π::mType
end

function MagnetostaticParameters()
    MagnetostaticParameters(1e-7)
end

function Base.show(stream::IO, pp::MagnetostaticParameters)
    println(stream, "Magnetostatic:")
    print(stream, "\tμ/4π:"); show(stream, pp.μ_4π); 
    println(stream)
end

struct SPCFwParameters{pType <: Real} <: PotentialParameters
    rOH::pType
    ∠HOH::pType
    kb::pType
    ka::pType
end