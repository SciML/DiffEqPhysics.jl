# This is an attempt to include the required for n-body simulations fields into structures
# Indeed, tools of ConcreteAbstractions.jl seem to be more suitable for the fields inhereting
DiffEqBase.@def position_velocity_mass begin
    r::SVector{3,cType}
    v::SVector{3,cType}
    m::mType
end

abstract type Body
end

struct MassBody{cType <: Real,mType <: Real} <: Body
    @position_velocity_mass
end

struct ChargedParticle{cType <: Real,mType <: Real,qType <: Real} <: Body
    @position_velocity_mass
    q::qType
end

struct MagneticParticle{cType <: Real,mType <: Real,mmType <: Real} <: Body
    @position_velocity_mass
    mm::SVector{3,mmType}
end
