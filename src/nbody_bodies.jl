# Fields for position, velocity and mass of a particle should be required for every descendant
abstract type Body
end

struct MassBody{mType <: AbstractFloat,cType <: AbstractFloat} <: Body
    r::SVector{3,cType} #required 
    v::SVector{3,cType} #required
    m::mType             #required
end

struct ChargedParticle{fType <: AbstractFloat} <: Body
    r::SVector{3,fType}
    v::SVector{3,fType}
    m::fType 
    q::fType     
end

struct SphericalBody{mType,qType,cType <: AbstractFloat} <: Body
    r::SVector{3,cType}  
    v::SVector{3,cType}  
    m::mType              
    R::cType
    q::qType 
    mm::SVector{3,qType}
end