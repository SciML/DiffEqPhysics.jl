abstract type BoundaryConditions
end

struct PeriodicBoundaryConditions{cType <: Real} <: BoundaryConditions
    boundary::SVector{6,cType}
end

PeriodicBoundaryConditions(L::Real) = PeriodicBoundaryConditions(SVector(0, L, 0, L, 0, L))

Base.start(::PeriodicBoundaryConditions) = 1

Base.done(pbc::PeriodicBoundaryConditions, state) = state > length(pbc.boundary)

function Base.next(pbc::PeriodicBoundaryConditions, state) 
    pbc.boundary[state], state + 1
end

function Base.getindex(pbc::PeriodicBoundaryConditions, i::Integer)
    1 <= i <= length(pbc.boundary) || throw(BoundsError(pbc, i))
    pbc.boundary[i]
end


struct InfiniteBox{cType <: Real} <: BoundaryConditions
    boundary::SVector{6,<:cType}
end

InfiniteBox() = InfiniteBox(SVector(-Inf, Inf, -Inf, Inf, -Inf, Inf))