# f = H(p,q)
# Should create a conversion to a PartitionedODEProblem by autodifferentiation
type HamiltonianProblem{uType,tType,isinplace,F,C} <: AbstractPartitionedODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function HamiltonianProblem(f,u0,tspan; iip = isinplace(f,3),callback=nothing)
  @assert typeof(u0) <: Tuple
  @assert length(u0) == 2
  HamiltonianProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end

# Defined by (L,q') where q' is a vector for the q-derivatives
# Should create a conversion to a HamiltonianProblem
# H(p,q) = dot(p,q') - L
# Or a direct conversion to a PartitionedODEProblem
type LagrangianProblem{uType,tType,isinplace,F,C} <: AbstractPartitionedODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function LagrangianProblem(f,u0,tspan; iip = isinplace(f[2],4),callback=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(u0) <: Tuple
  @assert length(u0) == 2
  LagrangianProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end

function NBodyProblem(V,mass_matrix,u0,v0,tspan; kwargs...)
    @assert length(u0) == length(v0)
    @assert isa(V(tspan[1], u0), Real)
    NBodyProblem(Val{isdiag(mass_matrix)},V,mass_matrix,u0,v0,tspan; kwargs...)
end

function NBodyProblem(::Type{Val{true}},V,mass_matrix,u0,v0,tspan; kwargs...)
    function acceleration!(t, x, v, dv)
        ForwardDiff.gradient!(dv, q->V(t, q), x)
        for i in eachindex(dv)
            dv[i] /= -mass_matrix.diag[i]
        end
    end
    SecondOrderODEProblem(acceleration!, u0, v0, tspan; kwargs...)
end

function NBodyProblem(::Type{Val{false}},V,mass_matrix,u0,v0,tspan; kwargs...)
    factorized_M = cholfact(mass_matrix)
    function acceleration!(t, x, v, dv)
        ForwardDiff.gradient!(v, q->V(t, q), x)
        A_ldiv_B!(dv, factorized_M, v)
        scale!(dv, -1)
    end
    SecondOrderODEProblem(acceleration!, u0, v0, tspan; kwargs...)
end

