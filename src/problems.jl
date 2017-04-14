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
