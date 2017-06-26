type AbstractPartitionedODEProblem  end
# f = H(p,q)
# Should create a conversion to a PartitionedODEProblem by autodifferentiation
type HamiltonianProblem{uType,tType,isinplace,F,C} <: AbstractPartitionedODEProblem{uType,tType,isinplace}
  f :: F
  u0:: uType
  tspan   :: Tuple{tType,tType}
  callback:: C
end

function HamiltonianProblem(f,u0,tspan; iip = isinplace(f,3),callback=nothing)
  @assert typeof(f)  <: Tuple
  @assert typeof(u0) <: Tuple
  @assert length(u0) == 2
  HamiltonianProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end

# Defined by (L,q') where q' is a vector for the q-derivatives
# Should create a conversion to a HamiltonianProblem
# H(p,q) = dot(p,q') - L
# Or a direct conversion to a PartitionedODEProblem
type LagrangianProblem{uType,tType,isinplace,F,C} <: AbstractPartitionedODEProblem{uType,tType,isinplace}
  f :: F
  u0:: uType
  tspan :: Tuple{tType,tType}
  callback :: C
end

function LagrangianProblem(f,u0,tspan; iip = isinplace(f[2],4),callback=nothing)
  @assert typeof(f)  <: Tuple
  @assert typeof(u0) <: Tuple
  @assert length(u0) == 2
  LagrangianProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end

# NBody Problem

type NBodyProblem{uType,tType,isinplace,F,C,NBodies} <: AbstractPartitionedODEProblem{uType,tType,isinplace}
  f  :: F
  u0 :: uType
  tspan    :: Tuple{tType, tType}
  callback :: C
end
function NBodyProblem{f,u0,tspan; iip = isinplace(f[2],4), callback=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(u0) <: Tuple
  @assert isa(Nbodies,Int)
  @assert length(u0) == NBodies*4
  NBodyProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
  end
u = [x;y;v;w]

pleiades = (t,u,du) -> begin
  x = view(u,1:NBodies)   # x
  y = view(u,NBodies+1:2*NBodies)  # y
  v = view(u,2*NBodies+1:3*NBodies) # x′
  w = view(u,3*NBodies+1:4*NBodies) # y′
  du[1:NBodies] .= v
  du[NBodies+1:2*NBodies].= w
  for i in 2*NBodies:3*NBodies
    du[i] = zero(u[1])
  end
  for i=1:NBodies,j=1:NBodies
    if i != j
      r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
      du[2*NBodies+i] += j*(x[j] - x[i])/r
      du[3*NBodies+i] += j*(y[j] - y[i])/r
    end
  end
end
tspan= (0.,3.)

prob = ODEProblem(pleiades,u,tspan)

sol = solve(prob)

plot(sol, vars=[ (i, i+NBodies) for i in 1:NBodies] )
