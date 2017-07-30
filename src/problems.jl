# f = H(p,q)
# Should create a conversion to a PartitionedODEProblem by autodifferentiation
type HamiltonianProblem{uType,tType,isinplace,F,C} <: AbstractODEProblem{uType,tType,isinplace}
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
type LagrangianProblem{uType,tType,isinplace,F,C} <: AbstractODEProblem{uType,tType,isinplace}
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

function NBodyProblem(potential,mass_matrix,u0,v0,tspan; kwargs...)
    # check number of particles
    @assert length(u0) == length(v0)
    @assert isdiag(mass_matrix)
    ind = i -> div(i-1, dim) + 1
    N   = length(u0)
    dim = div(N, size(mass_matrix, 1))
    @assert dim*size(mass_matrix, 1) == N
    @assert dim ∈ (2, 3)

    function f_wrapper2(t,u)
        x = @view u[1:2:end]
        y = @view u[2:2:end]
        potential(t,x,y,mass_matrix)
    end

    function f_wrapper3(t,u)
        x = @view u[1:3:end]
        y = @view u[2:3:end]
        z = @view u[3:3:end]
        potential(t,x,y,z,mass_matrix)
    end

    function acceleration!(t, x, v, dv)
        N = length(v)
        config_length = N > 10 ? 10 : N
        fun = dim == 2 ? q -> f_wrapper2(t, q) : q -> f_wrapper3(t, q)
        cfg = ForwardDiff.GradientConfig(fun, u0, ForwardDiff.Chunk{config_length}())
        ForwardDiff.gradient!(dv, fun, x, cfg)
        for i in eachindex(dv)
            dv[i] /= -mass_matrix.diag[ind(i)]
        end
    end
    SecondOrderODEProblem(acceleration!, u0, v0, tspan; kwargs...)
end

