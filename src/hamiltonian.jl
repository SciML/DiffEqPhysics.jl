using OrdinaryDiffEq, ForwardDiff

function HamiltonianProblem(
  H::F,q0::Real,p0::Real,tspan::Tuple{T,T}; kwargs...) where T<:Real

  D = length(q0)
  Ham = (x) -> H(x[1], x[2])
  function EOM(Ham)
    grad!(result, x) = ForwardDiff.gradient!(result, Ham, x::AbstractArray)
    #magic
  return dHdq, dHdp

  dHdq, dHdp = EOM(HAM)
  return ODEProblem((dHdp,-dHdq), (q0,p0), tspan; kwargs...)
end

function HamiltonianProblem(
  H::F,q0::AbstractVector,p0::AbstractVector,tspan::Tuple{T,T}; kwargs...) where T<:Real

  D = length(q0)
  Ham = (x) -> H(x[1:D], x[D+1:end])
  function EOM(Ham, D)
    grad!(result, x) = ForwardDiff.gradient!(result, Ham, x::AbstractArray)
    #magic
  return dHdq, dHdp

  dHdq, dHdp = EOM(HAM)
  return ODEProblem((dHdp,-dHdq), (q0,p0), tspan; kwargs...)
end
