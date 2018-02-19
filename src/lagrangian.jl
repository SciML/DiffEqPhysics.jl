using OrdinaryDiffEq, ForwardDiff

# Auxiliary functions
∂L∂q(L, t, q, dq) = ForwardDiff.gradient(a->L(t,a,dq), q)
∂L∂dq(L, t, q, dq) = ForwardDiff.gradient(a->L(t,q,a), dq)
Dtq(L, t, q, dq) = ForwardDiff.derivative(a->∂L∂q(L,a,q,dq), t)
Dqdq(L, t, q, dq) = ForwardDiff.jacobian(a->∂L∂dq(L,t,a,dq), q)
Ddqdq(L, t, q, dq) = ForwardDiff.hessian(a->L(t,q,a), dq)
mass_mat(L, t, q, dq) = Ddqdq(L, t, q, dq)
function generalized_acceleration!(ddq, L, mass_matrix, t, q, dq)
    F   = ∂L∂q(L, t, q, dq)
    rhs = F - Dqdq(L, t, q, dq)*dq - Dtq(L, t, q, dq)
    A_ldiv_B!(ddq, mass_matrix, rhs)
end
function lagrangian(xv, K, V, charts!)
  (t, q, dq) -> begin
    charts!(xv, q, dq)
    x, v = xv
    K(x, v) - V(x, v)
  end
end

struct LagrangianProblem{iip} <: AbstractDynamicalODEProblem end

function LagrangianProblem(K,V,dim,charts!,q0,dq0,tspan,p=nothing; kwargs...)
  iip = (typeof(q0) <: AbstractArray) && !(typeof(q0) <: SArray)
  LagrangianProblem{iip}(K, V, dim, charts!, q0, dq0, tspan, p; kwargs...)
end

function LagrangianProblem{T}(K,V,dim,charts!,q0,dq0,tspan,p=nothing; kwargs...) where T
  caches = [similar(promote(q0, dq0)[1], dim) for i in 1:3]
  L = lagrangian(caches[1:2], K, V, charts!)
  mass_matrix = mass_mat(L, promote(tspan...)[1], q0, dq0)
  fact_mass_matrix = cholfact!(mass_matrix)
  f! = function (ddq, dq, q, p, t)
    generalized_acceleration!(L, fact_mass_matrix, t, q, dq)
  end
  xv0 = charts!(caches[3], q0, dq0)
  SecondOrderODEProblem{T}(f!, xv0[2], xv0[1], tspan; kwargs...)
end
