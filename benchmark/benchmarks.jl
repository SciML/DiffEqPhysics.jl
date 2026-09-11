using DiffEqPhysics, OrdinaryDiffEq, OrdinaryDiffEqSymplecticRK, BenchmarkTools

const SUITE = BenchmarkGroup()

# Gravitational two-body Hamiltonian
function H(p, q, params)
    return sum(abs2, p) / 2 - 1 / sqrt(sum(abs2, q))
end
p0 = [0.0, 0.5]
q0 = [1.0, 0.0]
prob = HamiltonianProblem(H, p0, q0, (0.0, 100.0))

# Harmonic oscillator
function H_ho(p, q, params)
    return (sum(abs2, p) + sum(abs2, q)) / 2
end
prob_ho = HamiltonianProblem(H_ho, [0.0], [1.0], (0.0, 500.0))

# =============================================================================
# Hamiltonian solves (symplectic + adaptive)
# =============================================================================

SUITE["solve"] = BenchmarkGroup()

SUITE["solve"]["kepler_symplectic"] = @benchmarkable solve(
    $prob, VelocityVerlet(); dt = 0.01
)
SUITE["solve"]["kepler_rk"] = @benchmarkable solve($prob, Tsit5())
SUITE["solve"]["oscillator_symplectic"] = @benchmarkable solve(
    $prob_ho, SymplecticEuler(); dt = 0.05
)

# =============================================================================
# Construction
# =============================================================================

SUITE["construct"] = BenchmarkGroup()
SUITE["construct"]["hamiltonian_problem"] = @benchmarkable HamiltonianProblem(
    $H, $p0, $q0, (0.0, 100.0)
)
