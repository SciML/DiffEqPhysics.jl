__precompile__()

module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq, RecursiveArrayTools
using ForwardDiff, StaticArrays, RecipesBase

include("nbody_problem.jl")
include("hamiltonian.jl")
include("plot.jl")
include("nbody_simulation.jl")

export HamiltonianProblem, LagrangianProblem, NBodyProblem, plot_orbits

export NBodySimulation
export MassBody, ChargedParticle, MagneticParticle
export PotentialParameters, LennardJonesParameters, GravitationalParameters, 
       ElectrostaticParameters, MagnetostaticParameters, SPCFwParameters
export PotentialNBodySystem, ChargedParticles, GravitationalSystem, WaterSPCFw
export PeriodicBoundaryConditions, CubicPeriodicBoundaryConditions, InfiniteBox
export AndersenThermostat, BerendsenThermostat
export run_simulation, get_position, get_velocity, get_masses, temperature,
       initial_energy, kinetic_energy, potential_energy, total_energy, rdf, msd,
       generate_bodies_in_cell_nodes

end # module
