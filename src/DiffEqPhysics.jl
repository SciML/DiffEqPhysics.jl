module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, RecursiveArrayTools
using ForwardDiff, StaticArrays, RecipesBase
using Random, Printf, LinearAlgebra

include("nbody_problem.jl")
include("hamiltonian.jl")
include("plot.jl")

export HamiltonianProblem, NBodyProblem, plot_orbits

export NBodySimulation
export MassBody, ChargedParticle, MagneticParticle
export PotentialParameters, LennardJonesParameters, GravitationalParameters,
       ElectrostaticParameters, MagnetostaticParameters, SPCFwParameters
export PotentialNBodySystem, ChargedParticles, GravitationalSystem, WaterSPCFw
export PeriodicBoundaryConditions, CubicPeriodicBoundaryConditions, InfiniteBox
export AndersenThermostat, BerendsenThermostat, NoseHooverThermostat, LangevinThermostat
export run_simulation, get_position, get_velocity, get_masses, temperature,
       initial_energy, kinetic_energy, potential_energy, total_energy, rdf, msd,
       generate_bodies_in_cell_nodes, run_simulation_sde

end # module
