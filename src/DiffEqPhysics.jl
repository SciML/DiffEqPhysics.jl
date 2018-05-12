__precompile__()

module DiffEqPhysics

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq, RecursiveArrayTools
using ForwardDiff, StaticArrays, RecipesBase

include("nbody.jl")
include("hamiltonian.jl")
include("plot.jl")
include("nbody_gravitational")
include("animation_gravitational")

export HamiltonianProblem, LagrangianProblem, NBodyProblem, plot_orbits

export NBodyGravProblem, MassBody, GravImagingData, solve, 
	   plot_xy_trailing, plot_xy_scattering, plot_xy


end # module
