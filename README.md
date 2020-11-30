# DiffEqPhysics

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/SciML/DiffEqPhysics.jl/workflows/CI/badge.svg)](https://github.com/SciML/DiffEqPhysics.jl/actions?query=workflow%3ACI)

This package provides `HamiltonianProblem`. For `NBodyProblem`, please use package [NBodySimulator](https://github.com/SciML/NBodySimulator.jl)

<!--
## Simulation of gravitationally interacting bodies

In order to create bodies/particles for the problem, one needs to use the MassBody structure and its constructor, which accepts mass, initial coordinates and velocity of the body.

```julia
body1 = MassBody(2.0, SVector(0.0, 1.0, 0.0), SVector( 5.775e-6, 0.0, 0.0))
body2 = MassBody(2.0, SVector(0.0,-1.0, 0.0), SVector(-5.775e-6, 0.0, 0.0))
```

Usually we solve an n-body problem for a certain period of time:

```julia
tspan = (0.0, 1111150.0);
```

Solving gravitational problem one needs to specify the gravitational constant G.
```julia
G = 6.673e-11
```

In fact, now we have enough parameters to create an NBodyGravProblem object:

```julia
problem = NBodyGravProblem([body1,body2], G, tspan)
```

Solution to the problem might be evaluated using the standard `solve` function:
```julia
solution = solve(problem, Tsit5());
```

And, finally, we plot our solution showing two equal bodies rotating on the same orbit:
```julia
plot_xy_scattering(solution,"./anim_two_boddies_scattering.gif")
```

<img src="https://user-images.githubusercontent.com/16945627/39958539-d2cf779c-561d-11e8-96a8-ffc3a595be8b.gif" alt="Here should appear a gif of rotating bodies" width="350"/>

-->
