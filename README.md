# SemismoothQVIs

This repository implements the numerical examples found in:

(1) "A Globalized Inexact Semismooth Newton Method for Obstacle-type Quasi-variational Inequalities with Applications in Thermoforming", Amal Alphonse, Constantin Christof, Ioannis. P. A. Papadopoulos and Michael Hintermüller (2024).

## Installation
To install the package, try:

```julia
julia> ] add git@github.com:ioannisPApapadopoulos/SemismoothQVIs.jl.git
```

If that does not work try:

```julia
julia> ] add https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl.git
```

If that also does not work, contact John. 

Last resort: download the zip file, unzip, locate the folder and run

```julia
julia> ] 
pkg> activate /location to folder/SemismoothQVIs
  Activating project at `\location to folder\SemismoothQVIs.jl`
(SemismoothQVIs) pkg> initialize
```
If it says "Activating new project" then you are in the wrong folder.

Checkout the examples in the examples/ folder and have fun!

## Solvers

The solvers implemented are `fixed_point`, `semismoothnewton` &  `moreau_yosida_newton`.

`fixed_point` employs the algorithm uᵢ₊₁ = S(Φ(uᵢ)). This converges provided there is a contraction.
`semismoothnewton` is Algorithm 1 in the paper.
`moreau_yosida_newton` is a smoothed regularization of Algorithm 1.

The optional arguments in the solvers are:

|Optional arguments|Type|Description|Notes|
|:-:|:-:|:-:|:-:|
|`max_its`|Integer|Maximum number of iterations of the outer solver before it is terminated.|Default is `20`.|
|`tol`|Float|Absolute tolerance for inner and outer solvers|Default is `1e-10`.|
|`globalization`|Boolean (true/false)|Compute globalization iterate $x_B$ as well as the Newton iterate $x_N$|Default is `false`. Only for `semismoothnewton` and `moreau_yosida_newton`.|
|proj_rc|(Float, Float)|(radius, center) for projection.|Default is `(Inf, 0.0)`. Not yet implemented for `moreau_yosida_newton`. Nonzero center not yet implemented for `semismoothnewton`.|
|`bt`|Boolean (true/false)|Backtracking linesearch for Newton solver.|Default is `true`.|
|`PF`|Boolean (true/false)|Use a path-following Moreau-Yosida regularization to solve the inner obstacle problem.|Default is `true`. Only for `fixed_point` and `semismoothnewton`.|
|`ρ0`|Float|Initial parameter for the path-following Moreau-Yosida regularization.|Default is `1`. Only for `fixed_point`.|
|`FS`|Boolean (true/false)|Use the primal-dual active set strategy (possibly after a path-following).|Default is `true`. Only for `fixed_point` and `semismoothnewton`.|
|`ρ`|Float|Fixed Moreau-Yosida parameter for `moreau_yosida_newton`.|Default is `1e-5`. Only for `moreau_yosida_newton`.|
|`inner_max_its`|Integer|Max number of iterations for obstacle problem subsolver in `moreau_yosida_newton`.|Default is `400`. Only for `moreau_yosida_newton`.|
|`show_inner_trace`|Boolean (true/false)|Show the trace of the inner subsolvers.|Default is `true`.|
|`show_trace`|Boolean (true/false)|Show the trace of the outer solve.|Default is `true`.|


