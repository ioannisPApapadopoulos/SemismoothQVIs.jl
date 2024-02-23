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

The solvers implemented are fixed_point, semismoothnewton &  moreau_yosida_newton.

fixed_point employs the algorithm uᵢ₊₁ = S(Φ(uᵢ)). This converges provided there is a contraction.
semismoothnewton is Algorithm 1.1 in the paper.
moreau_yosida_newton is a smoothed regularization of Algorithm 1.1/

The optional arguments in the solvers are:

|Optional arguments|Type|Description|
|:-:|:-:|:-:|
|max_its|Integer|Maximum number of iterations of the outer solver before it is terminated.||
|tol|Float|.|
|globalization|Boolean (true/false)|Compute globalization iterate x_B as well as the Newton iterate x_N|Only for semismoothnewton and moreay_yosida_newton.|