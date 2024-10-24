# SemismoothQVIs

This repository implements the numerical examples found in:

"A globalized inexact semismooth Newton method for nonsmooth fixed-point equations involving variational inequalities", Amal Alphonse, Constantin Christof, Michael Hintermüller, and Ioannis. P. A. Papadopoulos (2024) [arXiv](https://arxiv.org/abs/2409.19637).

|Figure|File: examples/|
|:-:|:-:|
|1a, 1d|[test1/test1a.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test1/test1a.jl)|
|1b, 1c|[test1/test1b.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test1/test1b.jl)|
|2|[test2/test2a.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test2/test2a.jl)|
|3|[test3/test3.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test3/test3.jl)|
|4|[test4/test4.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test4/test4.jl)|

|Table|File: examples/|
|:-:|:-:|
|1|[test2/test2b.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test2/test2b.jl)|
|2|[test3/test3.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test3/test3.jl)|
|3|[test4/test4.jl](https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl/tree/main/examples/test4/test4.jl)|

## Installation

The package is not registered. Please install via

```pkg> add https://github.com/ioannisPApapadopoulos/SemismoothQVIs.jl.git```

## Solvers

The solvers implemented are `fixed_point`, `semismoothnewton` &  `moreau_yosida_newton`.

`fixed_point` employs the algorithm uᵢ₊₁ = S(Φ(uᵢ)). This converges provided there is a contraction.
`semismoothnewton` implements Algorithms 1 and 2 in the paper.
`moreau_yosida_newton` is a smoothed regularization of Algorithm 1.

The optional arguments in the solvers are:

|Optional arguments|Type|Description|Notes|
|:-:|:-:|:-:|:-:|
|`max_its`|Integer|Maximum number of iterations of the outer solver before it is terminated.|Default is `20`.|
|`out_tol`|Float|Absolute tolerance for outer solver|Default is `1e-10`.|
|`in_tol`|Float|Absolute tolerance for inner solvers|Default is `1e-10`.|
|`hik_tol`|Float|Absolute tolerance for inner HIK solver|Default is `1e-10`.|
|`linesearch`|Boolean (true/false)|Perform a backtracking linesearch for the SSN update|Default is `false`. Only for `semismoothnewton`.|
|`globalization`|Boolean (true/false)|Compute globalization iterate $x_B$ as well as the Newton iterate $x_N$|Default is `false`. Only for `semismoothnewton` and `moreau_yosida_newton`.|
|`proj_rc`|(Float, Float)|(radius, center) for projection.|Default is `(Inf, 0.0)`. Not yet implemented for `moreau_yosida_newton`. Nonzero center not yet implemented for `semismoothnewton`.|
|`bt`|Boolean (true/false)|Backtracking linesearch for (inner) Newton solver.|Default is `true`.|
|`PF`|Boolean (true/false)|Use a path-following Moreau-Yosida regularization to solve the inner obstacle problem.|Default is `true`. Only for `fixed_point` and `semismoothnewton`.|
|`ρ0`|Float|Initial parameter for the path-following Moreau-Yosida regularization.|Default is `1`. Only for `fixed_point`.|
|`ρ_min`|Float|Terminate path-following Moreau-Yosida regularization if ρ<ρ_min.|Default is `1e-6`. Only for `fixed_point` and `semismoothnewton`.|
|`FS`|Boolean (true/false)|Use the primal-dual active set strategy (possibly after a path-following).|Default is `true`. Only for `fixed_point` and `semismoothnewton`.|
|`ρ`|Float|Fixed Moreau-Yosida parameter for `moreau_yosida_newton`.|Default is `1e-5`. Only for `moreau_yosida_newton`.|
|`inner_max_its`|Integer|Max number of iterations for obstacle problem subsolver in `moreau_yosida_newton`.|Default is `400`. Only for `moreau_yosida_newton`.|
|`show_inner_trace`|Boolean (true/false)|Show the trace of the inner subsolvers.|Default is `true`.|
|`show_trace`|Boolean (true/false)|Show the trace of the outer solve.|Default is `true`.|


