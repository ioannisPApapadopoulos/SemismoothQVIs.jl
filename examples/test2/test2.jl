using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

n = 200 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = fem_model(n) # Piecewise linear FEM discretization on (0,1)

# Thermoforming parameter choices
k = π^2
Φ₀ = interpolate_everywhere(x->0.0, VT)
Ψ₀ = interpolate_everywhere(x->0.0, VT)
α₁ = 100.0
α₂ = 1.1
ϕ = interpolate_everywhere(x->α₂*10*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
ψ = interpolate_everywhere(x->5*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
f = interpolate_everywhere(x->α₁*π^2*sin(π*x[1]), Uu)

# Radius of (u,T) = (0,0) solution
Rα = 3α₁*  ( sqrt( ((5α₂ + 8)*π^2 + 8π) / (80α₂) )  - π/4   ) /  (10*(1+π)) 
50α₂/3/α₁ *(Rα + 20*(1+π)/3/π/α₁ * Rα^2) ≈ 1

g(s) = s ≥ 0.0 ? 0.0 : 4 / α₁ * s^2
dg(s)= s ≥ 0.0 ? 0.0 : 8 / α₁ * s

# Explicit solutions
(u₁, T₁) = (interpolate_everywhere(x->0.0, Vu), interpolate_everywhere(x->0.0, VT))
(u₂, T₂) = (interpolate_everywhere(x->α₁*sin(π*x[1]), Vu), interpolate_everywhere(x->α₁*(5-cos(2π*x[1]))/( 10*π^2 ), VT))

# Thermoforming struct
Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)

# Initial guess
T₀ = interpolate_everywhere(x->5.0, UT)
u₀ = interpolate_everywhere(x->400*(1.0 -x[1])*x[1], Uu)

# Solve via fixed point & semismooth Newton, both converge to maximal solution
(zhs1_max, h1_1, its_1_max) = fixed_point(Q, u₀, T₀; max_its=50, tol=1e-20, PF=false, bt=true, proj_rc=(Inf, 0.0), show_inner_trace=false);
(zhs3_max, h1_3, its_3_max, is_3) = semismoothnewton(Q, u₀, T₀; max_its=100, tol=1e-20, PF=false, globalization=true, proj_rc=(Inf,0.0), show_inner_trace=false);

# Solve via fixed point & semismooth Newton + projection, both converge to minimal solution
(zhs1_min, h1_1, its_1_min) = fixed_point(Q, u₀, T₀; max_its=50, tol=1e-20, PF=false, bt=true, proj_rc=(Rα/2, 0.0), show_inner_trace=false);
(zhs3_min, h1_3, its_3_min, is_3) = semismoothnewton(Q, u₀, T₀; max_its=100, tol=1e-20, PF=false, globalization=false, proj_rc=(Rα/2,0.0), show_inner_trace=false);



uh = zhs1_max[its_1_max[1]][1]; Th = zhs1_max[its_1_max[1]][2]; (h1(Q, u₂, uh), h1(Q, T₂, Th))
uh = zhs1_min[its_1_min[1]][1]; Th = zhs1_min[its_1_min[1]][2]; (h1(Q, u₁, uh), h1(Q, T₁, Th))

uh = zhs3_max[its_3_max[1]+1][1]; Th = zhs3_max[its_3_max[1]+1][2]; (h1(Q, u₂, uh), h1(Q, T₂, Th))
uh = zhs3_min[its_3_min[1]+1][1]; Th = zhs3_min[its_3_min[1]+1][2]; (h1(Q, u₁, uh), h1(Q, T₁, Th))


mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

# Plot solution
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mold(Point.(xx))],#  mold(Point.(xx)) Th(Point.(xx))],
    label=["Membrane" "Mold"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=2, xlabelfontsize=20)


# ls = [:solid, :dash, :dashdot, :dashdotdot]
# p = Plots.plot(1:length(h1_1), h1_1, linestyle=ls[1], linewidth=2, marker=:dot, label="Fixed Point", yscale=:log10,)
# Plots.plot!(1:length(h1_3), h1_3,
#     linewidth=2,
#     linestyle=ls[1],
#     marker=:square,
#     yscale=:log10,
#     label="Semismooth Newton",
#     xlabel="Iterations",
#     ylabel=L"\Vert R(u_{i}) \Vert_{H^1}",
#     yticks=[1e-15,1e-10,1e-5,1e0],
#     xtick=2:2:22,
#     ylim=[1e-18,1e3]
# )