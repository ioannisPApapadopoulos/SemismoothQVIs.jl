using Gridap, SemismoothQVIs
using Plots, LaTeXStrings


n = 100 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = FEM_1D_model(n) # Piecewise linear FEM discretization on (0,1)

# Thermoforming parameter choices
k = π^2
Φ₀ = interpolate_everywhere(x->0.0, VT)
Ψ₀ = interpolate_everywhere(x->0.0, VT)
α₁ = 1.0
α₂ = 1.2
ϕ = interpolate_everywhere(x->α₂*10*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
ψ = interpolate_everywhere(x->5*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
f = interpolate_everywhere(x->α₁*π^2*sin(π*x[1]), Uu)

g(s) = s ≥ 0.0 ? 0.0 : 4 / α₁ * s^2
dg(s)= s ≥ 0.0 ? 0.0 : 8 / α₁ * s

# Explicit solutions
(u₁, T₁) = (interpolate_everywhere(x->0.0, Vu), interpolate_everywhere(x->0.0, VT))
(u₂, T₂) = (interpolate_everywhere(x->α₁*sin(π*x[1]), Vu), interpolate_everywhere(x->α₁*(5-cos(2π*x[1]))/( 10*π^2 ), VT))

# Thermoforming struct
Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)

# Initial guess
T₀ = interpolate_everywhere(x->5.0, UT)
# Initial guess of u₀ = (-Δ)⁻¹ f
J = Gridap.Algebra.jacobian(FEOperator(Q.au0..., Uu, Vu), interpolate_everywhere(x->0.0, Uu))
u₀ = FEFunction(Vu, J \ f.free_values)


# fixed point iteration
(zhs1, h1_1, its_1) = fixed_point(Q, u₀, T₀; max_its=100, tol=1e-20,  ρ0=1, PF=false, proj_rc=(1.2, 0.0))

# Semismooth Newton method
(zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, u₀, T₀; max_its=10, tol=1e-20, PF=false, globalization=true, proj_rc=(1.2,0.0));

# Extract final solution
uh = zhs1[its_1[1]][1]
Th = zhs1[its_1[1]][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
# Check errors
(h1(Q, u₁, uh), h1(Q, T₁, Th))
(h1(Q, u₂, uh), h1(Q, T₂, Th))

uh = zhs3[its_3[1]+1][1]
Th = zhs3[its_3[1]+1][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
(h1(Q, u₁, uh), h1(Q, T₁, Th))
(h1(Q, u₂, uh), h1(Q, T₂, Th))

# Plot solution
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mold(Point.(xx))],#  mold(Point.(xx)) Th(Point.(xx))],
    label=["Membrane" "Mold"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=2, xlabelfontsize=20)


ls = [:solid, :dash, :dashdot, :dashdotdot]
p = Plots.plot(1:length(h1_1), h1_1, linestyle=ls[1], linewidth=2, marker=:dot, label="Fixed Point", yscale=:log10,)
Plots.plot!(1:length(h1_3), h1_3,
    linewidth=2,
    linestyle=ls[1],
    marker=:square,
    yscale=:log10,
    label="Semismooth Newton",
    xlabel="Iterations",
    ylabel=L"\Vert R(u_{i}) \Vert_{H^1}",
    yticks=[1e-15,1e-10,1e-5,1e0],
    xtick=2:2:22,
    ylim=[1e-18,1e3]
)