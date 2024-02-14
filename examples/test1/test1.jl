using Gridap, SemismoothQVIs
using Plots, LaTeXStrings


n = 400 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = FEM_1D_model(n) # Piecewise linear FEM discretization on (0,1)

# Thermoforming parameter choices
k = 1
Φ₀ = interpolate_everywhere(x->max(0, abs(x[1]-0.5)-0.25), VT)
Ψ₀ = interpolate_everywhere(x->1.0, VT)
α₁ = 1.0
ϕ = interpolate_everywhere(x->sin(π*x[1])/α₁, UT)
ψ = ϕ
f = interpolate_everywhere(x->π^2*sin(π*x[1]) + max(0, -abs(x[1] - 0.625) + 0.125), Uu)
α₂ = 1.0
g(s) = (1-s) ≤ 0.0 ? α₁ :  α₁ + (1-s)/α₂
dg(s)= (1-s) ≤ 0.0 ? 0.0 :  -1/α₂

# Explicit solution
(u₁, T₁) = (interpolate_everywhere(x->sin(π*x[1]), Vu), interpolate_everywhere(x->α₁, VT))

# Thermoforming struct
Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)

# Initial guess
u₀ = interpolate_everywhere(x->0.0, Uu)
T₀ = interpolate_everywhere(x->5.0, UT)

# fixed point iteration
(zhs1, h1_1, its_1) = fixed_point(Q, u₀, T₀; max_its=10, tol=1e-12,  ρ0=1, PF=false)

# Semismooth Newton method
(zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, u₀, T₀; max_its=10, tol=1e-12, PF=false, globalization=true);

# Extract final solution
uh = zhs1[its_1[1]][1]
Th = zhs1[its_1[1]][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
# Check error
(h1(Q, u₁, uh), h1(Q, T₁, Th))

# uh = zhs3[its_3[1]+1][1]
# Th = zhs3[its_3[1]+1][2]
# mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

# (h1(Q, u₁, uh), h1(Q, T₁, Th))

# Plot solution
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mold(Point.(xx))],#  mold(Point.(xx)) Th(Point.(xx))],
    label=["Membrane" "Mold"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=2, xlabelfontsize=20)