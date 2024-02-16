using Gridap, SemismoothQVIs
using Plots, LaTeXStrings


n = 200 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = FEM_1D_model(n) # Piecewise linear FEM discretization on (0,1)

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
# Initial guess of u₀ = (-Δ)⁻¹ f
J = Gridap.Algebra.jacobian(FEOperator(Q.au0..., Uu, Vu), interpolate_everywhere(x->0.0, Uu))
u₀ = FEFunction(Vu, J \ f.free_values)
u₀ = interpolate_everywhere(x->400*(1.0 -x[1])*x[1], Uu)

# fixed point iteration
aT = Q.aT[1]
xhs = []
is1s = []
is2s = []
hus = []
hTs = []
hms = []
RT = []
h1_norm(u) = sqrt(sum(∫(∇(u) ⋅ ∇(u) + u ⋅ u)*Q.dΩ))
# for m = [0.1, 1.0, 5.0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, Inf]
# for m = [100.0]
m = 9
    (zhs1, h1_1, its_1) = fixed_point(Q, u₀, T₀; max_its=50, tol=1e-20, PF=false, bt=true, proj_rc=(Inf, 0.0), show_inner_trace=false);
    append!(xhs, [zhs1[its_1[1]]])
    uh = zhs1[its_1[1]][1]
    Th = zhs1[its_1[1]][2]
    mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
    [h1(Q, u₁, uh), h1(Q, T₁, Th)]
    [h1(Q, u₂, uh), h1(Q, T₂, Th)]
    # append!(is1s, norm( [h1(Q, u₁, uh), h1(Q, T₁, Th)], Inf) < 1e-10)
    # append!(is2s, norm( [h1(Q, u₂, uh), h1(Q, T₂, Th)], Inf) < 1e-1)
    append!(hus, h1_norm(uh))
    append!(hTs, h1_norm(Th))
    append!(hms, h1_norm(mold))

    bT(Th, R) = aT(Th, R, uh)
    append!(RT, norm(Gridap.Algebra.residual(FEOperator(bT, UT, VT), Th)))

end

uh, Th = u₂, T₂
aT = Q.aT[1]
bT(Th, R) = aT(Th, R, uh)
norm(Gridap.Algebra.residual(FEOperator(bT, UT, VT), Th))

au0, ju0 = Q.au0
Uu, Vu = Q.fe_space_u
opu = FEOperator(au0, ju0, Uu, Vu)
lb = -1e10*ones(Vu.nfree)
ub = interpolate_everywhere(Q.Φ₀ + Q.ϕ ⋅ Th, Uu).free_values
uh = SemismoothQVIs.hik(opu, uh, lb, ub, max_iter=1, damping=1, tol=1e-10, show_trace=true);

hik()

# aT = Q.aT[1]
# bT(Th, R) = aT(Th, R, SemismoothQVIs.projectionB(Q, xhs[1][1], (m*Rα, 0.0)))
# norm(Gridap.Algebra.residual(FEOperator(bT, UT, VT), xhs[1][2]))

# Semismooth Newton method
(zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, u₂, T₂; max_its=100, tol=1e-20, PF=false, globalization=false, proj_rc=(Inf,0.0), show_inner_trace=false, X=J);

# # Extract final solution
# uh = zhs1[its_1[1]][1]
# Th = zhs1[its_1[1]][2]
# mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
# # Check errors
# (h1(Q, u₁, uh), h1(Q, T₁, Th))
# (h1(Q, u₂, uh), h1(Q, T₂, Th))

uh = zhs3[its_3[1]+1][1]
Th = zhs3[its_3[1]+1][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
(h1(Q, u₁, uh), h1(Q, T₁, Th))
(h1(Q, u₂, uh), h1(Q, T₂, Th))

Th = T₂
uh = u₂
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ T₂, VT)
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