using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

n = 100
dΩ, Uu, Vu, UT, VT = FEM_1D_model(n)

k = π^2
Φ₀ = interpolate_everywhere(x->0.0, VT)
Ψ₀ = interpolate_everywhere(x->0.0, VT)
α = 1
ϕ = interpolate_everywhere(x->10*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
ψ = interpolate_everywhere(x->5*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
f = interpolate_everywhere(x->α*π^2*sin(π*x[1]), Uu)

g(s) = s ≥ 0.0 ? 0.0 : 4 / α * s^2
dg(s)= s ≥ 0.0 ? 0.0 : 8 / α * s

(u₁, T₁) = (interpolate_everywhere(x->0.0, Vu), interpolate_everywhere(x->0.0, VT))
(u₂, T₂) = (interpolate_everywhere(x->α*sin(π*x[1]), Vu), interpolate_everywhere(x->α*(5-cos(2π*x[1]))/( 10*π^2 ), VT))


Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)

uᵢ = FEFunction(Vu, ones(Vu.nfree))
Tᵢ = FEFunction(VT, ones(VT.nfree))
(zhs1, h1_1, its_1) = fixed_point(Q, uᵢ, Tᵢ; max_its=50, min_its=0, tol=1e-20, proj_rc=(Inf, 0), bt=false, PF=false, ρ0=1)
(zhs3, h1_3, its_3, is_3) = newtonss(Q, uᵢ, Tᵢ; max_its=100, tol=1e-20, globalization=true, proj_rc=(Inf,0.0), PF=false, bt=false);


uh = zhs1[its_1[1]][1]
Th = zhs1[its_1[1]][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

(h1(Q, u₁, uh), h1(Q, T₁, Th))
(h1(Q, u₂, uh), h1(Q, T₂, Th))

uh = zhs3[its_3[1]+1][1]
Th = zhs3[its_3[1]+1][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

(h1(Q, u₁, uh), h1(Q, T₁, Th))
(h1(Q, u₂, uh), h1(Q, T₂, Th))

xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mold(Point.(xx))],#  mold(Point.(xx)) Th(Point.(xx))],
    label=["Membrane" "Mold"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=2, xlabelfontsize=20)