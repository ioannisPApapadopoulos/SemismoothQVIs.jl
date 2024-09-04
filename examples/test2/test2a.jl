using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

"""
Section 4.3 Test 2: a one-dimensional QVI with two known solutions

Find u ∈ H¹₀(Ω) that satisfies
    u ≤ Φ₀ + Φ(u), ⟨-Δu - f, v - u⟩≥0 ∀ v∈H¹₀(Ω), v≤ Φ₀ + Φ(u)
with Φ(u) given by Φ(u) := ϕT and T as the solution of
    kT - ΔT = g(Ψ₀ + ψT - u), ∂ν T = 0 on ∂Ω.
"""

n = 2000 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = fem_model(n) # Piecewise linear FEM discretization on (0,1)

# Thermoforming parameter choices
k = π^2
Φ₀ = interpolate_everywhere(x->0.0, VT)
Ψ₀ = interpolate_everywhere(x->0.0, VT)
α₁ = 10.0
α₂ = 1.0
ϕ = interpolate_everywhere(x->α₂*10*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
ψ = interpolate_everywhere(x->5*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
f = interpolate_everywhere(x->α₁*π^2*sin(π*x[1]), Uu)
g(s) = s ≥ 0.0 ? 0.0 : 4 / α₁ * s^2
dg(s)= s ≥ 0.0 ? 0.0 : 8 / α₁ * s

(u₁, T₁) = (interpolate_everywhere(x->0.0, Vu), interpolate_everywhere(x->0.0, VT))
(u₀, T₀) = (interpolate_everywhere(x->α₁*sin(π*x[1]), Vu), interpolate_everywhere(x->α₁*(5-cos(2π*x[1]))/( 10*π^2 ), VT))

# Thermoforming struct
Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)


(zhs1, h1_1, its_1) = fixed_point(Q, u₀, T₀; max_its=200, out_tol=1e-13, in_tol=1e-13, PF=true, bt=true, proj_rc=(Inf, 0.0), show_inner_trace=false);
uh, Th = zhs1[end]; (h1(Q, u₁, uh), h1(Q, T₁, Th))
err1a, _ = EOC(Q, first.(zhs1), u₁)
err1b, _ = EOC(Q, first.(zhs1), u₀)
(zhs2, h1_2, its_2, is_2) = semismoothnewton(Q, u₀, T₀; max_its=80, out_tol=1e-13, in_tol=1e-13, PF=true, globalization=true, proj_rc=(Inf,0.0), show_inner_trace=false);
err2a, _ = EOC(Q, first.(zhs2), u₁)
err2b, _ = EOC(Q, first.(zhs2), u₀)
(zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, u₀, T₀; max_its=80, out_tol=1e-13, in_tol=1e-13, PF=true, globalization=false, proj_rc=(Inf,0.0), show_inner_trace=false);
err3a, _ = EOC(Q, first.(zhs3), u₁)
err3b, _ = EOC(Q, first.(zhs3), u₀)
(zhs4, h1_4, its_4, is_4) = semismoothnewton(Q, u₀, T₀; max_its=80, out_tol=1e-13, in_tol=1e-13, PF=true, globalization=false, linesearch=true, proj_rc=(Inf,0.0), show_inner_trace=false);
err4a, _ = EOC(Q, first.(zhs4), u₁)
err4b, _ = EOC(Q, first.(zhs4), u₀)

ls = [:solid, :dash, :dashdot, :dashdotdot]
err1b[1] = NaN
err2b[1] = NaN
err3b[1] = NaN
err4b[1] = NaN
Plots.plot(0:lastindex(err1b)-1, err1b, linestyle=ls[1], linewidth=2, marker=:dot, label="Fixed Point", yscale=:log10,)
Plots.plot!(0:lastindex(err2b)-1, err2b, linestyle=ls[1], linewidth=2, marker=:dtriangle, label="Algorithm 2", yscale=:log10,)
p = Plots.plot!(0:lastindex(err3b)-1, err3b,
    linewidth=2,
    linestyle=ls[1],
    marker=:square,
    title=L"\alpha_1 = 10, \; \alpha_2 = 1, \; u_0 = I_{h} \bar u_2",
    label="Algorithm 1",
    xlabel="Iterations" * L" $i$",
    ylabel=L"\Vert u_{i, h} - \bar u_2 \, \Vert_{H^1_0(0,1)}",
    xlabelfontsize=15, ylabelfontsize=15, legendfontsize=8,xtickfontsize=10,ytickfontsize=10,
    legend=:right,
    xlim=[0,75],
    ylim=[1e-6, 5e1],
    yscale=:log10
)
Plots.plot!(0:lastindex(err4b)-1, err4b, linestyle=ls[1], linewidth=2, marker=:diamond, label="Algorithm 1 + backtracking")
Plots.savefig("test2-Fig1b.pdf")
