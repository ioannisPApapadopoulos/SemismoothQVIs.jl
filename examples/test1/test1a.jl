using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

"""
Section 4.2 Test 1: a one-dimensional QVI with a known solution

Find u ∈ H¹₀(Ω) that satisfies
    u ≤ Φ₀ + Φ(u), ⟨-Δu - f, v - u⟩≥0 ∀ v∈H¹₀(Ω), v≤ Φ₀ + Φ(u)
with Φ(u) given by Φ(u) := ϕT and T as the solution of
    kT - ΔT = g(Ψ₀ + ψT - u), ∂ν T = 0 on ∂Ω.
"""

n = 2000 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = fem_model(n) # Piecewise linear FEM discretization on (0,1)

# Thermoforming parameter choices
k = 1
Φ₀ = interpolate_everywhere(x->max(0, abs(x[1]-0.5)-0.25), VT)
Ψ₀ = interpolate_everywhere(x->1.0, VT)
α₁, α₂ = 1+1/π, 1.01
ϕ = interpolate_everywhere(x->sin(π*x[1])/α₁, UT)
ψ = ϕ
f = interpolate_everywhere(x->π^2*sin(π*x[1]) + 100*max(0, -abs(x[1] - 0.625) + 0.125), Uu)
g(s) = s ≤ 1.0 ? α₁ + atan( (1-s)/ (2α₂)) : α₁ + atan( (1-s)/ α₂)
dg(s)= s ≤ 1.0 ? -2α₂/(4*α₂^2 + (1-s)^2) :  -α₂/(α₂^2 + (1-s)^2)

# Explicit solution
(u₁, T₁) = (interpolate_everywhere(x->sin(π*x[1]), Vu), interpolate_everywhere(x->α₁, VT))

# Thermoforming struct
Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)

# Initial guess
T₀ = interpolate_everywhere(x->5.0, UT)

# Hack for obtaining stiffness matrix J ≈ -Δ
J = Gridap.Algebra.jacobian(FEOperator(Q.au0..., Uu, Vu), interpolate_everywhere(x->0.0, Uu))

errs1, errs2, errs3 =  Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
eocs1, eocs2, eocs3 = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
isxB_2, isxB_3 = Vector{Bool}[], Vector{Bool}[]
zhs1, zhs2, zhs3 = [], [], []
for (u₀, j) in zip([interpolate_everywhere(x->0.0, Uu), FEFunction(Vu, J \ f.free_values)], [18,13])
    # fixed point iteration
    (zhs1, h1_1, its_1) = fixed_point(Q, u₀, T₀; max_its=j, out_tol=1e-13, in_tol=1e-13, ρ0=1, PF=true, show_inner_trace=false);
    # Semismooth Newton method
    (zhs2, h1_2, its_2, is_2) = semismoothnewton(Q, u₀, T₀; max_its=30, out_tol=1e-12, in_tol=1e-12, PF=true, globalization=true, show_inner_trace=false);
    (zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, u₀, T₀; max_its=30, out_tol=1e-12, in_tol=1e-12, PF=true, globalization=false, linesearch=true,show_inner_trace=false);
    err1, eoc1 = EOC(Q, first.(zhs1), u₁)
    err2, eoc2 = EOC(Q, first.(zhs2), u₁)
    err3, eoc3 = EOC(Q, first.(zhs3), u₁)
    push!(errs1, err1); push!(errs2, err2); push!(errs3, err3); 
    push!(eocs1, eoc1); push!(eocs2, eoc2); push!(eocs3, eoc3); 
    push!(isxB_2, is_2[2]); push!(isxB_3, is_3[2])
end

ls = [:solid, :dash, :dashdot, :dashdotdot]
u0string = [L"u_0 = 0", L"u_0 = (-\Delta)^{-1} f"]
p = plot()
display(p)
for i in 1:2
    Plots.plot!(0:lastindex(errs1[i])-1, errs1[i], linestyle=ls[i], linewidth=2, marker=:dot, label="Fixed Point  "*u0string[i], yscale=:log10,)
    p = Plots.plot!(0:lastindex(errs2[i])-1, errs2[i],
        linewidth=2,
        linestyle=ls[i],
        marker=:square,
        title=L"\alpha_1 = 1+\pi^{-1}, \alpha_2 = 101/100",
        label="Algorithm 2 "*u0string[i],
        xlabel="Iterations" * L" $i$",
        ylabel=L"\Vert u_{i, h} - \bar u \, \Vert_{H^1_0(0,1)}",
        xlabelfontsize=15, ylabelfontsize=15, legendfontsize=8,xtickfontsize=10,ytickfontsize=10,
        yticks=10.0.^(-7:2:4), #[1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
        legend=:topright,
        ylim=[1e-8,1e4],
        xlim=[0,19],
        xticks=0:2:18,
        yscale=:log10
    )
end
display(p)
lim = 2; [annotate!(x+1, y+yp , Plots.text( "$(round.(eocs2[2][i+1], digits=2))", 12, theme_palette(:auto)[4])) for (i, x, y, yp) in zip(0:lim-1, 2:lim+1, errs2[2][3:lim+2], [-1e-2, -1e-5])]
lim = 2; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:lim-1, 2:lim+1, errs2[2][3:lim+2])]
lim = 2; [annotate!(x-1.1, y+yp , Plots.text( "$(round.(eocs2[1][i+1], digits=2))", 12, theme_palette(:auto)[2])) for (i, x, y, yp) in zip(0:lim-1, 2:lim+1, errs2[1][3:lim+2], [-1e-5, -1e-8])]
lim = 2; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:lim-1, 2:lim+1, errs2[1][3:lim+2])]
plot!()
Plots.savefig("test1-atan-Fig1a.pdf")

# Extract final solution
uh = zhs1[end][1]
Th = zhs1[end][2]
mould = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
# Check error
(h1(Q, u₁, uh), h1(Q, T₁, Th))

uh = zhs2[end][1]
Th = zhs2[end][2]
mould = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

(h1(Q, u₁, uh), h1(Q, T₁, Th))

# Plot solution
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mould(Point.(xx))],#  mould(Point.(xx)) Th(Point.(xx))],
    label=["Membrane  "*L"\bar u" "Mould  "*L"\Phi_0 + \varphi \bar T"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=3,
    legend=:bottom,
    xlabelfontsize=20,legendfontsize=12,xtickfontsize=10,ytickfontsize=10,)
Plots.savefig("test1-solutions.pdf")