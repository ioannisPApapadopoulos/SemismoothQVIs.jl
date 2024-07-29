using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

"""
Section 7.2 Test 1: a one-dimensional QVI with a known solution

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
α₁, α₂ = 1e-2, 1e-2
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
J = Gridap.Algebra.jacobian(FEOperator(Q.au0..., Uu, Vu), interpolate_everywhere(x->0.0, Uu))

errs1, errs2, errs3, errs4 =  Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
eocs1, eocs2, eocs3, eocs4 = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
isxB_2, isxB_3, isxB_4 = Vector{Bool}[], Vector{Bool}[], Vector{Bool}[]
zhs1, zhs2, zhs3, zhs4 = [], [], [], []

for (u₀, j) in zip([interpolate_everywhere(x->0.0, Uu), FEFunction(Vu, J \ f.free_values)], [18,18])
    # fixed point iteration
    (zhs1, h1_1, its_1) = fixed_point(Q, u₀, T₀; max_its=j, out_tol=1e-13, in_tol=1e-13,  ρ0=1, PF=true, show_inner_trace=false);
    # Algorithm 1
    (zhs2, h1_2, its_2, is_2) = semismoothnewton(Q, u₀, T₀; max_its=j-1, out_tol=1e-12, in_tol=1e-12, PF=true, globalization=true, show_inner_trace=false);
    # vanilla SSN
    (zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, u₀, T₀; max_its=j-1, out_tol=1e-12, in_tol=1e-12, PF=true, globalization=false, show_inner_trace=false);
    # Armijo backtracking linesearch
    (zhs4, h1_4, its_4, is_4) = semismoothnewton(Q, u₀, T₀; max_its=j-1, out_tol=1e-12, in_tol=1e-12, PF=true, globalization=false, show_inner_trace=false, linesearch=true);
    
    err1, eoc1 = EOC(Q, first.(zhs1), u₁)
    err2, eoc2 = EOC(Q, first.(zhs2), u₁)
    err3, eoc3 = EOC(Q, first.(zhs3), u₁)
    err4, eoc4 = EOC(Q, first.(zhs4), u₁)

    push!(errs1, err1); push!(errs2, err2); push!(errs3, err3); push!(errs4, err4); 
    push!(eocs1, eoc1); push!(eocs2, eoc2); push!(eocs3, eoc3); push!(eocs4, eoc4); 
    push!(isxB_2, is_2[2]); push!(isxB_3, is_3[2]); push!(isxB_4, is_4[2]);
end

ls = [:solid, :dash, :dashdot, :dashdotdot]
u0string = [L"u_0 = 0", L"u_0 = (-\Delta)^{-1} f"]
p = plot()
for i in 1:2
    Plots.plot!(0:lastindex(errs1[i])-1, errs1[i], linestyle=ls[i], linewidth=2, marker=:dot, label="Fixed Point  "*u0string[i], yscale=:log10,)
    p = Plots.plot!(0:lastindex(errs2[i])-1, errs2[i],
        linewidth=2,
        linestyle=ls[i],
        marker=:square,
        markeralpha= i == 1 ? 1.0 : 0.5,
        title=L"\alpha_1 = 10^{-2}, \alpha_2 = 10^{-2}",
        label="Algorithm 1  "*u0string[i],
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
lim = 18; [annotate!(x, y-0.55, Plots.text( "$(round(eocs2[2][i+1], digits=2))", 8, theme_palette(:auto)[4])) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs2[2][3:2:lim+2])]
lim = 18; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs2[2][3:2:lim+2])]
lim = 18; [annotate!(x, y-3e-5, Plots.text( "$(round(eocs2[1][i+1], digits=2))", 8, theme_palette(:auto)[2])) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs2[1][3:2:lim+2])]
lim = 18; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs2[1][3:2:lim+2])]
plot!()
Plots.savefig("test1-atan-Fig1b.pdf")

p = plot()
for i in 1:2
     p = Plots.plot!(0:lastindex(errs3[i])-1, errs3[i],
        linewidth=2,
        linestyle=ls[i],
        marker=:square,
        title=L"\alpha_1 = 10^{-2}, \alpha_2 = 10^{-2}",
        label="Vanilla SSN  "*u0string[i],
        xlabel="Iterations" * L" $i$",
        ylabel=L"\Vert u_{i, h} - \bar u \, \Vert_{H^1_0(0,1)}",
        xlabelfontsize=15, ylabelfontsize=15, legendfontsize=8,xtickfontsize=10,ytickfontsize=10,
        yticks=10.0.^(-7:2:4), #[1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
        legend=:topright,
        ylim=[1e-8,1e4],
        xlim=[0,19],
        xticks=0:2:18,
        yscale=:log10,
        color=theme_palette(:auto)[2*i]
    )
    p = Plots.plot!(0:lastindex(errs4[i])-1, errs4[i], linestyle=ls[i], 
        linewidth=2, marker=:dtriangle, 
        label="Backtracking SSN  "*u0string[i],
        color=theme_palette(:auto)[2*i+1])
end
# plot(1:5, rand(5); yscale=:log10)
display(p)
lim = 14; [annotate!(x, y-yp, Plots.text( "$(round(eocs3[2][i+1], digits=2))", 8, theme_palette(:auto)[4])) for (i, x, y, yp) in zip(0:2:lim-1, 2:2:lim+1, errs3[2][3:2:lim+2], errs3[2][3:2:lim+2]/1.2 )]
lim = 14; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs3[2][3:2:lim+2])]

lim = 18; [annotate!(x, y+5*yp, Plots.text( "$(round(eocs4[2][i+1], digits=2))", 8, theme_palette(:auto)[5])) for (i, x, y, yp) in zip(0:2:lim-1, 2:2:lim+1, errs4[2][3:2:lim+2], errs4[2][3:2:lim+2]/1.2 )]
lim = 18; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs4[2][3:2:lim+2])]

lim = 5; [annotate!(x, y-yp, Plots.text( "$(round(eocs3[1][i+1], digits=2))", 8, theme_palette(:auto)[2])) for (i, x, y, yp) in zip(0:2:lim-1, 2:2:lim+1, errs3[1][3:2:lim+2], errs3[1][3:2:lim+2]/1.2)]
lim = 5; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs3[1][3:2:lim+2])]

lim = 16; [annotate!(x, y+3*yp, Plots.text( "$(round(eocs4[1][i+1], digits=2))", 8, theme_palette(:auto)[3])) for (i, x, y, yp) in zip(0:2:lim-1, 2:2:lim+1, errs4[1][3:2:lim+2], errs4[1][3:2:lim+2]/1.2 )]
lim = 16; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:2:lim-1, 2:2:lim+1, errs4[1][3:2:lim+2])]

plot!()
Plots.savefig("test1-atan-Fig1c.pdf")

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

uh = zhs3[end][1]
Th = zhs3[end][2]
mould = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
(h1(Q, u₁, uh), h1(Q, T₁, Th))

uh = zhs4[end][1]
Th = zhs4[end][2]
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