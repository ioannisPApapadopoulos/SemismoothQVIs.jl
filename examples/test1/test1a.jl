using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

n = 2000 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = fem_model(n) # Piecewise linear FEM discretization on (0,1)

# Thermoforming parameter choices
k = 1
Φ₀ = interpolate_everywhere(x->max(0, abs(x[1]-0.5)-0.25), VT)
Ψ₀ = interpolate_everywhere(x->1.0, VT)
α₁, α₂ = 1+1/π, 1.01
ϕ = interpolate_everywhere(x->sin(π*x[1])/α₁, UT)
ψ = ϕ
f = interpolate_everywhere(x->π^2*sin(π*x[1]) + max(0, -abs(x[1] - 0.625) + 0.125), Uu)
g(s) = (1-s) ≤ 0.0 ? α₁ + (1-s)/α₂ : α₁ + (1-s)/(2α₂)
dg(s)= (1-s) ≤ 0.0 ? -1/α₂ :  -1/(2α₂)

# Explicit solution
(u₁, T₁) = (interpolate_everywhere(x->sin(π*x[1]), Vu), interpolate_everywhere(x->α₁, VT))

# Thermoforming struct
Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)

# Initial guess
T₀ = interpolate_everywhere(x->5.0, UT)
J = Gridap.Algebra.jacobian(FEOperator(Q.au0..., Uu, Vu), interpolate_everywhere(x->0.0, Uu))

errs1, errs3, eocs1, eocs3 = [], [], [], []
for (u₀, j) in zip([interpolate_everywhere(x->0.0, Uu), FEFunction(Vu, J \ f.free_values)], [18,13])
    # fixed point iteration #13
    (zhs1, h1_1, its_1) = fixed_point(Q, u₀, T₀; max_its=j, tol=1e-13,  ρ0=1, PF=false, show_inner_trace=false);
    # Semismooth Newton method
    (zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, u₀, T₀; max_its=30, tol=1e-12, PF=false, globalization=false, show_inner_trace=false);
    err1, eoc1 = EOC(Q, first.(zhs1), u₁)
    err3, eoc3 = EOC(Q, first.(zhs3), u₁)
    append!(errs1, [err1]); append!(errs3, [err3]); append!(eocs1, [eoc1]); append!(eocs3, [eoc3])
end

ls = [:solid, :dash, :dashdot, :dashdotdot]
u0string = [L"u_0 = 0", L"u_0 = (-\Delta)^{-1} f"]
p = plot()
for i in 1:2
    p = Plots.plot!(1:lastindex(errs1[i]), errs1[i], linestyle=ls[i], linewidth=2, marker=:dot, label="Fixed Point  "*u0string[i], yscale=:log10,)
    Plots.plot!(1:lastindex(errs3[i]), errs3[i],
        linewidth=2,
        linestyle=ls[i],
        marker=:square,
        yscale=:log10,
        title=L"\alpha_1 = 1+\pi^{-1}, \alpha_2 = 101/100",
        label="Semismooth Newton  "*u0string[i],
        xlabel="Iterations" * L" $i$",
        ylabel=L"\Vert u_{i, h} - \bar u \, \Vert_{H^1_0(0,1)}",
        xlabelfontsize=15, ylabelfontsize=15, legendfontsize=8,xtickfontsize=10,ytickfontsize=10,
        yticks=10.0.^(-7:0), #[1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
        legend=:topright,
        ylim=[3e-8,1e0],
        xticks=1:18,
    )
end
display(p)
lim = 1; [annotate!(x+1.5, y+1e-7, Plots.text( "EOC=$(round(eoc3[i], digits=2))", 12)) for (i, x, y) in zip(1:lim, 3:lim+2, errs3[2][3:lim+2])]
plot!()
Plots.savefig("test1-convergence.pdf")

# Extract final solution
uh = zhs1[its_1[1]][1]
Th = zhs1[its_1[1]][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
# Check error
(h1(Q, u₁, uh), h1(Q, T₁, Th))

uh = zhs3[its_3[1]+1][1]
Th = zhs3[its_3[1]+1][2]
mold = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

(h1(Q, u₁, uh), h1(Q, T₁, Th))

# Plot solution
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mold(Point.(xx))],#  mold(Point.(xx)) Th(Point.(xx))],
    label=["Membrane" "Mold"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=3,
    legend=:bottom,
    xlabelfontsize=20,legendfontsize=12,xtickfontsize=10,ytickfontsize=10,)
Plots.savefig("test1-solutions.pdf")