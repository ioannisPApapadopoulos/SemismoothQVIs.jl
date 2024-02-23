using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

n = 100 # dofs, h = 1/n
dΩ, U, V, _, _ = fem_model(n) # Piecewise linear FEM discretization on (0,1)

# Nonlinear VI parameter choices
α = 2.0
Φ(u) = u < 0 ? 0.0 : α*u*cos(u)
dΦ(u) = u < 0 ? 0.0 : α*cos(u) - α*u*sin(u)
φ = interpolate_everywhere(x->0.2, U)
f = interpolate_everywhere(x->10.0*sin(2π*x[1]), U)

# Nonlinear VI struct
VI = NonlinearVI(dΩ, Φ, dΦ, φ, f, U)

# Initial guess

u₀ = interpolate_everywhere(x->0.0, U)
(zhs1, h1_1, its_1) = fixed_point(VI, u₀; max_its=20, tol=1e-11,  ρ0=1, PF=true, show_inner_trace=false);
(zhs2, h1_2, its_2) = visolver(VI, u₀; max_its=20, tol=1e-11,  ρ0=1, PF=true, show_inner_trace=true);
(zhs3, h1_3, its_3, is_3) = semismoothnewton(VI, u₀; max_its=30, tol=1e-12, PF=true, globalization=false, show_inner_trace=false);

# Extract final solution
uh = zhs3[8]
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx))],#  mold(Point.(xx)) Th(Point.(xx))],
    label=L"u",
    linestyle=[:solid],
    xlabel=L"x",
    linewidth=3,
    legend=:bottom,
    xlabelfontsize=20,legendfontsize=12,xtickfontsize=10,ytickfontsize=10,)


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
    p = Plots.plot!(0:lastindex(errs1[i])-1, errs1[i], linestyle=ls[i], linewidth=2, marker=:dot, label="Fixed Point  "*u0string[i], yscale=:log10,)
    Plots.plot!(0:lastindex(errs3[i])-1, errs3[i],
        linewidth=2,
        linestyle=ls[i],
        marker=:square,
        yscale=:log10,
        title=L"\alpha_1 = 1+\pi^{-1}, \alpha_2 = 101/100",
        label="Semismooth Newton  "*u0string[i],
        xlabel="Iterations" * L" $i$",
        ylabel=L"\Vert u_{i, h} - \bar u \, \Vert_{H^1_0(0,1)}",
        xlabelfontsize=15, ylabelfontsize=15, legendfontsize=8,xtickfontsize=10,ytickfontsize=10,
        yticks=10.0.^(-7:2:4), #[1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
        legend=:topright,
        ylim=[3e-8,1e4],
        xticks=0:2:18,
        xlim=[0,18.5]
    )
end
display(p)
lim = 2; [annotate!(x+1, y+3e-7, Plots.text( "$(round.(eocs3[2][i+1], digits=2))", 12)) for (i, x, y) in zip(0:lim-1, 2:lim+1, errs3[2][3:lim+2])]
lim = 2; [annotate!(x, y, Plots.text( "O", 17)) for (i, x, y) in zip(0:lim-1, 2:lim+1, errs3[2][3:lim+2])]
plot!()
Plots.savefig("test1a-convergence.pdf")

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