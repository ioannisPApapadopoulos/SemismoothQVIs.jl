using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

n = 2000 # dofs, h = 1/n
dΩ, Uu, Vu, UT, VT = fem_model(n) # Piecewise linear FEM discretization on (0,1)

# Thermoforming parameter choices
k = π^2
Φ₀ = interpolate_everywhere(x->0.0, VT)
Ψ₀ = interpolate_everywhere(x->0.0, VT)
α₁ = 10.0
α₂ = 1.5
ϕ = interpolate_everywhere(x->α₂*10*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
ψ = interpolate_everywhere(x->5*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
f = interpolate_everywhere(x->α₁*π^2*sin(π*x[1]), Uu)
g(s) = s ≥ 0.0 ? 0.0 : 4 / α₁ * s^2
dg(s)= s ≥ 0.0 ? 0.0 : 8 / α₁ * s


# Radius of (u,T) = (0,0) solution
Rα = 3α₁*  ( sqrt( ((5α₂ + 8)*π^2 + 8π) / (80α₂) )  - π/4   ) /  (10*(1+π)) 
50α₂/3/α₁ *(Rα + 20*(1+π)/3/π/α₁ * Rα^2) ≈ 1


# Explicit solutions
(u₁, T₁) = (interpolate_everywhere(x->0.0, Vu), interpolate_everywhere(x->0.0, VT))
(u₂, T₂) = (interpolate_everywhere(x->α₁*sin(π*x[1]), Vu), interpolate_everywhere(x->α₁*(5-cos(2π*x[1]))/( 10*π^2 ), VT))

# Thermoforming struct
Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)

# Initial guess
T₀ = interpolate_everywhere(x->5.0, UT)
u₀ = interpolate_everywhere(x->1e2*(1.0 -x[1])*x[1], Uu)

u₂_norm = h10(Q, u₂)

# Solve via fixed point & semismooth Newton, both converge to maximal solution
(zhs1_max, h1_1, its_1_max) = fixed_point(Q, u₀, T₀; max_its=50, tol=1e-13, PF=true, bt=true, proj_rc=(Inf, 0.0), show_inner_trace=false);
(zhs2_max, h1_2, its_2_max, is_2) = semismoothnewton(Q, u₀, T₀; max_its=100, tol=1e-13, PF=true, globalization=true, proj_rc=(Inf,0.0), show_inner_trace=false);

errs1, errs2, errs3, eocs1, eocs2, eocs3, isxB, zhs1_min, zhs2_min, zhs3_min = [], [], [], [], [], [], [], [], [], []
Rs = range(0.01, 12.3, 10)
for R in Rs
# R = 0.1
    # Solve via fixed point & semismooth Newton + projection, both converge to minimal solution
    print("Considering R=$R.\n")
    (zhs1_min, h1_1, its_1_min) = fixed_point(Q, u₀, T₀; max_its=20, tol=1e-15, PF=true, bt=true, proj_rc=(R, 0.0), show_inner_trace=false);
    (zhs2_min, h1_2, its_2_min, is_2) = semismoothnewton(Q, u₀, T₀; max_its=20, tol=1e-15, PF=true, globalization=true, proj_rc=(R,0.0), show_inner_trace=false);
    (zhs3_min, h1_3, its_3_min, is_3) = semismoothnewton(Q, u₀, T₀; max_its=20, tol=1e-15, PF=true, globalization=false, proj_rc=(R,0.0), show_inner_trace=false);

    err1, eoc1 = EOC(Q, first.(zhs1_min), u₁)
    err2, eoc2 = EOC(Q, first.(zhs2_min), u₁)
    err3, eoc3 = EOC(Q, first.(zhs3_min), u₁)
    append!(errs1, [err1]); 
    append!(errs2, [err2]);
    append!(errs3, [err3]);
    append!(eocs1, [eoc1]); 
    append!(eocs2, [eoc2]);
    append!(eocs3, [eoc3]);
    append!(isxB, [is_3[2]])
end

its = [round.(Rs, digits=3) (length.(errs1).-1) (length.(errs2).-2) (length.(errs3).-2)]
brackets = round.([maximum.(eocs1) maximum.(eocs2) maximum.(eocs3)], digits=2)
tab = latex_table(its,cap=20, brackets=brackets)
open("test2_table.log", "w") do file
    write(file, tab)
end
# ls = [:solid, :dash, :dashdot, :dashdotdot]
Rs = 0.1:0.1:0.6
u0string = [L"$r=%$R$" for R in Rs]
p = plot()
for i in 1:6
    Plots.plot!(0:lastindex(errs1[i])-1, errs1[i],  linewidth=2, marker=:dot, label="Fixed Point  "*u0string[i], yscale=:log10,)
    p = Plots.plot!(0:lastindex(errs2[i])-1, errs2[i],
        linewidth=2,
        # linestyle=ls[i],
        marker=:square,
        title=L"\alpha_1 = 10, \alpha_2 = 3/2",
        label="Semismooth Newton  "*u0string[i],
        xlabel="Iterations" * L" $i$",
        ylabel=L"\Vert u_{i, h} - \bar u_1 \, \Vert_{H^1_0(0,1)}",
        xlabelfontsize=15, ylabelfontsize=15, legendfontsize=8,xtickfontsize=10,ytickfontsize=10,
        # yticks=10.0.^(-7:2:4), #[1e-7,1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
        legend=:topright,
        # ylim=[1e-8,1e4],
        # xlim=[0,19],
        # xticks=0:2:18,
        yscale=:log10,
    )
end
    
display(p)



uh, Th = zhs1_max[end]; (h1(Q, u₂, uh), h1(Q, T₂, Th))
uh, Th = zhs1_min[end]; (h1(Q, u₁, uh), h1(Q, T₁, Th))

uh, Th = zhs2_max[end]; (h1(Q, u₂, uh), h1(Q, T₂, Th))
uh, Th = zhs2_min[end]; (h1(Q, u₁, uh), h1(Q, T₁, Th))

mould = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

# Plot solution
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mould(Point.(xx))],#  mould(Point.(xx)) Th(Point.(xx))],
    label=["Membrane" "Mould"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=2, xlabelfontsize=20)


# ls = [:solid, :dash, :dashdot, :dashdotdot]
# p = Plots.plot(1:length(h1_1), h1_1, linestyle=ls[1], linewidth=2, marker=:dot, label="Fixed Point", yscale=:log10,)
# Plots.plot!(1:length(h1_2), h1_2,
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