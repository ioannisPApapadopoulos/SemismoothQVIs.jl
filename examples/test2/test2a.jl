using Gridap, SemismoothQVIs
using Plots, LaTeXStrings


# Thermoforming parameter choices
k = π^2
α₁ = 10.0
α₂ = 1.0
g(s) = s ≥ 0.0 ? 0.0 : 4 / α₁ * s^2
dg(s)= s ≥ 0.0 ? 0.0 : 8 / α₁ * s

RT = []
ns = 100:100:200
for n in ns
    print("Considering n=$n.\n")
    dΩ, Uu, Vu, UT, VT = fem_model(n) # Piecewise linear FEM discretization on (0,1)
 
    Φ₀ = interpolate_everywhere(x->0.0, VT)
    Ψ₀ = interpolate_everywhere(x->0.0, VT)
    ϕ = interpolate_everywhere(x->α₂*10*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
    ψ = interpolate_everywhere(x->5*π^2*sin(π*x[1])/(5-cos(2π*x[1])), UT)
    f = interpolate_everywhere(x->α₁*π^2*sin(π*x[1]), Uu)

    (u₀, T₀) = (interpolate_everywhere(x->α₁*sin(π*x[1]), Vu), interpolate_everywhere(x->α₁*(5-cos(2π*x[1]))/( 10*π^2 ), VT))

    # Thermoforming struct
    Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)
    aT = first(Q.aT)
    bT(Th, R) = aT(Th, R, u₀)
    append!(RT, norm(Gridap.Algebra.residual(FEOperator(bT, UT, VT), T₀)))
end

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