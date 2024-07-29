using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

"""
Section 7.4 Test 3: the thermoforming QVI in two dimensions

Find u ∈ H¹₀(Ω) that satisfies
    u ≤ Φ₀ + Φ(u), ⟨-Δu - f, v - u⟩≥0 ∀ v∈H¹₀(Ω), v≤ Φ₀ + Φ(u)
with Φ(u) given by Φ(u) := ϕT and T as the solution of
    kT - ΔT = g(Ψ₀ + ψT - u), ∂ν T = 0 on ∂Ω.
"""

n = 50
dΩ, Uu, Vu, UT, VT = fem_model(n, n)

include("test3_parameters_pyramid.jl")

Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)
uᵢ = FEFunction(Vu, zeros(Vu.nfree))
Tᵢ = FEFunction(VT, ones(VT.nfree))

(zhs1, h1_1, its_1) = fixed_point(Q, uᵢ, Tᵢ; max_its=14, min_its=0, in_tol=1e-20, out_tol=1e-20, bt=true, PF=true, ρ0=1, show_inner_trace=false)
(zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, uᵢ, Tᵢ; max_its=5, in_tol=1e-20, out_tol=1e-20, globalization=true, PF=true, bt=true, show_inner_trace=false);


h1_2s = []
for ρ in [1e-2, 1e-4, 1e-6, 1e-8]
    (zhs2, h1_2, its_2, is_2) = moreau_yosida_newton(Q, uᵢ, Tᵢ;max_its=9, inner_max_its=60, ρ=ρ, in_tol=1e-20, out_tol=1e-20, globalization=false, bt=true, show_inner_trace=false);
    append!(h1_2s, [h1_2])
end


uh = zhs3[end][1]
Th = zhs3[end][2]
mould = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)
deform = interpolate_everywhere(ϕ ⋅ Th, VT)

# writevtk(Ω,"u", cellfields=["membrane"=>zhs3[length(h1_3)][1]])

xx = range(0,1,100)
p = plot(xx, [uh(Point.(xx,0.5)) mould(Point.(xx,0.5)) Th(Point.(xx,0.5)) deform(Point.(xx,0.5))],
    linewidth=2,
    label=["Membrane" "Mould" "Temperature" "Mould Deformation"],
    linestyle=[:solid :dash],
    xlabel=L"x_1",
    xlabelfontsize=20,
    ylim=[0,1.25])
Plots.savefig("pyr-thermoforming-slice.pdf")

xx = range(0,1,100)
Plots.gr_cbar_offsets[] = (-0.05,-0.01)
Plots.gr_cbar_width[] = 0.03
p = surface(xx, xx, (x, y) -> uh(Point.(x,y)), 
    color=:redsblues, #:vik,
    xlabel=L"x_1", ylabel=L"x_2", zlabel=L"u(x_1,x_2)",
    # camera=(30,-30),
    title="Membrane  "*L"u",
    margin=(-6, :mm),
    zlim=[0,1.25],
)
Plots.savefig("pyr-thermoforming-membrane.pdf")

xx = range(0,1,150)
Plots.gr_cbar_offsets[] = (-0.05,-0.01)
Plots.gr_cbar_width[] = 0.03
p = surface(xx, xx, (x, y) -> mould(Point.(x,y)), 
    color=:redsblues, #:vik,
    xlabel=L"x", ylabel=L"y", zlabel=L"(\Phi_0 + \varphi T)(x,y)",
    # camera=(30,-30),
    title="Mould  "*L"\Phi_0 + \varphi T",
    margin=(-6, :mm),
    zlim=[0,1.25],
)
Plots.savefig("pyr-thermoforming-mould.pdf")


ls = [:solid, :dash, :dashdot, :dashdotdot]
p = Plots.plot(1:length(h1_1), h1_1, linestyle=ls[1], linewidth=2, marker=:dot, label="Fixed Point", yscale=:log10,)
Plots.plot!(1:length(h1_3), h1_3,
    linewidth=2,
    linestyle=ls[1],
    marker=:square,
    yscale=:log10,
    label="Algorithm 1",
    xlabel="Iterations",
    ylabel=L"\Vert R(u_{i}) \Vert_{H^1}",
    yticks=[1e-15,1e-10,1e-5,1e0],
    xtick=2:2:22,
    ylim=[1e-17,5e0])
for (i, ρ)  in zip(1:4, [1e-2, 1e-4, 1e-6, 1e-8])
    p = Plots.plot!(1:length(h1_2s[i]), h1_2s[i], linestyle=ls[i], linewidth=2, marker=:dtriangle, label="MY-Newton (ρ = $ρ)", yscale=:log10,)
    display(p)
end

Plots.savefig("pyr-tf-convergence.pdf")


h1s = []
its = zeros(7, 4)
its2 = zeros(7, 4)
ns = [25,50,100,150,200,250,300]

for (n, i) in zip(ns, 1:7)
    print("Considering n=$n.\n")
    dΩ, Uu, Vu, UT, VT = fem_model(n,n)
    include("test3_parameters_pyramid.jl")

    uᵢ = FEFunction(Vu, zeros(Vu.nfree))
    Tᵢ = FEFunction(VT, ones(VT.nfree))

    Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)
    (zhs, h1_, its_, is_) = semismoothnewton(Q, uᵢ, Tᵢ; max_its=4, in_tol=1e-15, out_tol=1e-12, globalization=true, show_inner_trace=false);
    (_, _, its2_, _) = semismoothnewton(Q, uᵢ, Tᵢ; max_its=4, in_tol=1e-15, out_tol=1e-12, globalization=false, show_inner_trace=true);


    its[i,:] .= its_
    its2[i,:] .= its2_
    append!(h1s, [h1_])
end

its_ns1 = hcat(round.( 1 ./ ns, digits=5), its)
its_ns2 = hcat(round.( 1 ./ ns, digits=5), its2)
open("test3_table.log", "w") do file
    write(file, latex_table(its_ns1, its_ns2))
end