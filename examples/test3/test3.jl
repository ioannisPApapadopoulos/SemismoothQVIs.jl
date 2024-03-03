using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

n = 50
dΩ, Uu, Vu, UT, VT = fem_model(n, n)

include("test3_parameters_pyramid.jl")

Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)
uᵢ = FEFunction(Vu, zeros(Vu.nfree))
Tᵢ = FEFunction(VT, ones(VT.nfree))

(zhs1, h1_1, its_1) = fixed_point(Q, uᵢ, Tᵢ; max_its=14, min_its=0, tol=1e-15, bt=true, PF=true, ρ0=1, show_inner_trace=false)
(zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, uᵢ, Tᵢ; max_its=6, tol=1e-15, globalization=false, PF=true, bt=false, show_inner_trace=false);


h1_2s = []
for ρ in [1e-2, 1e-4, 1e-6, 1e-8]
    (zhs2, h1_2, its_2, is_2) = moreau_yosida_newton(Q, uᵢ, Tᵢ;max_its=9, inner_max_its=60, ρ=ρ, tol=1e-20, globalization=true, bt=true, show_inner_trace=false);
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
    # title="Slice at  "*L"x_2=1/2",
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
p = surface(xx, xx, (x, y) -> mould(Point.(x,y)), 
    color=:redsblues, #:vik,
    xlabel=L"x", ylabel=L"y", zlabel=L"(\Phi_0 + \varphi T)(x,y)",
    # camera=(30,-30),
    title="Mold  "*L"\Phi_0 + \varphi T",
    margin=(-6, :mm),
    zlim=[0,1.25],
)
Plots.savefig("pyr-thermoforming-mould.pdf")




# using DelimitedFiles
# h1_1s = [readdlm("h1_1s.log")[:]]
# h1_3s = [readdlm("h1_3s.log")[:]]

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
    ylim=[1e-16,5e0])
for (i, ρ)  in zip(1:4, [1e-2, 1e-4, 1e-6, 1e-8])
    p = Plots.plot!(1:length(h1_2s[i]), h1_2s[i], linestyle=ls[i], linewidth=2, marker=:dtriangle, label="MY-Newton (ρ = $ρ)", yscale=:log10,)
    display(p)
end

Plots.savefig("tf-convergence.pdf")


h1_3s = []
its = zeros(7, 4)
ns = [25,50,100,150,200,250,300]
for (n, i) in zip(ns, 1:7)
    dΩ, Uu, Vu, UT, VT = FEM_2D_model(n)
    include("test3_parameters_pyramid.jl")

    uᵢ = FEFunction(Vu, zeros(Vu.nfree))
    Tᵢ = FEFunction(VT, ones(VT.nfree))

    Q = GeneralizedThermoformingQVI(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, Uu, UT)
    (zhs3, h1_3, its_3, is_3) = semismoothnewton(Q, uᵢ, Tᵢ; max_its=4, tol=1e-20, globalization=false, show_inner_trace=false);

    # append!(its, its_3)
    its[i,:] .= its_3
    append!(h1_3s, [h1_3])
end

its_ns = hcat(round.( 1 ./ ns, digits=4), its)
text = SemismoothQVIs.latex_table(its_ns)
open("test3_table.log", "w") do file
    write(file, text)
end