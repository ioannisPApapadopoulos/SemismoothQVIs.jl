using Gridap, SemismoothQVIs
using Plots, LaTeXStrings

"""
Section 8 An alternative application: nonlinear VIs and PDEs

Given an obstacle Φ₀, find u ∈ H¹₀(Ω) that satisfies
    u ≤ Φ₀, ⟨-Δu - f - Φ(u), v - u⟩≥0 ∀ v∈H¹₀(Ω), v≤ Φ₀.

"""

# Nonlinear VI parameter choices
α = 1.0
b₁(s) = s < 0 ? 0.0 : α*s
db₁(s) = s < 0 ? 0.0 : α

b₂(s) = s + cos(s)
db₂(s) = 1.0 - sin(s)

Φ(u) =  b₁(u) * b₂(u)
dΦ(u) = b₁(u)*db₂(u) + db₁(u)*b₂(u)


h1s, isxBs, zhs = [], [], []
its = zeros(8, 4)
its2 = zeros(8, 4)
ns = [50, 100, 200, 400, 800, 1600, 3200, 6400]
for (n, i) in zip(ns, 1:8)
    print("Considering n=$n.\n")
    dΩ, U, V, _, _ = fem_model(n) # Piecewise linear FEM discretization on (0,1)

    f = interpolate_everywhere(x->50.0*sin(2π*x[1]), U)
    Φ₀ = interpolate_everywhere(x->1.0, U)
    u₀ = interpolate_everywhere(x->0.0, U)

    # Nonlinear VI struct
    VI = NonlinearVI(dΩ, Φ, dΦ, Φ₀, f, U)
    (zhs, h1_, its_, is_) = semismoothnewton(VI, u₀; max_its=50, in_tol=1e-10, out_tol=1e-10, PF=true, ρ_min=1e-8, globalization=true, show_inner_trace=false);
    (_, _, its2_, _) = semismoothnewton(VI, u₀; max_its=50, in_tol=1e-10, out_tol=1e-10, PF=true, ρ_min=1e-8, globalization=false, show_inner_trace=false);

    its[i,:] .= its_
    its2[i,:] .=its2_
    append!(h1s, [h1_])
    append!(isxBs, [is_[2]])
end

its_ns1 = hcat(round.( 1 ./ ns, digits=5), [its[:,1] its[:,3:end]])
its_ns2 = hcat(round.( 1 ./ ns, digits=5), [its2[:,1] its2[:,3:end]])
open("test4_table.log", "w") do file
    write(file, latex_table(its_ns1, its_ns2))
end

# Extract final solution
uh = zhs[end]
xx = range(0,1,100)
p = plot(xx, [uh(Point.(xx)) Φ₀(Point.(xx))],#  mold(Point.(xx)) Th(Point.(xx))],
    label=[L"u(x)" L"\Phi_0(x)"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=3,
    legend=:bottomleft,
    xlabelfontsize=20,legendfontsize=12,xtickfontsize=10,ytickfontsize=10,)
savefig("test4-solutions.pdf")