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
(zhs1_max, h1_1, its_1_max) = fixed_point(Q, u₀, T₀; max_its=50, in_tol=1e-13, out_tol=1e-13, PF=true, bt=true, proj_rc=(Inf, 0.0), show_inner_trace=false);
(zhs2_max, h1_2, its_2_max, is_2) = semismoothnewton(Q, u₀, T₀; max_its=100, in_tol=1e-13, out_tol=1e-13, PF=true, globalization=false, proj_rc=(Inf,0.0), show_inner_trace=false);

errs1, errs2, errs3, errs4 = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
eocs1, eocs2, eocs3, eocs4 = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
isxB = Vector{Bool}[]
zhs1_min, zhs2_min, zhs3_min, zhs4_min = [], [], [], []
Rs = range(0.01, h10(Q, u₂), 10)
for R in Rs
    # Solve via fixed point & semismooth Newton + projection, both converge to minimal solution
    print("Considering R=$R.\n")
    (zhs1_min, h1_1, its_1_min) = fixed_point(Q, u₀, T₀; max_its=20, in_tol=1e-15, out_tol=1e-15, PF=true, bt=true, proj_rc=(R, 0.0), show_inner_trace=false);
    (zhs2_min, h1_2, its_2_min, is_2) = semismoothnewton(Q, u₀, T₀; max_its=20, in_tol=1e-15, out_tol=1e-15, PF=true, globalization=true, proj_rc=(R,0.0), show_inner_trace=false);
    (zhs3_min, h1_3, its_3_min, is_3) = semismoothnewton(Q, u₀, T₀; max_its=20, in_tol=1e-15, out_tol=1e-15, PF=true, globalization=false, proj_rc=(R,0.0), show_inner_trace=false);
    (zhs4_min, h1_4, its_4_min, is_4) = semismoothnewton(Q, u₀, T₀; max_its=20, in_tol=1e-15, out_tol=1e-15, PF=true, globalization=false, proj_rc=(R,0.0), linesearch=true,show_inner_trace=false);

    err1, eoc1 = EOC(Q, first.(zhs1_min), u₁)
    err2, eoc2 = EOC(Q, first.(zhs2_min), u₁)
    err3, eoc3 = EOC(Q, first.(zhs3_min), u₁)
    err4, eoc4 = EOC(Q, first.(zhs4_min), u₁)
    push!(errs1, err1); push!(errs2, err2); push!(errs3, err3); push!(errs4, err4)
    push!(eocs1, eoc1); push!(eocs2, eoc2); push!(eocs3, eoc3); push!(eocs4, eoc4)
    push!(isxB, is_3[2])
end

its = [round.(Rs, digits=3) (length.(errs1).-1) (length.(errs2).-2) (length.(errs3).-2) (length.(errs4).-2)]'
eoc = round.([maximum.(eocs1) maximum.(eocs2) maximum.(eocs3) maximum.(eocs4)], digits=2)'
converge = [last.(errs1) .< 1e-12 last.(errs2) .< 1e-12  last.(errs3) .< 1e-12 last.(errs4) .< 1e-12 ]'
tab = test2_latex_table(its,eoc,converge)
open("test2_table.log", "w") do file
    write(file, tab)
end

uh, Th = zhs1_max[end]; (h1(Q, u₂, uh), h1(Q, T₂, Th))
uh, Th = zhs1_min[end]; (h1(Q, u₁, uh), h1(Q, T₁, Th))

uh, Th = zhs2_max[end]; (h1(Q, u₂, uh), h1(Q, T₂, Th))
uh, Th = zhs2_min[end]; (h1(Q, u₁, uh), h1(Q, T₁, Th))

mould = interpolate_everywhere(Φ₀ + ϕ ⋅ Th, VT)

# Plot solution
xx = range(0,1,50)
p = plot(xx, [uh(Point.(xx)) mould(Point.(xx))],
    label=["Membrane" "Mould"],
    linestyle=[:solid :dash],
    xlabel=L"x",
    linewidth=2, xlabelfontsize=20)
