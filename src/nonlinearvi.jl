"""
    NonlinearVI()

Given an obstacle Φ₀, find u ∈ H¹₀(Ω) that satisfies
    u ≤ Φ₀, ⟨-Δu - f - Φ(u), v - u⟩≥0 ∀ v∈H¹₀(Ω), v≤ Φ₀.

"""
# 1e-10 = 1e-10
struct NonlinearVI{T}
    dΩ::Gridap.CellData.GenericMeasure
    Φ::Function
    dΦ::Function
    Φ₀::Gridap.FESpaces.SingleFieldFEFunction
    f::Gridap.FESpaces.SingleFieldFEFunction
    au::Tuple{Function, Function, Function}
    au0::Tuple{Function, Function, Function}
    asn::Tuple{Function, Function}
    asnp::Tuple{Function, Function}
    amy::Tuple{Function, Function}
    fe_space_u::Tuple{FESpace, FESpace}
end

Xinner(VI::NonlinearVI{T}, u, v) where T = sum(∫(∇(u) ⋅ ∇(v))*VI.dΩ)
Xnorm(VI::NonlinearVI{T}, u) where T = sqrt(Xinner(VI, u, u))
# _inner_X(dΩ, u, h) = sum(∫(∇(u) ⋅ ∇(h))*dΩ)
# _dpB = (dΩ, h, u, ur) -> ur[2]/ur[1] * ( h - ((h->_inner_X(dΩ, u, h)) ∘ h) ⋅ u /ur[1]^2)

function NonlinearVI(dΩ::Gridap.CellData.GenericMeasure, Φ::Function, dΦ::Function,
    Φ₀::Gridap.FESpaces.SingleFieldFEFunction, f::Gridap.FESpaces.SingleFieldFEFunction, Uu::FESpace)


    
    # au0(uh, v, u_) =∫(∇(uh) ⋅ ∇(v) - (Φ ∘ u_) ⋅ v - f ⋅ v) * dΩ
    # 

    au(uh, v, s, u_) =∫( ∇(uh) ⋅ ∇(v) +  (s ∘ (uh - Φ₀)) ⋅ v - f ⋅ v - (Φ ∘ u_) ⋅ v) * dΩ
    ju(uh, duh, v, ds) =∫( ∇(duh) ⋅ ∇(v) +  (ds ∘ (uh - Φ₀)) ⋅ duh ⋅ v) * dΩ
    juu(uh, duh, v, ds, u_) =∫( ∇(duh) ⋅ ∇(v) +  (ds ∘ (uh - Φ₀)) ⋅ duh ⋅ v - (dΦ ∘ u_) ⋅ duh ⋅ v) * dΩ
    
    au0(uh, v, u_) =∫(∇(uh) ⋅ ∇(v) - (Φ ∘ u_) ⋅ v - f ⋅ v) * dΩ
    ju0(uh, duh, v) =∫(∇(duh) ⋅ ∇(v)) * dΩ
    ju0u(uh, duh, v, u_) =∫(∇(duh) ⋅ ∇(v) - (dΦ ∘ u_) ⋅ duh ⋅ v) * dΩ

    asn((du, w), (v, q), R) = ∫(
            R ⋅ v
    )*dΩ
    jsn((du, w), (ddu, dw), (v, q), dΦu_)= ∫(
        ddu ⋅ v - dw ⋅ v
        + ∇(dw) ⋅ ∇(q) - dΦu_ ⋅ ddu ⋅ q
        # + ∇(dw) ⋅ ∇(q) - dΦu_ ⋅ ∇(ddu) ⋅ ∇(q) -  ddu ⋅ ∇(dΦu_) ⋅  ∇(q)
        # ddu ⋅ v - dξ ⋅ v - dw ⋅ v
        # + (dΦ ∘ uh) ⋅ ddu ⋅ ζ - dξ ⋅ ζ
        # + ∇(dw) ⋅ ∇(q) + ∇(dξ) ⋅ ∇(q)
    )*dΩ   

    asnp() = error("Not implemented.")
    jsnp() = error("Not implemented.")
    jmy() = error("Not implemented.")
    # asnp((du, pu, ξ, w), (v, pv, ζ, q), R) = ∫(
    #         R ⋅ v
    # )*dΩ
    # jsnp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q), uh, puh, Th, ur)= ∫(
    #     ddu ⋅ v - ψ ⋅ dξ ⋅ v - dw ⋅ v
    #     + _dpB(dΩ, ddu, uh, ur) ⋅ pv - dpu ⋅ pv
    #     + AΦ(dξ, ζ, puh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - puh)) ⋅ dpu ⋅ ζ
    #     + ∇(dw) ⋅ ∇(q) + ψ ⋅ ∇(dξ) ⋅ ∇(q) + dξ ⋅ ∇(ψ) ⋅ ∇(q)
    # )*dΩ

    # AS(w, q, uh, Th, ds) = ∇(w) ⋅ ∇(q) + (ds ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ w ⋅ q
    # jmy((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uh, Th, ds)= ∫(
    #     ddu ⋅ v - dw ⋅ v
    #     + AΦ(dξ, ζ, uh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ ddu ⋅ ζ
    #     + AS(dw, q, uh, Th, ds) - (ds ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ ϕ ⋅ dξ ⋅ q
    # )*dΩ

    Vu = Uu.space
    NonlinearVI{Float64}(dΩ, Φ, dΦ, Φ₀, f, (au, ju, juu), (au0, ju0, ju0u), (asn, jsn), (asnp, jsnp), (asn, jmy), (Uu, Vu))
end

function projectionB(VI::NonlinearVI{T}, u, proj_rc; show_trace=true) where T
    Uu = first(VI.fe_space_u)
    r, c = proj_rc
    u_norm = Xnorm(VI, u-c)
    if u_norm ≤ r
        show_trace && print("No projection\n")
        return u
    else
        show_trace && print("Projection, u_norm=$u_norm\n")
        return FEFunction(Uu, r * (u.free_values[:] .- c) / u_norm)
    end
end

function Moreau_Yosida_it(VI::NonlinearVI; u₀=[], uᵢ=[], ρ=1e-5, NL=false, bt=true, tol=1e-10, max_iter=400, show_trace=true)

    s(u) = σ(u,ρ)/ρ 
    ds(u) = dσ(u,ρ)/ρ

    au, ju, juu = VI.au
    Uu, Vu = VI.fe_space_u

    bu(uh, v) = NL ? au(uh, v, s, uh) : au(uh, v, s, uᵢ)
    jbu(uh, duh, v) = NL ? juu(uh, duh, v, ds, uh) : ju(uh, duh, v, ds)

    opu = FEOperator(bu, jbu, Uu, Vu)
    if bt == false
        uB, its_u = newton(opu, u₀, max_iter=max_iter, damping=1, tol=tol, info=true, show_trace=show_trace);
    else
        nls = NLSolver(show_trace=show_trace, method=:newton, linesearch=LineSearches.BackTracking(), ftol=tol, xtol=10*eps())
        solver = FESolver(nls)
        uB, its_u = solve!(u₀,solver,opu)
        its_u = its_u.result.iterations
    end
    return (uB, its_u)
end

function Path_Following_S(VI::NonlinearVI; u₀=[], ρ0=1, ρ_min=1e-6, NL=false, max_its=20, tol=1e-10, bt=true, show_trace=true)
    Uu = first(VI.fe_space_u)
    uh = FEFunction(Uu, u₀.free_values[:])
    its = 0
    for ρ in [ρ0*10.0^(-i) for i in 0:20]
        show_trace && print("\n Considering ρ = $ρ.\n")
        (uh, it) = Moreau_Yosida_it(VI, u₀=uh, uᵢ=FEFunction(Uu,u₀.free_values[:]), ρ=ρ, NL=NL, bt=bt, tol=tol, show_trace=show_trace)
        its += it
        if ρ ≤ ρ_min
            break
        end
    end
    return (uh, its)
end

function hik_S(VI::NonlinearVI; u₀=[], uᵢ=[], tol=1e-10, NL=false, show_trace=true)
    au0, ju0, ju0u = VI.au0
    Uu, Vu = VI.fe_space_u
    bu(uh, v) = NL ? au0(uh, v, uh) : au0(uh, v, uᵢ)
    ju(uh, duh, v) = NL ? ju0u(uh, duh, v, uh) : ju0(uh, duh, v)

    opu = FEOperator(bu, ju, Uu, Vu)
    lb = -1e10*ones(Vu.nfree)
    ub = VI.Φ₀.free_values
    uB, its_u = hik(opu, u₀, lb, ub, max_iter=800, damping=1, tol=tol, info=true, show_trace=show_trace);
    return (uB, its_u)
end

function inner_solve(VI::NonlinearVI, u, proj_rc::Tuple{Number, Number}, tol::Number, hik_tol::Number, bt::Bool, PF::Bool, FS::Bool, ρ0::Number, ρ_min::Number, NL::Bool, newton_its::Int, pf_its::Int, hik_its::Int; show_trace=true)
    show_trace && print("\n   Project u.\n")
    pu = projectionB(VI, u, proj_rc, show_trace=show_trace)

    if PF==true
        show_trace && print("\n   Path-following MY for u.\n")
        (S, it) = Path_Following_S(VI, u₀=pu, tol=tol, bt=bt, ρ0=ρ0, ρ_min=ρ_min, NL=NL, show_trace=show_trace); pf_its+=it;
    else
        S = pu
    end

    if FS==true
        show_trace && print("\n   HIK feasibility step for u.\n")
        (S, it) = hik_S(VI, u₀=S, uᵢ=pu, tol=hik_tol, NL=NL, show_trace=show_trace); hik_its+=it;
    end

    return pu, S, newton_its, pf_its, hik_its
end


h1(VI::NonlinearVI, u, v) = sqrt(sum(∫((u-v) ⋅ (u-v) + ∇(u-v) ⋅ ∇(u-v))*VI.dΩ))

function visolver(VI::NonlinearVI, uᵢ; max_its=20, min_its=0, tol=1e-10, hik_tol=1e-10, proj_rc=(Inf, 0), bt=true, PF=true, FS=true, ρ0=1, ρ_min=1e-6, show_trace=true, show_inner_trace=true)

    Uu, Vu = VI.fe_space_u
    h1_1, zhs = [], []
    newton_its, pf_its, hik_its, outer_its =  0, 0, 0, 0
    append!(zhs, [uᵢ])
    while outer_its < max_its

        puᵢ, uB, newton_its, pf_its, hik_its = inner_solve(VI, uᵢ, proj_rc, tol, hik_tol, bt, PF, FS, ρ0, ρ_min, true, newton_its, pf_its, hik_its, show_trace=show_inner_trace)

        # TODO: is this the correct error to record?
        append!(h1_1, h1(VI, uB, uᵢ)+h1(VI, puᵢ, uᵢ))
        uᵢ = FEFunction(Vu, uB.free_values[:])
  
        outer_its += 1
        show_trace && print("VI Solver: Iteration $outer_its, ‖uB - uᵢ‖ + ‖P ∘ uᵢ - uᵢ‖ = $(last(h1_1))\n")

        append!(zhs, [uᵢ])
        if last(h1_1) ≤ tol && outer_its ≥ min_its
            break
        end


    end
    its = (outer_its, newton_its, pf_its, hik_its)
    return (zhs, h1_1, its)
end

function fixed_point(VI::NonlinearVI, uᵢ; max_its=20, min_its=0, tol=1e-10, hik_tol=1e-10, proj_rc=(Inf, 0), bt=true, PF=true, FS=true, ρ0=1, ρ_min=1e-6, show_trace=true, show_inner_trace=true)

    Uu, Vu = VI.fe_space_u
    h1_1, zhs = [], []
    newton_its, pf_its, hik_its, outer_its =  0, 0, 0, 0
    append!(zhs, [uᵢ])
    while outer_its < max_its

        puᵢ, uB, newton_its, pf_its, hik_its = inner_solve(VI, uᵢ, proj_rc, tol, hik_tol, bt, PF, FS, ρ0, ρ_min, false, newton_its, pf_its, hik_its, show_trace=show_inner_trace)

        # TODO: is this the correct error to record?
        append!(h1_1, h1(VI, uB, uᵢ)+h1(VI, puᵢ, uᵢ))
        uᵢ = FEFunction(Vu, uB.free_values[:])
  
        outer_its += 1
        show_trace && print("Fixed point: Iteration $outer_its, ‖uB - uᵢ‖ + ‖P ∘ uᵢ - uᵢ‖ = $(last(h1_1))\n")

        append!(zhs, [uᵢ])
        if last(h1_1) ≤ tol && outer_its ≥ min_its
            break
        end


    end
    its = (outer_its, newton_its, pf_its, hik_its)
    return (zhs, h1_1, its)
end


function semismoothnewton(VI::NonlinearVI, uᵢ; max_its=10, tol=1e-10, hik_tol=1e-10, globalization=false, proj_rc=(Inf, 0.0), ρ_min=1e-6, bt=true, PF=true, FS=true, show_trace=true, show_inner_trace=true, X=[])

    FS == false && @warn("Are you sure you want FS=false? This will prevent superlinear convergence.")
    Uu, Vu = VI.fe_space_u

    asn, jsn = VI.asn
    asnp, jsnp = VI.asnp

    U, V = MultiFieldFESpace([Uu, Uu]), MultiFieldFESpace([Vu, Vu])
    Up, Vp =  MultiFieldFESpace([Uu, Uu, Uu]), MultiFieldFESpace([Vu, Vu, Vu])

    zhs, h1s = [], []
    append!(zhs, [uᵢ])
    h1c, outer_its, hik_its, pf_its, newton_its = 1,0,0,0,0
    is_proj, is_xBs = [], []

    puᵢ, uB, newton_its, pf_its, hik_its = inner_solve(VI, uᵢ, proj_rc, tol, hik_tol, bt, PF, FS, 1, ρ_min, false, newton_its, pf_its, hik_its, show_trace=show_inner_trace)
    h1c = h1(VI, uB, uᵢ)
    append!(h1s, h1c)

    r, c = proj_rc
    c != 0.0 && error("Projection not implemented for centre not at c=0.0, currently c=$c.")
    while outer_its < max_its && h1c > tol

        R = uᵢ - uB
        dΦu_ = interpolate_everywhere(x->VI.dΦ(uᵢ(Point(x...))), Uu)
        # dΦu_ = VI.dΦ ∘ uᵢ
        u_norm = Xnorm(VI, uᵢ)

        if u_norm ≤ r
            show_trace && print("‖u‖ ≤ r.  ")
            b((du, w), (v, q)) = asn((du, w), (v, q), R)
            jb((du, w), (ddu, dw), (v, q)) = jsn((du, w), (ddu, dw), (v, q), dΦu_)
            append!(is_proj, false)
        else
            error("Not implemented.")
            # print("‖u‖ > r.  ")
            # bp((du, pu, ξ, w), (v, pv, ζ, q)) = asnp((du, pu, ξ, w), (v, pv, ζ, q), R)
            # jbp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q)) = jsnp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q), uᵢ, puᵢ, Tᵢ, (u_norm, r))
            # append!(is_proj, true)
        end

        m, Ũ, Ṽ, b̃, jb̃ = u_norm ≤ r ? (1, U, V, b, jb) : (2, Up, Vp, bp, jbp)

        # show_trace && print("Semismooth Newton step.\n")
        Ah = interpolate_everywhere(VI.Φ₀ - uB, Vu)
        A = findall(Ah.free_values .≤ 0) .+ (m*Vu.nfree)
        inac = setdiff(1:((m+1)*Vu.nfree), A)

        op = FEOperator(b̃, jb̃, Ũ, Ṽ)
        zh = FEFunction(Ũ, [uᵢ.free_values[:]; zeros(m*Vu.nfree)])

        res, J  = Gridap.Algebra.residual_and_jacobian(op, zh);
        dz = zeros((m+1)*Vu.nfree)
        dz[inac] = -J[inac,inac] \ res[inac]
        dzh = FEFunction(zh.fe_space, dz)
        δuN = dzh.single_fe_functions[1]

        # τ = defl_τ(uᵢ.free_values[:], δuN.free_values[:], [zeros(Vu.nfree)], X)

        uN = FEFunction(Vu, uᵢ.free_values[:] + δuN.free_values[:])
        puN, SN, newton_its, pf_its, hik_its = inner_solve(VI, uN, proj_rc, tol, hik_tol, bt, PF, FS, 1e-2, ρ_min, false, newton_its, pf_its, hik_its, show_trace=show_inner_trace)

        # TODO: is this the correct error to be tracking?
        h1N = h1(VI, uN, SN) + h1(VI, uN, puN)

        if globalization == true
            puB, SB, newton_its, pf_its, hik_its = inner_solve(VI, uB, proj_rc, tol, hik_tol, bt, PF, FS, 1e-2, ρ_min, false, newton_its, pf_its, hik_its, show_trace=show_inner_trace)

            # TODO: is this the correct error to be tracking?
            h1B = h1(VI, uB, SB) + h1(VI, uB, puB)
            show_trace && h1B < h1N ? print("H¹ norms = ($h1B, $h1N), uB superior.\n") : print("H¹ norms = ($h1B, $h1N), uN superior.\n")

            is_xB = h1B < h1N ? true : false
            append!(is_xBs, is_xB)

            h1c, uᵢ, puᵢ, uB = is_xB ? (h1B, uB, puB, SB) : (h1N, uN, puN, SN)
        else
            h1c, uᵢ, puᵢ, uB = h1N, uN, puN, SN
        end

        append!(h1s, h1c)
        append!(zhs, [uᵢ])
        outer_its += 1
        show_trace && print("Semismooth Newton: Iteration $outer_its, ‖uᵢ₊₁ - uᵢ‖ + ‖P ∘ uᵢ - uᵢ‖ = $h1c\n")

    end
    append!(zhs, [uB])
    its = (outer_its, newton_its, pf_its, hik_its)
    is = (is_proj, is_xBs)
    return (zhs, h1s, its, is)
end

# function moreau_yosida_newton(VI::NonlinearVI, uᵢ, Tᵢ; ρ=1e-5, max_its=10, inner_max_its=400, tol=1e-10, globalization=false, proj_rc=(Inf, 0.0), bt=true, PF=true)
    
#     Uu, Vu = VI.fe_space_u
#     UT, VT = VI.fe_space_T

#     amy, jmy = VI.amy

#     U, V = MultiFieldFESpace([Uu, UT, Uu]), MultiFieldFESpace([Vu, VT, Vu])
    
#     zhs, h1s = [], []
#     h1c, outer_its, hik_its, pf_its, newton_its = 1,0,0,0,0
#     is_proj, is_xBs = [], []

#     # puᵢ, Tᵢ, uB, newton_its, pf_its, hik_its = inner_solve(VI, uᵢ, Tᵢ, proj_rc, tol, bt, PF, false, ρ, newton_its, pf_its, hik_its)
    
#     (Tᵢ, it) = Φ(VI, uᵢ, T₀=Tᵢ, bt=true, tol=tol); newton_its+=it;
#     uᵢ_ = FEFunction(Vu, uᵢ.free_values[:])
#     (uB, it) = Moreau_Yosida_it(VI, Tᵢ, u₀=uᵢ_, ρ=ρ, bt=true, tol=tol, max_iter=inner_max_its); pf_its+=it;

#     h1c = h1(VI, uB, uᵢ)
#     append!(h1s, h1c)

#     ds(u) = dσ(u,ρ)/ρ

#     while outer_its < max_its && h1c > tol
#         R = uᵢ - uB

#         u_norm = Xnorm(VI, uᵢ)
#         r, c = proj_rc

#         if u_norm ≤ r
#             print("\n|u| ≤ r.\n")
#             b((du, ξ, w), (v, ζ, q)) = amy((du, ξ, w), (v, ζ, q), R)
#             jb((du, ξ, w), (ddu, dξ, dw), (v, ζ, q)) = jmy((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uᵢ, Tᵢ, ds)
#             append!(is_proj, false)
#         else
#             error("Projection not implemented for Moreau-Yosida Newton.")
#         end

#         print("\nMoreau-Yosida Newton step.\n")
#         op = FEOperator(b, jb, U, V)
#         # zh = FEFunction(U, [uᵢ.free_values; Tᵢ.free_values; uᵢ.free_values])
#         # zh = FEFunction(U, [uᵢ.free_values; Tᵢ.free_values; uᵢ.free_values])
#         # zN = newton(op, zh, max_iter=100, damping=1, tol=1e-7);
#         # δuN = zN.single_fe_functions[1]
#         # uN = FEFunction(Vu, uᵢ.free_values + δuN.free_values)

#         zh = FEFunction(U, [uᵢ.free_values[:]; zeros(VT.nfree+Vu.nfree)])

#         res, J  = Gridap.Algebra.residual_and_jacobian(op, zh);
#         dz = -J \ res
#         dzh = FEFunction(zh.fe_space, dz)
#         δuN = dzh.single_fe_functions[1]
#         uN = FEFunction(Vu, uᵢ.free_values[:] + δuN.free_values[:])

#         # puN, TN, SN, newton_its, pf_its, hik_its = inner_solve(VI, uN, Tᵢ, proj_rc, tol, bt, PF, false, ρ, newton_its, pf_its, hik_its)
#         (TN, it) = Φ(VI, uN, T₀=Tᵢ, bt=true, tol=tol); newton_its+=it;
#         uN_ = FEFunction(Vu, uN.free_values[:])
#         (SN, it) = Moreau_Yosida_it(VI, TN, u₀=uN_, ρ=ρ, bt=true, tol=tol, max_iter=inner_max_its); pf_its+=it;
#         puN = uN

#         h1N = h1(VI, uN, SN)

#         if globalization == true
#             # puB, TB, SB, newton_its, pf_its, hik_its = inner_solve(VI, uB, Tᵢ, proj_rc, tol, bt, PF, false, ρ, newton_its, pf_its, hik_its)

#             (TB, it) = Φ(VI, uB, T₀=Tᵢ, bt=true, tol=tol); newton_its+=it;
#             uB_ = FEFunction(Vu, uB.free_values[:])
#             (SB, it) = Moreau_Yosida_it(VI, TB, u₀=uB_, ρ=ρ, bt=true, tol=tol, max_iter=inner_max_its); pf_its+=it;
#             puB = uB

#             h1B = h1(VI, uB, SB)
#             h1B < h1N ? print("H¹ norms = ($h1B, $h1N), uB superior.\n\n") : print("H¹ norms = ($h1B, $h1N), uN superior.\n\n")

#             is_xB = h1B < h1N ? true : false
#             append!(is_xBs, is_xB)

#             h1c, uᵢ, puᵢ, Tᵢ, uB = is_xB ? (h1B, uB, puB, TB, SB) : (h1N, uN, puN, TN, SN)
#         else
#             h1c, uᵢ, puᵢ, Tᵢ, uB = h1N, uN, puN, TN, SN
#         end

#         append!(h1s, h1c)
#         append!(zhs, [[uᵢ, Tᵢ]])
#         outer_its += 1

#     end
#     append!(zhs, [[uB, Tᵢ]])
#     its = (outer_its, newton_its, pf_its, hik_its)
#     is = (is_proj, is_xBs)
#     return (zhs, h1s, its, is)
# end

function EOC(VI::NonlinearVI, us::AbstractVector, u::Gridap.FESpaces.SingleFieldFEFunction)
    errs = zeros(length(us))
    for i in 1:lastindex(us)
        errs[i] = Xnorm(VI, us[i]-u) 
    end
    errs, log.( errs[3:end] ./ errs[2:end-1] ) ./ log.( errs[2:end-1] ./ errs[1:end-2] )
end