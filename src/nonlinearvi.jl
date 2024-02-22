"""
    NonlinearVI()

Given an obstacle φ, find u ∈ H¹₀(Ω) that satisfies
    u ≤ φ, ⟨-Δu - f - Φ(u), v - u⟩≥0 ∀ v∈H¹₀(Ω), v≤ φ.

"""
IN_TOL = 1e-10
struct NonlinearVI{T}
    dΩ::Gridap.CellData.GenericMeasure
    k::T
    Φ₀::Gridap.FESpaces.SingleFieldFEFunction
    ϕ::Gridap.FESpaces.SingleFieldFEFunction
    Ψ₀::Gridap.FESpaces.SingleFieldFEFunction
    ψ::Gridap.FESpaces.SingleFieldFEFunction
    g::Function
    dg::Function
    f::Gridap.FESpaces.SingleFieldFEFunction
    au::Tuple{Function, Function}
    au0::Tuple{Function, Function}
    aT::Tuple{Function, Function}
    asn::Tuple{Function, Function}
    asnp::Tuple{Function, Function}
    amy::Tuple{Function, Function}
    fe_space_u::Tuple{FESpace, FESpace}
    fe_space_T::Tuple{FESpace, FESpace}
end

Xinner(V::NonlinearVI{T}, u, v) where T = sum(∫(∇(u) ⋅ ∇(v))*V.dΩ)
Xnorm(V::NonlinearVI{T}, u) where T = sqrt(Xinner(V, u, u))
_inner_X(dΩ, u, h) = sum(∫(∇(u) ⋅ ∇(h))*dΩ)
_dpB = (dΩ, h, u, ur) -> ur[2]/ur[1] * ( h - ((h->_inner_X(dΩ, u, h)) ∘ h) ⋅ u /ur[1]^2)

function NonlinearVI(dΩ::Gridap.CellData.GenericMeasure, k::T, Φ₀::Gridap.FESpaces.SingleFieldFEFunction,
    ϕ::Gridap.FESpaces.SingleFieldFEFunction,  Ψ₀::Gridap.FESpaces.SingleFieldFEFunction,
    ψ::Gridap.FESpaces.SingleFieldFEFunction,  g::Function, dg::Function, f::Gridap.FESpaces.SingleFieldFEFunction, Uu::FESpace, UT::FESpace) where T

    au(uh, v, s, Th) =∫( ∇(uh) ⋅ ∇(v) +  (s ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ v - f ⋅ v) * dΩ
    ju(uh, duh, v, ds, Th) =∫( ∇(duh) ⋅ ∇(v) +  (ds ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ duh ⋅ v) * dΩ
    
    aT(Th, R, uh) = ∫( ∇(Th) ⋅ ∇(R) + k*Th ⋅ R  - (g ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ R) * dΩ
    jT(Th, dTh, R, uh) =∫( ∇(dTh) ⋅ ∇(R) +  k*dTh ⋅ R - (dg ∘ (Ψ₀ + ψ ⋅ Th - uh) ⋅ (ψ ⋅ dTh) ⋅ R) ) * dΩ
    
    au0(uh, v) =∫(∇(uh) ⋅ ∇(v) - f ⋅ v) * dΩ
    ju0(uh, duh, v) =∫(∇(duh) ⋅ ∇(v)) * dΩ

    AΦ(ξ, ζ, uh, Th) = ∇(ξ) ⋅ ∇(ζ) + k*ξ ⋅ ζ - (dg ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ (ψ ⋅ ξ) ⋅ ζ
    asn((du, ξ, w), (v, ζ, q), R) = ∫(
            R ⋅ v
    )*dΩ
    jsn((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uh, Th)= ∫(
        ddu ⋅ v - ψ ⋅ dξ ⋅ v - dw ⋅ v
        + AΦ(dξ, ζ, uh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ ddu ⋅ ζ
        + ∇(dw) ⋅ ∇(q) + ψ ⋅ ∇(dξ) ⋅ ∇(q) + dξ ⋅ ∇(ψ) ⋅ ∇(q)
    )*dΩ   

    asnp((du, pu, ξ, w), (v, pv, ζ, q), R) = ∫(
            R ⋅ v
    )*dΩ
    jsnp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q), uh, puh, Th, ur)= ∫(
        ddu ⋅ v - ψ ⋅ dξ ⋅ v - dw ⋅ v
        + _dpB(dΩ, ddu, uh, ur) ⋅ pv - dpu ⋅ pv
        + AΦ(dξ, ζ, puh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - puh)) ⋅ dpu ⋅ ζ
        + ∇(dw) ⋅ ∇(q) + ψ ⋅ ∇(dξ) ⋅ ∇(q) + dξ ⋅ ∇(ψ) ⋅ ∇(q)
    )*dΩ

    AS(w, q, uh, Th, ds) = ∇(w) ⋅ ∇(q) + (ds ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ w ⋅ q
    # amy((du, ξ, w), (v, ζ, q), R, uh, Th, ds)= ∫(
    #     du ⋅ v - w ⋅ v + R ⋅ v
    #     + AΦ(ξ, ζ, uh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ du ⋅ ζ
    #     + AS(w, q, uh, Th, ds) - (ds ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ ϕ ⋅ ξ ⋅ q
    # )*dΩ
    jmy((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uh, Th, ds)= ∫(
        ddu ⋅ v - dw ⋅ v
        + AΦ(dξ, ζ, uh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ ddu ⋅ ζ
        + AS(dw, q, uh, Th, ds) - (ds ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ ϕ ⋅ dξ ⋅ q
    )*dΩ

    Vu, VT = Uu.space, UT
    NonlinearVI{T}(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, (au, ju), (au0, ju0), (aT, jT), (asn, jsn), (asnp, jsnp), (asn, jmy), (Uu, Vu), (UT, VT))
end

# Smoothed Moreau-Yosida penalisation and its derivative
function σ(u, ρ)
    if u ≤ 0
        return 0.0
    elseif 0 < u < ρ
        return u^2/(2ρ)
    else
        return u - ρ/2
    end
end
function dσ(u, ρ)
    if u ≤ 0
        return 0.0
    elseif 0 < u < ρ
        return u/ρ
    else
        return 1.0
    end
end

function projectionB(V::NonlinearVI{T}, u, proj_rc; show_trace=true) where T
    Uu = first(V.fe_space_u)
    r, c = proj_rc
    u_norm = Xnorm(V, u-c)
    if u_norm ≤ r
        show_trace && print("No projection\n")
        return u
    else
        show_trace && print("Projection, u_norm=$u_norm\n")
        return FEFunction(Uu, r * (u.free_values[:] .- c) / u_norm)
    end
end

function Φ(V::NonlinearVI, uh; T₀=[], bt=true, tol=IN_TOL, show_trace=true)
    aT, jT = V.aT
    UT, VT = V.fe_space_T

    bT(Th, R) = aT(Th, R, uh)
    jbT(Th, dTh, R) = jT(Th, dTh, R, uh)
    opT = FEOperator(bT, jbT, UT, VT)
    if bt == false
        T, its_T = newton(opT, T₀, max_iter=400, damping=1, tol=IN_TOL, info=true, show_trace=show_trace)
    else
        nls = NLSolver(show_trace=show_trace, method=:newton, linesearch=LineSearches.BackTracking(), ftol=tol, xtol=10*eps())
        solver = FESolver(nls)
        T, its_T = solve!(FEFunction(VT, T₀.free_values[:]),solver,opT)
        its_T = its_T.result.iterations
    end
    return (T, its_T)
end

function Moreau_Yosida_it(V::NonlinearVI, T; u₀=[], ρ=1e-5, bt=true, tol=IN_TOL, max_iter=400, show_trace=true)

    s(u) = σ(u,ρ)/ρ 
    ds(u) = dσ(u,ρ)/ρ

    au, ju = V.au
    Uu, Vu = V.fe_space_u

    bu(uh, v) = au(uh, v, s, T)
    jbu(uh, duh, v) = ju(uh, duh, v, ds, T)
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

function Path_Following_S(V::NonlinearVI, Tᵢ; u₀=[], ρ0=1, max_its=20, tol=IN_TOL, bt=true, show_trace=true)
    Uu = first(V.fe_space_u)
    uh = FEFunction(Uu, u₀.free_values[:])
    its = 0
    for ρ in [ρ0*10.0^(-i) for i in 0:5]
        show_trace && print("\n Considering ρ = $ρ.\n")
        (uh, it) = Moreau_Yosida_it(V, Tᵢ, u₀=uh, ρ=ρ, bt=bt, tol=tol, show_trace=show_trace)
        its += it
        if ρ < 1e-6
            break
        end
    end
    return (uh, its)
end

function bm_S(V::NonlinearVI,Tᵢ; u₀=[],show_trace=true)
    au0, ju0 = V.au0
    Uu, Vu = V.fe_space_u
    opu = FEOperator(au0, ju0, Uu, Vu)
    lb = -1e10*ones(Vu.nfree)
    ub = interpolate_everywhere(V.Φ₀ + V.ϕ ⋅ Tᵢ, Uu).free_values
    uB, its_u = bm(opu, u₀, lb, ub, max_iter=800, damping=1, tol=IN_TOL, info=true, show_trace=show_trace);
    return (uB, its_u)
end

function hik_S(V::NonlinearVI,Tᵢ; u₀=[], show_trace=true)
    au0, ju0 = V.au0
    Uu, Vu = V.fe_space_u
    opu = FEOperator(au0, ju0, Uu, Vu)
    lb = -1e10*ones(Vu.nfree)
    ub = interpolate_everywhere(V.Φ₀ + V.ϕ ⋅ Tᵢ, Uu).free_values
    uB, its_u = hik(opu, u₀, lb, ub, max_iter=800, damping=1, tol=IN_TOL, info=true, show_trace=show_trace);
    return (uB, its_u)
end

function inner_solve(V::NonlinearVI, u, T_, proj_rc::Tuple{Number, Number}, tol::Number, bt::Bool, PF::Bool, FS::Bool, ρ0::Number, newton_its::Int, pf_its::Int, hik_its::Int; show_trace=true)
    show_trace && print("\n   Project u.\n")
    pu = projectionB(V, u, proj_rc, show_trace=show_trace)
    show_trace && print("\n   Solve for T.\n")
    (T, it) = Φ(V, pu, T₀=T_,tol=tol, bt=bt, show_trace=show_trace); newton_its+=it;

    if PF==true
        show_trace && print("\n   Path-following MY for u.\n")
        (S, it) = Path_Following_S(V, T, u₀=pu, tol=IN_TOL, bt=bt, ρ0=ρ0, show_trace=show_trace); pf_its+=it;
    else
        S = pu
    end

    if FS==true
        show_trace && print("\n   HIK feasibility step for u.\n")
        (S, it) = hik_S(V, T, u₀=S, show_trace=show_trace); hik_its+=it;
    end

    return pu, T, S, newton_its, pf_its, hik_its
end


h1(V::NonlinearVI, u, v) = sqrt(sum(∫((u-v) ⋅ (u-v) + ∇(u-v) ⋅ ∇(u-v))*V.dΩ))

function fixed_point(V::NonlinearVI, uᵢ, Tᵢ; max_its=20, min_its=0, tol=IN_TOL, proj_rc=(Inf, 0), bt=true, PF=true, ρ0=1, show_trace=true, show_inner_trace=true)

    Uu, Vu = V.fe_space_u
    h1_1, zhs = [], []
    newton_its, pf_its, hik_its, outer_its =  0, 0, 0, 0
    append!(zhs, [[uᵢ, Tᵢ]])
    while outer_its < max_its

        puᵢ, Tᵢ, uB, newton_its, pf_its, hik_its = inner_solve(V, uᵢ, Tᵢ, proj_rc, tol, bt, PF, true, ρ0, newton_its, pf_its, hik_its, show_trace=show_inner_trace)

        # TODO: is this the correct error to record?
        append!(h1_1, h1(V, uB, uᵢ)+h1(V, puᵢ, uᵢ))
        uᵢ = FEFunction(Vu, uB.free_values[:])
  
        outer_its += 1
        show_trace && print("Fixed point: Iteration $outer_its, ‖uB - uᵢ‖ + ‖P ∘ uᵢ - uᵢ‖ = $(last(h1_1))\n")

        append!(zhs, [[uᵢ, Tᵢ]])
        if last(h1_1) ≤ tol && outer_its ≥ min_its
            break
        end


    end
    its = (outer_its, newton_its, pf_its, hik_its)
    return (zhs, h1_1, its)
end

function semismoothnewton(V::NonlinearVI, uᵢ, Tᵢ; max_its=10, tol=IN_TOL, globalization=false, proj_rc=(Inf, 0.0), bt=true, PF=true, FS=true, show_trace=true, show_inner_trace=true, X=[])

    FS == false && @warn("Are you sure you want FS=false? This will prevent superlinear convergence.")
    Uu, Vu = V.fe_space_u
    UT, VT = V.fe_space_T

    asn, jsn = V.asn
    asnp, jsnp = V.asnp

    U, V = MultiFieldFESpace([Uu, UT, Uu]), MultiFieldFESpace([Vu, VT, Vu])
    Up, Vp =  MultiFieldFESpace([Uu, Uu, UT, Uu]), MultiFieldFESpace([Vu, Vu, VT, Vu])

    zhs, h1s = [], []
    append!(zhs, [[uᵢ, Tᵢ]])
    h1c, outer_its, hik_its, pf_its, newton_its = 1,0,0,0,0
    is_proj, is_xBs = [], []

    puᵢ, Tᵢ, uB, newton_its, pf_its, hik_its = inner_solve(V, uᵢ, Tᵢ, proj_rc, tol, bt, PF, FS, 1, newton_its, pf_its, hik_its, show_trace=show_inner_trace)
    h1c = h1(V, uB, uᵢ)
    append!(h1s, h1c)

    r, c = proj_rc
    c != 0.0 && error("Projection not implemented for centre not at c=0.0, currently c=$c.")
    while outer_its < max_its && h1c > tol

        R = uᵢ - uB

        u_norm = Xnorm(V, uᵢ)

        if u_norm ≤ r
            print("‖u‖ ≤ r.  ")
            b((du, ξ, w), (v, ζ, q)) = asn((du, ξ, w), (v, ζ, q), R)
            jb((du, ξ, w), (ddu, dξ, dw), (v, ζ, q)) = jsn((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uᵢ, Tᵢ)
            append!(is_proj, false)
        else
            print("‖u‖ > r.  ")
            bp((du, pu, ξ, w), (v, pv, ζ, q)) = asnp((du, pu, ξ, w), (v, pv, ζ, q), R)
            jbp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q)) = jsnp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q), uᵢ, puᵢ, Tᵢ, (u_norm, r))
            append!(is_proj, true)
        end

        m, Ũ, Ṽ, b̃, jb̃ = u_norm ≤ r ? (1, U, V, b, jb) : (2, Up, Vp, bp, jbp)

        # show_trace && print("Semismooth Newton step.\n")
        Ah = interpolate_everywhere(V.Φ₀ + V.ϕ ⋅ Tᵢ - uB, Vu)
        A = findall(Ah.free_values .≤ 0) .+ (m*Vu.nfree + VT.nfree)
        inac = setdiff(1:(VT.nfree+(m+1)*Vu.nfree), A)

        op = FEOperator(b̃, jb̃, Ũ, Ṽ)
        zh = FEFunction(Ũ, [uᵢ.free_values[:]; zeros(VT.nfree+m*Vu.nfree)])

        res, J  = Gridap.Algebra.residual_and_jacobian(op, zh);
        dz = zeros(VT.nfree+(m+1)*Vu.nfree)
        dz[inac] = -J[inac,inac] \ res[inac]
        dzh = FEFunction(zh.fe_space, dz)
        δuN = dzh.single_fe_functions[1]

        # τ = defl_τ(uᵢ.free_values[:], δuN.free_values[:], [zeros(Vu.nfree)], X)

        uN = FEFunction(Vu, uᵢ.free_values[:] + δuN.free_values[:])

        puN, TN, SN, newton_its, pf_its, hik_its = inner_solve(V, uN, Tᵢ, proj_rc, tol, bt, PF, FS, 1e-2, newton_its, pf_its, hik_its, show_trace=show_inner_trace)

        # TODO: is this the correct error to be tracking?
        h1N = h1(V, uN, SN) + h1(V, uN, puN)

        if globalization == true
            puB, TB, SB, newton_its, pf_its, hik_its = inner_solve(V, uB, Tᵢ, proj_rc, tol, bt, PF, FS, 1e-2, newton_its, pf_its, hik_its, show_trace=show_inner_trace)

            # TODO: is this the correct error to be tracking?
            h1B = h1(V, uB, SB) + h1(V, uB, puB)
            show_trace && h1B < h1N ? print("H¹ norms = ($h1B, $h1N), uB superior.\n") : print("H¹ norms = ($h1B, $h1N), uN superior.\n")

            is_xB = h1B < h1N ? true : false
            append!(is_xBs, is_xB)

            h1c, uᵢ, puᵢ, Tᵢ, uB = is_xB ? (h1B, uB, puB, TB, SB) : (h1N, uN, puN, TN, SN)
        else
            h1c, uᵢ, puᵢ, Tᵢ, uB = h1N, uN, puN, TN, SN
        end

        append!(h1s, h1c)
        append!(zhs, [[uᵢ, Tᵢ]])
        outer_its += 1
        show_trace && print("Semismooth Newton: Iteration $outer_its, ‖uᵢ₊₁ - uᵢ‖ + ‖P ∘ uᵢ - uᵢ‖ = $h1c\n")

    end
    append!(zhs, [[uB, Tᵢ]])
    its = (outer_its, newton_its, pf_its, hik_its)
    is = (is_proj, is_xBs)
    return (zhs, h1s, its, is)
end

function moreau_yosida_newton(V::NonlinearVI, uᵢ, Tᵢ; ρ=1e-5, max_its=10, inner_max_its=400, tol=IN_TOL, globalization=false, proj_rc=(Inf, 0.0), bt=true, PF=true)
    
    Uu, Vu = V.fe_space_u
    UT, VT = V.fe_space_T

    amy, jmy = V.amy

    U, V = MultiFieldFESpace([Uu, UT, Uu]), MultiFieldFESpace([Vu, VT, Vu])
    
    zhs, h1s = [], []
    h1c, outer_its, hik_its, pf_its, newton_its = 1,0,0,0,0
    is_proj, is_xBs = [], []

    # puᵢ, Tᵢ, uB, newton_its, pf_its, hik_its = inner_solve(V, uᵢ, Tᵢ, proj_rc, tol, bt, PF, false, ρ, newton_its, pf_its, hik_its)
    
    (Tᵢ, it) = Φ(V, uᵢ, T₀=Tᵢ, bt=true, tol=tol); newton_its+=it;
    uᵢ_ = FEFunction(Vu, uᵢ.free_values[:])
    (uB, it) = Moreau_Yosida_it(V, Tᵢ, u₀=uᵢ_, ρ=ρ, bt=true, tol=tol, max_iter=inner_max_its); pf_its+=it;

    h1c = h1(V, uB, uᵢ)
    append!(h1s, h1c)

    ds(u) = dσ(u,ρ)/ρ

    while outer_its < max_its && h1c > tol
        R = uᵢ - uB

        u_norm = Xnorm(V, uᵢ)
        r, c = proj_rc

        if u_norm ≤ r
            print("\n|u| ≤ r.\n")
            b((du, ξ, w), (v, ζ, q)) = amy((du, ξ, w), (v, ζ, q), R)
            jb((du, ξ, w), (ddu, dξ, dw), (v, ζ, q)) = jmy((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uᵢ, Tᵢ, ds)
            append!(is_proj, false)
        else
            error("Projection not implemented for Moreau-Yosida Newton.")
        end

        print("\nMoreau-Yosida Newton step.\n")
        op = FEOperator(b, jb, U, V)
        # zh = FEFunction(U, [uᵢ.free_values; Tᵢ.free_values; uᵢ.free_values])
        # zh = FEFunction(U, [uᵢ.free_values; Tᵢ.free_values; uᵢ.free_values])
        # zN = newton(op, zh, max_iter=100, damping=1, tol=1e-7);
        # δuN = zN.single_fe_functions[1]
        # uN = FEFunction(Vu, uᵢ.free_values + δuN.free_values)

        zh = FEFunction(U, [uᵢ.free_values[:]; zeros(VT.nfree+Vu.nfree)])

        res, J  = Gridap.Algebra.residual_and_jacobian(op, zh);
        dz = -J \ res
        dzh = FEFunction(zh.fe_space, dz)
        δuN = dzh.single_fe_functions[1]
        uN = FEFunction(Vu, uᵢ.free_values[:] + δuN.free_values[:])

        # puN, TN, SN, newton_its, pf_its, hik_its = inner_solve(V, uN, Tᵢ, proj_rc, tol, bt, PF, false, ρ, newton_its, pf_its, hik_its)
        (TN, it) = Φ(V, uN, T₀=Tᵢ, bt=true, tol=tol); newton_its+=it;
        uN_ = FEFunction(Vu, uN.free_values[:])
        (SN, it) = Moreau_Yosida_it(V, TN, u₀=uN_, ρ=ρ, bt=true, tol=tol, max_iter=inner_max_its); pf_its+=it;
        puN = uN

        h1N = h1(V, uN, SN)

        if globalization == true
            # puB, TB, SB, newton_its, pf_its, hik_its = inner_solve(V, uB, Tᵢ, proj_rc, tol, bt, PF, false, ρ, newton_its, pf_its, hik_its)

            (TB, it) = Φ(V, uB, T₀=Tᵢ, bt=true, tol=tol); newton_its+=it;
            uB_ = FEFunction(Vu, uB.free_values[:])
            (SB, it) = Moreau_Yosida_it(V, TB, u₀=uB_, ρ=ρ, bt=true, tol=tol, max_iter=inner_max_its); pf_its+=it;
            puB = uB

            h1B = h1(V, uB, SB)
            h1B < h1N ? print("H¹ norms = ($h1B, $h1N), uB superior.\n\n") : print("H¹ norms = ($h1B, $h1N), uN superior.\n\n")

            is_xB = h1B < h1N ? true : false
            append!(is_xBs, is_xB)

            h1c, uᵢ, puᵢ, Tᵢ, uB = is_xB ? (h1B, uB, puB, TB, SB) : (h1N, uN, puN, TN, SN)
        else
            h1c, uᵢ, puᵢ, Tᵢ, uB = h1N, uN, puN, TN, SN
        end

        append!(h1s, h1c)
        append!(zhs, [[uᵢ, Tᵢ]])
        outer_its += 1

    end
    append!(zhs, [[uB, Tᵢ]])
    its = (outer_its, newton_its, pf_its, hik_its)
    is = (is_proj, is_xBs)
    return (zhs, h1s, its, is)
end

function EOC(V::NonlinearVI, us::AbstractVector{<:Gridap.FESpaces.SingleFieldFEFunction}, u::Gridap.FESpaces.SingleFieldFEFunction)
    errs = zeros(length(us))
    for i in 1:lastindex(us)
        errs[i] = Xnorm(V, us[i]-u) 
    end
    errs, log.( errs[3:end] ./ errs[2:end-1] ) ./ log.( errs[2:end-1] ./ errs[1:end-2] )
end