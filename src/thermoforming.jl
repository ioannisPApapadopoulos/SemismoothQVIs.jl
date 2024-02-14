"""
    GeneralizedThermoformingQVI()

Find u ∈ H¹₀(Ω) that satisfies
    u ≤ Φ₀ + Φ(u), ⟨-Δu - f, v - u⟩≥0 ∀ v∈H¹₀(Ω), v≤ Φ₀ + Φ(u)
with Φ(u) given by Φ(u) := ϕT and T as the solution of
    kT - ΔT = g(Ψ₀ + ψT - u), ∂ν T = 0 on ∂Ω.

"""
IN_TOL = 1e-10
struct GeneralizedThermoformingQVI{T}
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
    fe_space_u::Tuple{FESpace, FESpace}
    fe_space_T::Tuple{FESpace, FESpace}
end

_inner_X(dΩ, u, h) = sum(∫(∇(u) ⋅ ∇(h) + u ⋅ h)*dΩ)
_dpB = (dΩ, h, u, ur) -> ur[2]/ur[1] * ( h - ((h->_inner_X(dΩ, u, h)) ∘ h) ⋅ u /ur[1])

function GeneralizedThermoformingQVI(dΩ::Gridap.CellData.GenericMeasure, k::T, Φ₀::Gridap.FESpaces.SingleFieldFEFunction,
    ϕ::Gridap.FESpaces.SingleFieldFEFunction,  Ψ₀::Gridap.FESpaces.SingleFieldFEFunction,
    ψ::Gridap.FESpaces.SingleFieldFEFunction,  g::Function, dg::Function, f::Gridap.FESpaces.SingleFieldFEFunction, Uu::FESpace, UT::FESpace) where T

    au(uh, v, s, Th) =∫( ∇(uh) ⋅ ∇(v) +  (s ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ v - f ⋅ v) * dΩ
    ju(uh, duh, v, ds, Th) =∫( ∇(duh) ⋅ ∇(v) +  (ds ∘ (uh - Φ₀ - ϕ ⋅ Th)) ⋅ duh ⋅ v) * dΩ
    
    aT(Th, R, uh) = ∫( ∇(Th) ⋅ ∇(R) + k*Th ⋅ R  - (g ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ R) * dΩ
    jT(Th, dTh, R, uh) =∫( ∇(dTh) ⋅ ∇(R) +  k*dTh ⋅ R - (dg ∘ (Ψ₀ + ψ ⋅ Th - uh) ⋅ (ψ ⋅ dTh) ⋅ R) ) * dΩ
    
    au0(uh, v) =∫(∇(uh) ⋅ ∇(v) - f ⋅ v) * dΩ
    ju0(uh, duh, v) =∫(∇(duh) ⋅ ∇(v)) * dΩ

    AΦ(ξ, ζ, uh, Th) = ∇(ξ) ⋅ ∇(ζ) + k*ξ ⋅ ζ - (dg ∘ (ψ ⋅ Th - uh)) ⋅ (ψ ⋅ ξ) ⋅ ζ
    asn((du, ξ, w), (v, ζ, q), R) = ∫(
            R ⋅ v
    )*dΩ
    jsn((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uh, Th)= ∫(
        ddu ⋅ v - ψ ⋅ dξ ⋅ v - dw ⋅ v
        + AΦ(dξ, ζ, uh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - uh)) ⋅ ddu ⋅ ζ
        + ∇(dw) ⋅ ∇(q) + ϕ ⋅ ∇(dξ) ⋅ ∇(q) + dξ ⋅ ∇(ϕ) ⋅ ∇(q)
    )*dΩ   

    asnp((du, pu, ξ, w), (v, pv, ζ, q), R) = ∫(
            R ⋅ v
    )*dΩ
    jsnp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q), uh, puh, Th, ur)= ∫(
        ddu ⋅ v - ψ ⋅ dξ ⋅ v - dw ⋅ v
        + _dpB(dΩ, ddu, uh, ur) ⋅ pv - dpu ⋅ pv
        + AΦ(dξ, ζ, puh, Th) + (dg ∘ (Ψ₀ + ψ ⋅ Th - puh)) ⋅ dpu ⋅ ζ
        + ∇(dw) ⋅ ∇(q) + ϕ ⋅ ∇(dξ) ⋅ ∇(q) + dξ ⋅ ∇(ϕ) ⋅ ∇(q)
    )*dΩ

    Vu, VT = Uu.space, UT
    GeneralizedThermoformingQVI{T}(dΩ, k, Φ₀, ϕ, Ψ₀, ψ, g, dg, f, (au, ju), (au0, ju0), (aT, jT), (asn, jsn), (asnp, jsnp), (Uu, Vu), (UT, VT))
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

Xinner(Q::GeneralizedThermoformingQVI{T}, u, v) where T = sum(∫(∇(u) ⋅ ∇(v) + u ⋅ v)*Q.dΩ)
Xnorm(Q::GeneralizedThermoformingQVI{T}, u) where T = sqrt(Xinner(Q, u, u))

function projectionB(Q::GeneralizedThermoformingQVI{T}, u, proj_rc) where T
    Uu = first(Q.fe_space_u)
    r, c = proj_rc
    u_norm = Xnorm(Q, u-c)
    if u_norm ≤ r
        print("No projection\n")
        return u
    else
        print("Projection, u_norm=$u_norm\n")
        return FEFunction(Uu, r * (u.free_values[:] .- c) / u_norm)
    end
end

function Φ(Q::GeneralizedThermoformingQVI, uh; Tic=[], bt=true, tol=IN_TOL)
    aT, jT = Q.aT
    UT, VT = Q.fe_space_T

    bT(Th, R) = aT(Th, R, uh)
    jbT(Th, dTh, R) = jT(Th, dTh, R, uh)
    opT = FEOperator(bT, jbT, UT, VT)
    if bt == false
        T, its_T = newton(opT, Tic, max_iter=400, damping=1, tol=IN_TOL, info=true)
    else
        nls = NLSolver(show_trace=true, method=:newton, linesearch=LineSearches.BackTracking(), ftol=tol, xtol=10*eps())
        solver = FESolver(nls)
        T, its_T = solve!(FEFunction(VT, Tic.free_values[:]),solver,opT)
        its_T = its_T.result.iterations
    end
    return (T, its_T)
end

function Moreau_Yosida_it(Q::GeneralizedThermoformingQVI, T; uic=[], ρ=1e-5, bt=true, tol=IN_TOL, max_iter=400)

    s(u) = σ(u,ρ)/ρ 
    ds(u) = dσ(u,ρ)/ρ

    au, ju = Q.au
    Uu, Vu = Q.fe_space_u

    bu(uh, v) = au(uh, v, s, T)
    jbu(uh, duh, v) = ju(uh, duh, v, ds, T)
    opu = FEOperator(bu, jbu, Uu, Vu)
    if bt == false
        uB, its_u = newton(opu, uic, max_iter=max_iter, damping=1, tol=tol, info=true);
    else
        nls = NLSolver(show_trace=true, method=:newton, linesearch=LineSearches.BackTracking(), ftol=tol, xtol=10*eps())
        solver = FESolver(nls)
        uB, its_u = solve!(uic,solver,opu)
        its_u = its_u.result.iterations
    end
    return (uB, its_u)
end

function Path_Following_S(Q::GeneralizedThermoformingQVI, Tᵢ; uic=[], ρ0=1, max_its=20, tol=IN_TOL, bt=true)
    Uu = first(Q.fe_space_u)
    uh = FEFunction(Uu, uic.free_values[:])
    its = 0
    for ρ in [ρ0*10.0^(-i) for i in 0:5]
        print("\n Considering ρ = $ρ.\n")
        (uh, it) = Moreau_Yosida_it(Q, Tᵢ, uic=uh, ρ=ρ, bt=bt, tol=tol)
        its += it
        if ρ < 1e-6
            break
        end
    end
    return (uh, its)
end

function bm_S(Q::GeneralizedThermoformingQVI,Tᵢ; uic=[])
    au0, ju0 = Q.au0
    Uu, Vu = Q.fe_space_u
    opu = FEOperator(au0, ju0, Uu, Vu)
    lb = -1e10*ones(Vu.nfree)
    ub = interpolate_everywhere(Q.Φ₀ + Q.ϕ ⋅ Tᵢ, Uu).free_values
    uB, its_u = bm(opu, uic, lb, ub, max_iter=800, damping=1, tol=IN_TOL, info=true);
    return (uB, its_u)
end

function hik_S(Q::GeneralizedThermoformingQVI,Tᵢ; uic=[])
    au0, ju0 = Q.au0
    Uu, Vu = Q.fe_space_u
    opu = FEOperator(au0, ju0, Uu, Vu)
    lb = -1e10*ones(Vu.nfree)
    ub = interpolate_everywhere(Q.Φ₀ + Q.ϕ ⋅ Tᵢ, Uu).free_values
    uB, its_u = hik(opu, uic, lb, ub, max_iter=800, damping=1, tol=IN_TOL, info=true);
    return (uB, its_u)
end

function inner_solve(Q::GeneralizedThermoformingQVI, u, T_, proj_rc::Tuple{Number, Number}, tol::Number, bt::Bool, PF::Bool, ρ0::Number, newton_its::Int, pf_its::Int, hik_its::Int)
    print("\n   Project u.\n")
    pu = projectionB(Q, u, proj_rc)
    print("\n   Solve for T.\n")
    (T, it) = Φ(Q, pu, Tic=T_,tol=tol, bt=bt); newton_its+=it;

    print("\n   Path-following MY for u.\n")
    if PF==true
        (S, it) = Path_Following_S(Q, T, uic=pu, tol=tol, bt=bt, ρ0=ρ0); pf_its+=it;
    else
        S = pu
    end

    print("\n   HIK feasibility step for u.\n")
    (S, it) = hik_S(Q, T, uic=S); hik_its+=it;

    return pu, T, S, newton_its, pf_its, hik_its
end


h1(Q::GeneralizedThermoformingQVI, u, v) = sqrt(sum(∫((u-v) ⋅ (u-v) + ∇(u-v) ⋅ ∇(u-v))*Q.dΩ))

function fixed_point(Q::GeneralizedThermoformingQVI, uᵢ, Tᵢ; max_its=20, min_its=0, tol=IN_TOL, proj_rc=(Inf, 0), bt=true, PF=true, ρ0=1)

    Uu, Vu = Q.fe_space_u
    h1_1, zhs = [], []
    newton_its, pf_its, hik_its, outer_its =  0, 0, 0, 0
    while outer_its < max_its

        _, Tᵢ, uB, newton_its, pf_its, hik_its = inner_solve(Q, uᵢ, Tᵢ, proj_rc, tol, bt, PF, ρ0, newton_its, pf_its, hik_its)

        append!(h1_1, h1(Q, uB, uᵢ))
        uᵢ = FEFunction(Vu, uB.free_values[:])
  
        outer_its += 1
        append!(zhs, [[uᵢ, Tᵢ]])
        if last(h1_1) ≤ tol && outer_its ≥ min_its
            break
        end


    end
    its = (outer_its, newton_its, pf_its, hik_its)
    return (zhs, h1_1, its)
end

function newtonss(Q::GeneralizedThermoformingQVI, uᵢ, Tᵢ; max_its=10, tol=IN_TOL, globalization=false, proj_rc=(Inf, 0.0), bt=true, PF=true)

    Uu, Vu = Q.fe_space_u
    UT, VT = Q.fe_space_T

    asn, jsn = Q.asn
    asnp, jsnp = Q.asnp

    U, V = MultiFieldFESpace([Uu, UT, Uu]), MultiFieldFESpace([Vu, VT, Vu])
    Up, Vp =  MultiFieldFESpace([Uu, Uu, UT, Uu]), MultiFieldFESpace([Vu, Vu, VT, Vu])

    zhs, h1s = [], []
    h1c = 1
    outer_its, hik_its, pf_its, newton_its = 0,0,0,0
    is_proj, is_xBs = [], []

    puᵢ, Tᵢ, uB, newton_its, pf_its, hik_its = inner_solve(Q, uᵢ, Tᵢ, proj_rc, tol, bt, PF, 1, newton_its, pf_its, hik_its)
    h1c = h1(Q, uB, uᵢ)
    append!(h1s, h1c)


    while outer_its < max_its && h1c > tol

        R = uᵢ - uB

        u_norm = Xnorm(Q, uᵢ)
        r, c = proj_rc

        if u_norm ≤ r
            print("\n|u| ≤ r.\n")
            b((du, ξ, w), (v, ζ, q)) = asn((du, ξ, w), (v, ζ, q), R)
            jb((du, ξ, w), (ddu, dξ, dw), (v, ζ, q)) = jsn((du, ξ, w), (ddu, dξ, dw), (v, ζ, q), uᵢ, Tᵢ)
            append!(is_proj, false)
        else
            print("\n|u| > r.\n")
            bp((du, pu, ξ, w), (v, pv, ζ, q)) = asnp((du, pu, ξ, w), (v, pv, ζ, q), R)
            jbp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q)) = jsnp((du, pu, ξ, w), (ddu, dpu, dξ, dw), (v, pv, ζ, q), uᵢ, puᵢ, Tᵢ, (u_norm, r))
            append!(is_proj, true)
        end

        m, Ũ, Ṽ, b̃, jb̃ = u_norm ≤ r ? (1, U, V, b, jb) : (2, Up, Vp, bp, jbp)

        print("Semismooth Newton step.\n")
        Ah = interpolate_everywhere(Q.Φ₀ + Q.ϕ ⋅ Tᵢ - uB, Vu)
        A = findall(Ah.free_values .≤ 0) .+ (m*Vu.nfree + VT.nfree)
        inac = setdiff(1:(VT.nfree+(m+1)*Vu.nfree), A)

        op = FEOperator(b̃, jb̃, Ũ, Ṽ)
        zh = FEFunction(Ũ, [uᵢ.free_values[:]; zeros(VT.nfree+m*Vu.nfree)])

        res, J  = Gridap.Algebra.residual_and_jacobian(op, zh);
        dz = zeros(VT.nfree+(m+1)*Vu.nfree)
        dz[inac] = -J[inac,inac] \ res[inac]
        dzh = FEFunction(zh.fe_space, dz)
        δuN = dzh.single_fe_functions[1]
        uN = FEFunction(Vu, uᵢ.free_values[:] + δuN.free_values[:])

        puN, TN, SN, newton_its, pf_its, hik_its = inner_solve(Q, uN, Tᵢ, proj_rc, tol, bt, PF, 1e-2, newton_its, pf_its, hik_its)

        h1N = h1(Q, uN, SN)

        if globalization == true
            puB, TB, SB, newton_its, pf_its, hik_its = inner_solve(Q, uB, Tᵢ, proj_rc, tol, bt, PF, 1e-2, newton_its, pf_its, hik_its)

            h1B = h1(Q, uB, SB)
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