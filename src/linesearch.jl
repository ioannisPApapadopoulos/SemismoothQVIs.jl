function merit_fcn(Q::GeneralizedThermoformingQVI, uN, SN)
    h1(Q, uN, SN) #+ h1(Q, uN, puN)
end

function ls_fϕ(Q::GeneralizedThermoformingQVI, u, du, ls_α::T, Tᵢ, proj_rc, in_tol, hik_tol, bt, PF, FS, ρ_min, newton_its, pf_its, hik_its) where T
    Uu, Vu = Q.fe_space_u
    uN = FEFunction(Vu, u.free_values[:] + ls_α*du.free_values[:])
    puN, TN, SN, newton_its, pf_its, hik_its = inner_solve(Q, uN, Tᵢ, proj_rc, in_tol, hik_tol, bt, PF, FS, 1e-2, ρ_min, newton_its, pf_its, hik_its, show_trace=false)
    merit_fcn(Q, uN, SN)
end

function bt_linesearch(Q::GeneralizedThermoformingQVI, u, du, Tᵢ, proj_rc, in_tol, hik_tol, bt, PF, FS, ρ_min, newton_its, pf_its, hik_its)
    ls = BackTracking()
    Uu, Vu = Q.fe_space_u
    
    ls_ϕ = ls_α -> ls_fϕ(Q,u,du,ls_α, Tᵢ, proj_rc, in_tol, hik_tol, bt, PF, FS, ρ_min, newton_its, pf_its, hik_its)

    uN = FEFunction(Vu, u.free_values[:])
    puN, TN, SN, newton_its, pf_its, hik_its = inner_solve(Q, uN, Tᵢ, proj_rc, in_tol, hik_tol, bt, PF, FS, 1e-2, ρ_min, newton_its, pf_its, hik_its, show_trace=false)

    ls_ϕ0 = merit_fcn(Q, uN, SN)
    ls_dϕ0 = sum(∫((uN-SN) ⋅ (du) + ∇(uN-SN) ⋅ ∇(du))*Q.dΩ)/h1(Q, uN, SN)  #+ sum(∫((uN-puN) ⋅ (du) + ∇(uN-puN) ⋅ ∇(du))*Q.dΩ)/h1(Q, uN, puN) 
    # print("ls_ϕ0: $ls_ϕ0, ls_dϕ0: $ls_dϕ0 \n")

    ls_α, _ =  ls(ls_ϕ, 1.0, ls_ϕ0, ls_dϕ0)
    ls_α
end