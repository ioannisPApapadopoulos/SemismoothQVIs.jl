import LinearAlgebra: Diagonal
# include("deflation.jl")
# import LinearAlgebra: cond
## Newton

function zero_dirichlet_bcs(r, J)
    r[1] = r[end] = 0
    J = Matrix(J)
    J[:,1] .= 0; J[1,:] .= 0; J[1,1] = 1;
    J[:,end] .= 0; J[end, :] .= 0; J[end, end] = 1;
    return r, J
end

# function defl_τ(u, du, us, M)
#     deflation_step_adjustment(u, du, us, M)
# end

function newton(op, uh; tol=1e-12, max_iter=1000, damping=1, us=[], M=[], info=false, show_trace=true)
    its, res = 0, 1.

    us_cfs = isempty(us) ? [] : [u.free_values[:] for u in us]
    
    res, J  = Gridap.Algebra.residual_and_jacobian(op, uh)
    show_trace && print("Newton: Iteration $its, residual norm = $(norm(res))\n")
    ucfs = uh.free_values[:]
    vh = FEFunction(uh.fe_space,ucfs)

    while norm(res) > tol && its ≤ max_iter
        its += 1
        # κ = cond(Array(J), 2)
        # print("κ(J) = $κ. \n")
        du = - J \ res
        ucfs = vh.free_values

        τ = isempty(us) ? 1 : defl_τ(ucfs, du, us_cfs, M)
        # show_trace && print("tau = $τ\n")
        cfs = ucfs + damping*τ*du

        vh = FEFunction(uh.fe_space,cfs)

        res, J  = Gridap.Algebra.residual_and_jacobian(op, vh)
        show_trace && print("Newton: Iteration $its, residual norm = $(norm(res))\n")
    end

    if info == true
        return vh, its
    else
        return vh
    end
end


## Benson-Munson


function reduced_residual(r::AbstractVector{T}, x::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) where T
    rr = r[:];
    rr[x .<= lb] = min.(rr[x .<= lb], zero(T))
    rr[x .>= ub] = max.(rr[x .>= ub], zero(T))
    return rr
end

function project!(x, lb, ub)
    b = x;
    b[x .< lb] = lb[x .< lb]
    b[x .> ub] = ub[x .> ub]
    return b
end

function project(x, lb, ub)
    b = x[:];
    b[x .< lb] = lb[x .< lb]
    b[x .> ub] = ub[x .> ub]
    return b
end

function bm(op, uh, lb::AbstractVector{T}, ub::AbstractVector{T}; us::AbstractVector=[], M=[], tol::T=1e-9, max_iter::Int=1000, damping=1, info=false, show_trace=true) where T

    x = uh.free_values[:]
    known_roots = [u.free_values[:] for u in us]

    index = Vector(1:lastindex(x))
    iter = 0

    project!(x, lb, ub)

    r, J  = Gridap.Algebra.residual_and_jacobian(op, uh)

    active_lb = findall(x .≈ lb) ∩ findall(r .> 0)
    active_ub = findall(x .≈ ub) ∩ findall(r .< 0)
    active = vcat(active_lb, active_ub)

    tmp_index = index[:]
    tmp_index[active] .= 0
    inactive  = findall(x->x!=0, tmp_index)

    norm_residual_Ω = norm(reduced_residual(r, x, lb, ub))
    show_trace && print("BM: Iteration 0, residual norm = $norm_residual_Ω\n")

    n = length(x)

    while norm_residual_Ω > tol && iter < max_iter
        update = zeros(T, n)
        update[inactive] = -J[inactive, inactive] \ r[inactive]
        # project!(update,lb-x,ub-x)

        # if norm(update) < 1e-10
        #     update = -r
        # end

        # if ~isempty(known_roots)
        #     τ = deflation_step_adjustment(x, update, known_roots, M)
        #     update = τ * update
        # end

        # update = nls.LineSearch.adjust(x, update, nls.damping)
        # for α in 1.0 ./ (1:100)
        #     y = x + α * damping * update
        #     project!(y,lb,ub)
        #     vh = FEFunction(uh.fe_space,x)
        #     rr = Gridap.Algebra.residual(op, vh)
        #     normres = norm(reduced_residual(rr, y, lb, ub))
        #     if normres ≤ norm_residual_Ω || α ≈ 0.01
        #         x = y
        #         r = rr
        #         norm_residual_Ω = normres
        #         uh = vh
        #         break
        #     end
        # end

        x = x + damping * update
        project!(x,lb,ub)
        uh = FEFunction(uh.fe_space,x)

        # xx = range(0,1,100)
        # p = Plots.plot(xx, uh.(Point.(xx)))
        # display(p)

        r, J  = Gridap.Algebra.residual_and_jacobian(op, uh)
        norm_residual_Ω = norm(reduced_residual(r, x, lb, ub))

        active_lb = index[x .<= lb]
        active_lb = intersect(active_lb, index[r .> zero(T)])
        active_ub = index[x .>=ub]
        active_ub = intersect(active_ub, index[r .< zero(T)])
        active = vcat(active_lb, active_ub)
        index2 = index[:]
        index2[active] .= 0
        inactive  = findall(x->x!=0, index2)

        iter += 1
        show_trace && print("BM: Iteration $iter, residual norm = $norm_residual_Ω\n")
    end
    
    if iter == max_iter
        show_trace && print("BM: Iteration max reached")
        # disp(normResidualOmega)
    end
    if info == true
        return uh, iter
    else
        return uh
    end
end

## HIK

function hik(op, uh, lb::AbstractArray{T}, ub; us::AbstractVector=[], M=[], tol::T=1e-9, max_iter::Int=1000, damping=1, info=false, show_trace=true) where T
    
    x = uh.free_values[:]

    known_roots = [u.free_values[:] for u in us]

    index = Vector(1:lastindex(x))
    iter = 0
    inactive = index

    active_lb = []
    active_ub = []
    active    = []
    
    project!(x, lb, ub)
    vh = FEFunction(uh.fe_space,x)
    n = length(x)
    
    r, J  = Gridap.Algebra.residual_and_jacobian(op, vh)

    active_lb = findall(x .≈ lb) ∩ findall(r .> 0)
    active_ub = findall(x .≈ ub) ∩ findall(r .< 0)
    active = vcat(active_lb, active_ub)

    # print(active)

    dual = zeros(n)
    dual[active_lb] = r[active_lb]
    dual[active_ub] = -r[active_ub]
    tmp_index = index[:]
    tmp_index[active] .= 0
    inactive  = findall(x->x!=0, tmp_index)
    
    norm_residual_Ω = norm(reduced_residual(r, x, lb, ub))
    # %normResidualOmega = norm(dual - max(zeros(n,1), dual+x-ub) - max(zeros(n,1), dual-x-lb));
    # %normResidualOmega = normResidualOmega + norm(evaluatedResidual+dual);
    show_trace && print("HIK: Iteration 0, residual norm = $norm_residual_Ω\n")
    
    
    while (norm_residual_Ω) > tol && (iter < max_iter)
        
        update = zeros(T, n);
        update[active_lb] = lb[active_lb] - x[active_lb]
        update[active_ub] = ub[active_ub] - x[active_ub]

        if ~isempty(active)
            cr = r[inactive] + J[inactive, active]*update[active]
        else
            cr = r[inactive]
        end
        update[inactive] = -J[inactive,inactive] \ cr

        # if ~isempty(known_roots)
        #     τ = deflation_step_adjustment(x, update, known_roots, M)
        #     update = τ * update
        # end

        x = x + damping*update;
        vh = FEFunction(uh.fe_space,x)

        # xx = range(0,1,100)
        # p = Plots.plot(xx, vh.(Point.(xx)))
        # display(p)
        # sleep(1)
    
        # which way should the sign be?
        dual[inactive] .= zero(T);
        dual[active_lb] = r[active_lb]
        dual[active_ub] = -r[active_ub]
        
        active_lb = findall((dual .- x .+ lb).>0)
        active_ub = findall((dual .-ub .+x) .>0)
        active = vcat(active_lb, active_ub)

        # print(active)

        tmp_index = index[:]
        tmp_index[active] .= 0
        inactive  = findall(x->x!=0, tmp_index)
        
        # norm_residual_Ω = norm(reduced_residual(r, x, lb, ub))
        project!(x,lb,ub)
        vh = FEFunction(uh.fe_space,x)
        r, J  = Gridap.Algebra.residual_and_jacobian(op, vh)
        # norm_residual_Ω  = norm(dual - max.(zeros(n,1), dual+x-ub) - max.(zeros(n,1), dual-x+lb))
        # norm_residual_Ω  = norm_residual_Ω  + norm(r+dual)
        norm_residual_Ω = norm(reduced_residual(r, x, lb, ub))


        iter += 1
        show_trace && print("HIK: Iteration $iter, residual norm = $norm_residual_Ω\n")
    end
    
    if iter == max_iter
        show_trace && print("HIK: Iteration max reached")
        @warn("HIK: Iteration max reached.")
        # disp(normResidualOmega)
    end
    if info == true
        return vh, iter
    else
        return vh
    end
end


## SSLS

function Φ(a::AbstractVector{T}, b::AbstractVector{T}) where T
    @assert length(a) == length(b)
    a + b - sqrt.(a.^2 + b.^2)     
end

function dΦ(a::AbstractVector{T}, b::AbstractVector{T}) where T
    @assert length(a) == length(b)
    if any(abs.(a) .> 1e-6) || any(abs.(b) .> 1e-6)
        return one(T) .- a./sqrt.(a.^2 + b.^2)
    else
        return ones(T, length(a)) ./ 2
    end
end

function FB(x, r, lb, ub, wherenoconstraint, wherelbconstraint, whereubconstraint, whereequalconstraint, wherebothconstraint)
    T = eltype(x)
    out = zeros(T,length(x))
    #  FIXME add a check for all indices here
    out[wherenoconstraint] = r[wherenoconstraint]
    
    idx = whereubconstraint
    out[idx] = Φ(ub[idx] - x[idx], -r[idx])
    
    idx = wherelbconstraint
    out[idx] = Φ(x[idx]-lb[idx],r[idx])
    
    idx = wherebothconstraint
    out[idx] = Φ(x[idx]-lb[idx],Φ(ub[idx]-x[idx],-r[idx]))
    
    idx = whereequalconstraint
    out[idx] = lb[idx] - x[idx]
    return out
end

function computeScaleAndShift(x, r, lb, ub, wherenoconstraint, wherelbconstraint, whereubconstraint, whereequalconstraint, wherebothconstraint)
    T = eltype(x)
    n = length(x)
    dshift = ones(T, n)
    dscale = ones(T, n)
    
    dshift[wherenoconstraint] .= zero(T)
    dscale[wherenoconstraint] .= one(T)
    
    idx = whereubconstraint
    dshift[idx] = dΦ(ub[idx] - x[idx], -r[idx])
    dscale[idx] = dΦ(-r[idx],ub[idx]-x[idx])
    
    idx = wherelbconstraint
    dshift[idx] = dΦ(x[idx] - lb[idx], r[idx])
    dscale[idx] = dΦ(r[idx], x[idx] - lb[idx])
   
    
    idx = wherebothconstraint;
    dshift1 = dΦ(x[idx] - lb[idx], -Φ(ub[idx] - x[idx], -r[idx]))
    dscale1 = dΦ(-Φ(ub[idx]-x[idx],-r[idx]), x[idx] - lb[idx])
    dshift2 = dΦ(ub[idx]-x[idx],-r[idx])
    dscale2 = dΦ(-r[idx],ub[idx] - x[idx])
    dshift[idx] = dshift1 + dscale1.*dshift2
    dscale[idx] = dscale1.*dscale2
    
    idx = whereequalconstraint
    dshift[idx] .= one(T)
    dscale[idx] .= zero(T)
    return (dshift, dscale)
end

function ssls(op, uh, lb::AbstractVector{T}, ub::AbstractVector{T};
     us::AbstractVector=[], tol::T=1e-9, max_iter::Int=1000, damping=1) where T

    iter = 0
    x = uh.free_values[:]
    known_roots = [u.free_values for u in us]

    project!(x,lb,ub)
    
    bound_tol = 1e10
    wherenoconstraint = intersect(findall(x->x<-bound_tol, lb), findall(x->x>bound_tol, ub))
    wherelbconstraint = intersect(findall(x->x>=-bound_tol, lb), findall(x->x>bound_tol, ub))
    whereubconstraint = intersect(findall(x->x<-bound_tol,lb), findall(x->x<=bound_tol, ub))
    whereequalconstraint = findall(x->x==ub, lb)
    wherebothconstraint = intersect(findall(x->x>=-bound_tol,lb), findall(x->x<=bound_tol,ub))
    
    r, J  = Gridap.Algebra.residual_and_jacobian(op, uh)
             
    fb = FB(x, r, lb, ub, wherenoconstraint,wherelbconstraint,whereubconstraint,whereequalconstraint,wherebothconstraint)
    normFB = norm(fb)
    print("Iteration 0, residual norm = $normFB\n")
    
    while normFB > tol && iter < max_iter
        
        (dshift, dscale) = computeScaleAndShift(x, r, lb, ub, wherenoconstraint,wherelbconstraint,whereubconstraint,whereequalconstraint,wherebothconstraint)
        shiftedJacobian = Diagonal(dshift) + Diagonal(dscale) * J

        update = -shiftedJacobian \ fb
        
        # if ~isempty(known_roots)
        #     τ = deflation_step_adjustment(x, update, known_roots, M, p, α)
        #     update =  τ * update
        # end
        
        # update = nls.LineSearch.adjust(x, update, nls.damping)
        x = x + damping*update
        uh = FEFunction(uh.fe_space,x)
        
        r, J  = Gridap.Algebra.residual_and_jacobian(op, uh)  
      
        fb = FB(x, r, lb, ub, wherenoconstraint,wherelbconstraint,whereubconstraint,whereequalconstraint,wherebothconstraint)
        normFB = norm(fb)               
        iter += 1
        print("Iteration $iter, residual norm = $normFB\n")
    end
    
    if iter == max_iter
        print("Iteration max reached")
    end
    return uh
end