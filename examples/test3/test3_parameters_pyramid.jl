k = 1

Φ₀ = interpolate_everywhere(x->1.0 - 2*max(abs(x[1]-1/2),  abs(x[2]-1/2)), VT)
Ψ₀ = Φ₀
ϕ = interpolate_everywhere(x->sin(π*x[1])*sin(π*x[2]), UT)
ψ = ϕ

f = interpolate_everywhere(x->25, Uu)

# g and its derivative
function g(s)
    if s ≤ 0
        return 1/5
    elseif 0 < s < 1
        return (1-s)/5
    else
        return 0.0
    end
end

function dg(s)
    if s ≤ 0
        return 0.0
    elseif 0 < s < 1
        return -1/5
    else
        return 0.0
    end
end
