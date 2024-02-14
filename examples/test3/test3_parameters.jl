k = 1

function w(r)
	if 0.1 ≤ r ≤ 0.3
		return 5.0 * r - 0.5
	elseif 0.3 < r < 0.7
		return 1.0
	elseif 0.7 <= r <= 0.9
		return 4.5 - 5.0 * r
	else
		return 0.0
	end
end

## initial mould
function Φ0(x)
	return w(x[1]) * w(x[2])
end

Φ₀ = interpolate_everywhere(x->Φ0(x), VT)
Ψ₀ = interpolate_everywhere(x->1.0, VT)

bump_1d(x) = (0.0 <= x[1] <= 1.0) ? exp(-0.25 / (x[1] - x[1]^2)) : 0.0
bump(x) = begin
	r = sqrt((x[1] - 0.5)^2 + (x[2] - 0.5)^2)
	return bump_1d(0.5 + r)
end
ϕ = interpolate_everywhere(x->0.1*bump(x), UT)
ψ = ϕ

f = interpolate_everywhere(x->40, Uu)
# function g̃(r,κ,s)
# 	if r <= 0.0
# 		return κ
# 	elseif 0 < r <= s / 4
# 		return κ - 8.0 * κ * r^2 / (3.0 * s^2)
# 	elseif s/4 < r ≤ 0.75 * s
# 		return 7.0 * κ / 6.0 - 4.0 / 3.0 * κ * r / s
# 	elseif 0.75*s < r ≤ s
# 		return 8.0 / 3.0 * (s - r)^2 / s^2
# 	else
# 		return 0.0
# 	end
# end

# function dg̃(r,κ,s)
# 	if r <= 0.0
# 		return 0.0
# 	elseif 0 < r <= s / 4
# 		return -16.0 * κ * r / (3.0 * s^2)
# 	elseif s/4 < r ≤ 0.75 * s
# 		return - 4.0 / 3.0 * κ / s
# 	elseif 0.75*s < r ≤ s
# 		return - 16.0 / 3.0 * (s - r) / s^2
# 	else
# 		return 0.0
# 	end
# end

# g(s) = 0.1*g̃(s,50.0,1.0)
# dg(s) = 0.1*dg̃(s, 50.0, 1.0)

# g and its derivative
function g(x)
    if x ≤ 0
        return 0.1
    elseif 0 < x < 1
        return (1-x)/10
    else
        return 0.0
    end
end

function dg(x)
    if x ≤ 0
        return 0.0
    elseif 0 < x < 1
        return -1/10
    else
        return 0.0
    end
end
