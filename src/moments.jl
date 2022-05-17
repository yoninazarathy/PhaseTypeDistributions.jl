"""
QQQQ
"""
function mean(d::PHDist)
    return -d.α * inv(d.T) * ones(model_size(d))
end

"""
QQQQ
"""
function var(d::PHDist)
    second_moment = 2d.α*inv(d.T)^2*ones(model_size(d))
    return second_moment - mean(d)^2
end

"""
Returns a vector of conditional scvs (conditional on abosrbing state)
"""
function scv(d::PHDist)
    return var(d)/mean(d)^2
end

"""
Returns a vector of absorption probabilities for the MAPH dist
"""
function absorption_probs(d::MAPHDist)
    d = non_degenerate_maph(d)
    D = Diagonal(d.T)
    return -d.α*inv(d.T)*d.T0
end

"""
Returns a vector of conditional means (conditional on absorbing state)
"""


function mean(d::MAPHDist)
    d = non_degenerate_maph(d)
    return (d.α * inv(d.T)^2* d.T0)./absorption_probs(d)
end

"""
Returns a vector of conditional variances (conditional on absorbing state)
"""
function var(d::MAPHDist)
    d = non_degenerate_maph(d)
    second_moment = -2*d.α*inv(d.T)^3*d.T0
    return second_moment./absorption_probs(d) .- (mean(d).^2)
end

"""
Returns a vector of conditional scvs (conditional on absorbing state)
"""
function scv(d::MAPHDist)
    return var(d) ./ (mean(d).^2) 
end