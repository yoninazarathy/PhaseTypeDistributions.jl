"""
QQQQ
"""
function mean(d::PHDist)
    return -d.α*inv(d.T)*ones(model_size(d))
end

"""
QQQQ
"""
function var(d::PHDist)
    second_moment = 2d.α*inv(d.T)^2*ones(model_size(d))
    return second_moment - mean(d)^2
end

"""
QQQQ
"""
function scv(d::PHDist)
    return var(d)/mean(d)^2
end

"""
Returns a vector of absorption probabilities for the MAPH dist
"""
function absorption_probs(d::MAPHDist)
    return 0 #QQQQ
end

"""
Returns a vector of conditional means (conditional on absorbing state)
"""
function mean(d::MAPHDist)
    return 0 #QQQQ
end

"""
Returns a vector of conditional variances (conditional on absorbing state)
"""
function var(d::MAPHDist)
    return 0 #QQQQ
end

"""
Returns a vector of conditional scvs (conditional on absorbing state)
"""
function scv(d::MAPHDist)
    return 0 #QQQQ
end