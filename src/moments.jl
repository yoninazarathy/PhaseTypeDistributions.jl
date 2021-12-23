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
Returns a vector of conditional scvs (conditional on abosrbing state)
"""
function scv(d::PHDist)
    return var(d)/mean(d)^2
end

"""
Returns a vector of absorption probabilities for the MAPH dist
"""
function absorption_probs(d::MAPHDist)
    D = Diagonal(d.T)
    PT  = I-inv(D)*d.T
    PT0 = -inv(D)*d.T0
    return d.α*inv(I-PT)*PT0 
end

"""
Returns a vector of conditional means (conditional on absorbing state)
"""
function mean(d::MAPHDist)
    return -d.α*inv(d.T)*d.T0
end

"""
Returns a vector of conditional variances (conditional on absorbing state)
"""
function var(d::MAPHDist)
    second_moment = 2d.α*inv(d.T)^2*d.T0
    return second_moment.-(mean(d).^2)
end

"""
Returns a vector of conditional scvs (conditional on absorbing state)
"""
function scv(d::MAPHDist)
    return var(d)./(mean(d).^2) #QQQQ
end