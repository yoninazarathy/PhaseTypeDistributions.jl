"""
$(TYPEDEF)

Mutable struct defines a phasetype distribution
"""

mutable struct PHDist

    "m × 1 Initial probability distribution across the phases. Elements should sum to 1"
    α::Matrix{<:Real}

    "m x m transition rate matrix for the MAPH states"    
    T::Matrix{<:Real}
    
end

function ph_constructor(α::Matrix{<:Real}, T::Matrix{<:Real})
    @assert size(T,1) == size(T, 2) "T must be a square matrix"
    @assert size(α, 2) == size(T, 1) "The length of α must be equal to the number of rows in T"
    return PHDist(α, T)
end


function hyper_exp_dist(mean_desired::Real, scv_desired::Real)
    @assert scv_desired > 1.0 "SCV must be greater than 1"
    μ1 = 1/(scv_desired+1) #mean parameter 
    p = (scv_desired-1)/(scv_desired+1+2/(μ1^2)-4/μ1)
    μ2 = (1-p)/(1-p/μ1) #mean parameter
    α = zeros(1, 2)
    α[1, 1] = p
    α[1, 2] = 1-p
    T = zeros(2,2)
    T[1,1] = -μ1
    T[2,2] = -μ2
    return ph_constructor(α, (1/mean_desired)*T)
end

function exp_dist(mean_desired::Real)
    T = zeros(1,1)
    α = ones(1)
    T[1,1] = -1/mean_desired
    return ph_constructor(α,T)
end

function hypo_exp_dist(mean::Real,scv::Real)
    @assert scv < 1.0 "SCV must be less than 1"
    n = Int(ceil(1/scv))
    ν1 = n/(1+sqrt((n-1)*(n*scv-1)))
    ν2 = ν1*(n-1)/(ν1-1)
    α = zeros(1, n)
    α[1, 1] = 1
    T = zeros(n,n)
    T[1,1] = -ν1
    T[1,2] = ν1
    for i = 2:(n-1)
        T[i,i] = -ν2
        T[i,i+1] = ν2
    end

    T[n,n] = -ν2

    return ph_constructor(α, (1/mean) * T) 
end

function get_absorbing_vector(ph::PHDist)
    return -1.0 * ph.T * ones(size(ph.T, 1))
end
# function mixed_expo_dist(p1::Int64,ω::Real)
#     α = ones(1, p1) / p 
#     α = (ones(p1)/p)'
#     return PHDist(α, Matrix(-ω*I(p1)))
# end












