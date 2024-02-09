"""
$(TYPEDFIELDS)
The mutable struct which representing the MAPH distribution in different parameterizations
$(TYPEDEF)
"""
mutable struct MAPHDist
    #first parameterization
    "m × 1 Initial probability distribution across the phases. Elements should sum to 1"
    α::Union{Matrix{Float64}, Vector{Float64}}

    "m x m transition rate matrix for the MAPH states"    
    T::Matrix{Float64}

    "m x n abosrbing rate matrix for the MAPH states."
    D::Matrix{Float64}

    #second parameterization
    "m × 1 vector of transition rates from the MAPH states to the absorbing states."
    q::Vector{Float64}

    "m x n absorbing probabilities matrix."
    R::Matrix{Float64}

    " Vector of n copies of m x m matrices, representing the transition probabilities between the MAPH states given absorbing states."
    P::Vector{Matrix{Float64}}

    MAPHDist(α::Union{M, Vector{Float64}}, T::M, D::M, q::Vector{Float64}, R::M, P::Vector{M})  where {M <: Matrix{Float64}}= new(α, T, D, q, R, P)

end

function is_degenerate_maph(T::Matrix{Float64})
    return !all(abs.(diag(T)) .> sqrt(eps()))
end

function MAPH_constructor(α::Union{M, Vector{Float64}}, T::M, D::M, q::Vector{Float64}, R::M, P::Vector{M}) where {M <: Matrix{Float64}}
    @assert size(α, 2) == size(T,1) "must have the same num of phases"
    @assert size(T, 1) == size(T, 2) "T must be a square matrix"
    @assert size(q, 1) == size(T,1) "q must be the same length as the number of rows in T"
    @assert size(D, 1) == size(T, 1) "T0 must be a mxn matrix"
    @assert length(α) == size(T, 1) "The length of α must be equal to the number of rows in T"
    @assert size(D) == size(R) "prob and rate matrix must have the same dimension"
    @assert sum(α) == 1 "initial prob must sum to 1"
    #adding more model constraints here.. e.g. the diagonal must be dominant
    return MAPHDist(α, T, D, q, R, P)
end

"""
Takes T, D of first parameterization and returns matrices of second parameterization.
"""
function R_P_from_T_D(T::M, D::M) where {M <: Matrix{Float64}}
    m, n  = size(D)
    q = -diag(T)
    R = -inv(T)*D
    # P = map(k-> [i != j ? -T[i,j] * R[j,k] / T[i,i] / R[i,k] : 0 for i in 1:m, j in 1:m], 1:n)
    P = [[i != j ? -T[i,j] * R[j,k] / T[i,i] / R[i,k] : 0.0 for i in 1:m, j in 1:m] for k in 1:n]
    return q, R, P
end

"""
Takes R, P, q of second parameterization and returns matrices of first parameterization.
"""
function T_D_from_R_P_q(q::Vector{Float64}, R::M, P::Vector{M}) where {M <: Matrix{Float64}}
    m = length(q)
    k = 1 #QQQQ this is some fixed k (abosrbing state - as conversion will work the same for all k)
    T = [i==j ? -q[i] : q[i] * P[k][i,j] * R[i,k] / R[j,k] for i in 1:m, j in 1:m]
    D = -T * R
    return T, D
end

function MAPH_constructor(α::Union{M, Vector{Float64}}, T::M, D::M)  where {M <: Matrix{Float64}}
    if is_degenerate_maph(T)
        valid_states = abs.(diag(T)) .> sqrt(eps())
        new_α = α[:, valid_states]
        new_T = T[valid_states, valid_states]
        new_D = D[valid_states, valid_states]
        return MAPHDist(new_α, new_T, new_D, R_P_from_T_D(new_T, new_D)...)
    end    
    return MAPHDist(α, T, D, R_P_from_T_D(T, D)...)
end

MAPH_constructor(α::Union{M, Vector{Float64}}, q::Vector{Float64}, R::M, P::Vector{M}) where {M <: Matrix{Float64}} = MAPHDist(α, T_D_from_R_P_q(q, R, P)..., q, R, P)

function update!(maph::MAPHDist, q::Vector{Float64}, R::M, P::Vector{M}) where {M <: Matrix{Float64}}
    T, D = T_D_from_R_P_q(q, R, P)
    maph.T = T
    maph.D = D
    maph.q = q
    maph.R = R
    maph.P = P
end



function update!(maph::MAPHDist, T::Matrix{Float64}, D::Matrix{Float64})
    q, R, P = R_P_from_T_D(T, D)
    maph.T = T
    maph.D = D
    maph.q = q
    maph.R = R
    maph.P = P
end

