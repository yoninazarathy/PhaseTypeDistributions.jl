"""
$(TYPEDEF)
The mutable struct which representing the MAPH distribution in different parameterizations
"""
mutable struct MAPHDist 
    """m × 1 Initial probability distribution across the phases. Elements should sum to 1"""
    α::Matrix{<:Real}

    """m x m transition rate matrix for the MAPH states"""    
    T::Matrix{<:Real}

    """m x n abosrbing rate matrix for the MAPH states."""
    D::Matrix{<:Real}

    """m × 1 vector of transition rates from the MAPH states to the absorbing states."""
    q::Vector{<:Real}

    """m x n absorbing probabilities matrix."""
    R::Matrix{<:Real}

    """Vector of n copies of m x m matrices, representing the transition probabilities between the MAPH states given absorbing states."""
    P::Vector{Matrix{<:Real}}
end

function is_degenerate_maph(T::Matrix{<:Real})
    return !all(abs.(diag(T)) .> sqrt(eps()))
end

function MAPH_constructor(α::Matrix{<:Real}, T::Matrix{<:Real}, D::Matrix{<:Real}, q::Vector{<:Real}, R::Matrix{<:Real}, P::Vector{ Matrix{<:Real}})
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

function R_P_from_T_D(T::Matrix{<:Real}, D::Matrix{<:Real})
    m, n  = size(D)
    q = -diag(T)
    R = -inv(T)*D
    # P = map(k-> [i != j ? -T[i,j] * R[j,k] / T[i,i] / R[i,k] : 0 for i in 1:m, j in 1:m], 1:n)
    k = 1
    P1 = [i != j ? -T[i,j] * R[j,k] / T[i,i] / R[i,k] : 0.0 for i in 1:m, j in 1:m]
    P = P_update_by_R(R, P1)
    replace!.(P, NaN => 0)
    return q, R, P
end

function P_update_by_R(R::Matrix{<:Real}, P1::Matrix{<:Real})
    m, n = size(R)
    new_P = map(k -> [ P1[i,j] * R[i,1]*R[j,k]/(R[j,1]*R[i,k]) for i ∈ 1:m, j ∈ 1:m], 2:n)
    pushfirst!(new_P, P1)
    return new_P
end


function is_valid_R_P(R::Matrix{<:Real}, P::Vector{ Matrix{<:Real}})
    new_P = P_update_by_R(R, P[1])
    @show new_P
    return all(new_P .≈  P)
end

function maph_random_parameters(m::T, n::T, parameterization = :RP) where {T <: Int}
    if parameterization == :TD
        error("not implmented")
    elseif parameterization == :RP
        v = rand(m)
        α = v / sum(v)
        α = reshape(α, 1, m)
        q = rand(Exponential(1), m)
        R = rand(m, n)
        R = R ./ sum(R, dims = 2)
        P1 = rand(m, m)
        P1 = P1 - Diagonal(P1)
        P1 = P1 ./ sum(P1, dims = 2) .* rand(m)
        P_collection = P_update_by_R(R, P1)

        return MAPH_constructor(α, q, R, P_collection)
    else
        error("only RP and TD")
    end

end


function T_D_from_R_P_q(q::Vector{<:Real}, R::Matrix{<:Real}, P::Vector{<:Matrix{<:Real}})
    #add equation 7 check the conditon 
    #add a helper function to check the conditions and fix P 
    m = length(q)
    valid_P = findall(map(k -> !any(isnan.(P[k])), 1:length(P)))
    @show valid_P
    k = 1
    T = [i==j ? -q[i] : q[i] * P[k][i,j] * R[i,k] / R[j,k] for i in 1:m, j in 1:m]
    D = -T * R
    return T, D
end

# function is_valid_α_R_P_q(maph::MAPHDist)
#     P, R  = maph.P, maph.R
#     m, n = length(maph.α), length(P)
#     R_factors = map(k -> [ P[1][i,j] * R[i,1]*R[j,k]/(R[j,1]*R[i,k]) for i ∈ 1:m, j ∈ 1:m], 2:n)
    
#     new_P = reduce(vcat, [P[1], [P[k] * R_factors[k] for k = 2:n]])

#     isvalid = all(map(k -> P[k] ≈ R_factors[k-1], 2:n))

    

#     if !isvalid
#         diff = map(k -> P[k] - R_factors[k-1], 2:n)
#         display(diff[2])
#     end
#     return isvalid

# end

# function is_valid_α_T_D(maph::MAPHDistributions)


function MAPH_constructor(α::Matrix{<:Real}, T:: Matrix{<:Real}, D::Matrix{<:Real})
    if is_degenerate_maph(T)
        valid_states = abs.(diag(T)) .> sqrt(eps())
        new_α = α[:, valid_states]
        new_T = T[valid_states, valid_states]
        new_D = D[valid_states, valid_states]
        return MAPHDist(new_α, new_T, new_D, R_P_from_T_D(new_T, new_D)...)
    end    
    return MAPHDist(α, T, D, R_P_from_T_D(T, D)...)
end

MAPH_constructor(α::Matrix{<:Real}, q::Vector{<:Real}, R::Matrix{<:Real}, P::Vector{<:Matrix{<:Real}})= MAPHDist(α, T_D_from_R_P_q(q, R, P)..., q, R, P)

function update!(maph::MAPHDist, q::Vector{<:Real}, R::Matrix{<:Real}, P::Vector{<:Matrix{<:Real}})
    T, D = T_D_from_R_P_q(q, R, P)
    maph.T = T
    maph.D = D
    maph.q = q
    maph.R = R
    maph.P = P
end

function update!(maph::MAPHDist, T::Matrix{<:Real}, D::Matrix{<:Real})
    q, R, P = R_P_from_T_D(T, D)
    maph.T = T
    maph.D = D
    maph.q = q
    maph.R = R
    maph.P = P
end

