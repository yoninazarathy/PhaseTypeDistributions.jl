"""
A mutable struct representing a probability distribution with parameters α and a transition matrix T.
    α::Adjoint{Float64, Vector{Float64}}`: The parameters α of the distribution.
    T::Matrix{Float64}`: The transition matrix T of the distribution.
"""

mutable struct PHDist
    α::Adjoint{Float64, Vector{Float64}}
    T::Matrix{Float64}
    PHDist(α::Adjoint{Float64, Vector{Float64}}, T::Matrix{Float64}) = new(α, T)
end

function build_ph(α::Adjoint{Float64, Vector{Float64}}, T::Matrix{Float64})
    if size(T, 1) != size(T, 2)
        error("T must be a square matrix")
    end
    
    if length(α) != size(T, 1)
        error("The length of α must be equal to the number of rows in T")
    end

    return PHDist(α, T)
end


"""
MAPH distribution
# Fields
- `m::Integer`: Number of phases in the MAPH distribution.
- `n::Integer`: Number of absorbing states ain the MAPH.
- `α::Adjoint{Float64, Vector{Float64}}`: Initial probability distribution across the phases. Elements should sum to 1.
- `T::Matrix{Float64}`: `m x m` transition rate matrix for the MAPH states.
- `T0::Matrix{Float64}`: `m x m` abosrbing rate matrix for the MAPH states.
-  q::Vector{Float64}: `m` vector of transition rates from the MAPH states to the absorbing states.
- `R::Matrix{Float64}`: `m x n` absorbing probabilities matrix.
- `P::Vector{Matrix{Float64}}`: Vector of `n` `m x m` matrices, representing the transition probabilities between the MAPH states given absorbing states.
"""
mutable struct MAPHDist
    α::Adjoint{Float64, Vector{Float64}} 
    T::Matrix{Float64} 
    T0::Matrix{Float64}
    q:: Vector{Float64}
    #second parameterization
    R::Matrix{Float64} #dimension at mxn
    P::Vector{Matrix{Float64}} #a vector of length n, with elements at dimension mxm
    MAPHDist(α::Adjoint{Float64, Vector{Float64}}, T::Matrix{Float64}, T0::Matrix{Float64}, q::Vector{Float64}, R::Matrix{Float64}, P::Vector{Matrix{Float64}}) = new(α, T, T0, q, R, P)
end

"""
Takes T, T0 of first parameterization and returns matrices of second parameterization.
"""
function R_P_from_T_T0(T, T0)
    m, n  = size(T0)
    q = -diag(T)
    R = -inv(T)*T0
    P = [[i != j ? -T[i,j] * R[j,k] / T[i,i] / R[i,k] : 0 for i in 1:m, j in 1:m] for k in 1:n]

    return q, R, P
end

"""
Takes R, P, q of second parameterization and returns matrices of first parameterization.
"""
function T_T0_from_R_P_q(q::Vector{Float64}, R::Matrix{Float64}, P::Vector{Matrix{Float64}})
    m = length(q)
    k = 1 #QQQQ this is some fixed k (abosrbing state - as conversion will work the same for all k)
    T = [i==j ? -q[i] : q[i] * P[k][i,j] * R[i,k] / R[j,k] for i in 1:m, j in 1:m]
    T0 = -T * R
    return T, T0
end

function build_maph(α::Adjoint{Float64, Vector{Float64}}, T::Matrix{Float64}, T0::Matrix{Float64})
    if size(T, 1) != size(T, 2)
        error("T must be a square matrix")
    end

    if size(q,1) != size(T,1)
        error("q must be the same length as the number of rows in T")
    end
    
    if size(T0, 1) != size(T, 1)
        error("T0 must be a mxn matrix")
    end
    
    if length(α) != size(T, 1)
        error("The length of α must be equal to the number of rows in T")
    end
    
    return MAPHDist(α, T, T0, R_P_from_T_T0(T, T0)...)
end

build_maph(α::Adjoint{Float64, Vector{Float64}}, R::Matrix{Float64}, P::Matrix{Float64})= MAPHDist(α, T_T0_from_R_P_q(q, R, P)..., q, R, P)



"""
QQQQ -doc to do
"""
function update_params_1to2!(dist::MAPHDist)
    dist.R, dist.P = R_P_from_T_T0(dist.T, dist.T0)
end

"""
QQQQ -doc to do
"""
function update_params_2to1!(dist::MAPHDist)
    dist.T, dist.T0 = T_T0_from_R_P_q(diag(dist.T), dist.R, dist.P)
end


"""
QQQQ
"""
mutable struct MAPHSufficientStats
    B::Vector{Float64} #initial starts
    Z::Vector{Float64} #time spent
    M::Matrix{Float64} #transitions between transient states
    N::Matrix{Float64} #transitions between transient to abosrbing states

    MAPHSufficientStats(B::Vector{Float64}, Z::Vector{Float64}, M::Matrix{Float64},  N::Matrix{Float64}) = new(B,Z,M,N)
    function MAPHSufficientStats(maph::MAPHDist) 
        m, n = model_size(maph)
        new(zeros(m), zeros(m), zeros(m,m),zeros(m,n))
    end
end

"""
QQQQ
"""
SingleObs = NamedTuple{(:y, :a), Tuple{Float64, Int64}}

"""
QQQQ
"""
MAPHObsData = Vector{SingleObs}

+(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B+ss2.B, ss1.Z+ss2.Z, ss1.M+ ss2.M, ss1.N+ss2.N)
/(ss::MAPHSufficientStats,n::Real) = MAPHSufficientStats(ss.B/n, ss.Z/n, ss.M/n, ss.N/n)
/(ss1::MAPHSufficientStats,ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B ./ ss2.B, ss1.Z ./  ss2.Z, ss1.M ./  ss2.M,  ss1.N ./  ss2.N)
-(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B-ss2.B, ss1.Z-ss2.Z, ss1.M - ss2.M, ss1.N - ss2.N)


""""
QQQQ - put doc string
"""
model_size(ph::PHDist) = length(ph.α)


"""
QQQQ - put doc string
"""
model_size(maph::MAPHDist) = (p = size(maph.T,1), q = size(maph.T0,2)) #transient is p and abosrbing is q

"""
QQQQ
"""
function non_degenerate_maph(maph::MAPHDist)::MAPHDist
    valid_states = abs.(diag(maph.T)) .> sqrt(eps())
    return MAPHDist(maph.α[valid_states]', maph.T[valid_states,valid_states], maph.T0[valid_states,:])
end