"""
QQQQ
"""
mutable struct PHDist
    α::Adjoint{Float64, Vector{Float64}}
    T::Matrix{Float64}
end


"""
Takes T, T0 of first parameterization and returns matrices of second parameterization.
"""
function R_P_from_T_T0(T, T0)
    m, n  = size(T0)
    R = -inv(T)*T0
    P = [[i != j ? -T[i,j] * R[j,k] / T[i,i] / R[i,k] : 0 for i in 1:m, j in 1:m] for k in 1:n]
                    #(T[i,j]/(-T[i,i]))*(R[j,k]/R[i,k])
    return R, P
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

"""
QQQQ
"""
mutable struct MAPHDist
    
    #dimensions
    m::Int #number of transient phases
    n::Int #number of abosrbing states 

    #used for both parameterizations
    α::Adjoint{Float64, Vector{Float64}} #dimension at m

    #First parameterization
    T::Matrix{Float64} #dimension at mxm. The negative of the diagonal elements 
                       #of this matrix are "q" of the second parameterization
    T0::Matrix{Float64} #dimension at mxn

    #second parameterization
    R::Matrix{Float64} #dimension at mxn
    P::Vector{Matrix{Float64}} #a vector of length n, with elements at dimension mxm

    """
    Constructor using the first parameterization.
    """
    function MAPHDist(α::Adjoint{Float64, Vector{Float64}}, T::Matrix{Float64}, T0::Matrix{Float64})
        m = length(α)
        n = size(T0)[2]
        size(T) != (m, m) && error("Wrong size for T")
        size(T0) != (m, n) && error("Wrong size for T0")
        !isapprox(sum(sum(T, dims=2) + sum(T0, dims=2)),0,atol = 10e-5) && error("Non-generator matrix")

        return new(m, n, α, T, T0, R_P_from_T_T0(T, T0)...)
    end

    """
    Constructor using the second parameterization.
    """
    function MAPHDist(α::Adjoint{Float64, Vector{Float64}}, q::Vector{Float64}, R::Matrix{Float64}, P::Vector{Matrix{Float64}})
        m = length(α)
        n = size(R)[2]
        size(R) != (m, n) && error("Wrong size for R")
        length(q) != m && error("Wrong length for q")
        length(P) != n && error("Wrong length for P")
        #QQQQ todo - more checks for stochastic matrix, etc...
        T, T0 = T_T0_from_R_P_q(q, R, P)
        return new(m, n, α, T_T0_from_R_P_q(q, R, P)...,R, P)
    end
end


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
    dist.T, dist.T0 = T_T0_from_R_P_q(-diag(dist.T), dist.R, dist.P)
end

"""
QQQQ
"""
mutable struct MAPHSufficientStats
    B::Vector{Float64} #initial starts
    Z::Vector{Float64} #time spent
    M::Matrix{Float64} #transitions between transient states
    N::Matrix{Float64} #transitions between transient to abosrbing states

    MAPHSufficientStats(B::Vector{Float64}, Z::Vector{Float64}, N::Matrix{Float64}) = new(B,Z,N)
    function MAPHSufficientStats(maph::MAPHDist) 
        p, q = model_size(maph)
        new(zeros(p), zeros(p), zeros(p,p+q))
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

+(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B+ss2.B, ss1.Z+ss2.Z, ss1.N+ ss2.N)
/(ss::MAPHSufficientStats,n::Real) = MAPHSufficientStats(ss.B/n,ss.Z/n,ss.N/n)
/(ss1::MAPHSufficientStats,ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B ./ ss2.B, ss1.Z ./  ss2.Z, ss1.N ./  ss2.N)
-(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B-ss2.B, ss1.Z-ss2.Z, ss1.N - ss2.N)

"""
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