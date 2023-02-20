"""
QQQQ
"""
mutable struct PHDist
    α::Adjoint{Float64, Vector{Float64}}
    T::Matrix{Float64}
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
    T::Matrix{Float64} #dimension at mxm. The diagonal elements of this matrix are "q" of the second parameterization
    T0::Matrix{Float64} #dimension at mxn

    #second parameterization
    R::Matrix{Float64} #dimension at mxn
    P::Vector{Matrix{Float64}} #a vector of length n, with elements at dimension mxm

    function MAPHDist(α::Adjoint{Float64, Vector{Float64}}, T::Matrix{Float64},T0::Matrix{Float64})
        @assert size(T)[1] == size(T0)[1]
        @assert isapprox(sum(sum(T, dims=2) + sum(T0, dims=2)),0,atol = 10e-5)
        return new(α, T, T0)
    end
end

"""
QQQQ -doc to do
"""
function update_params_1to2(dist::MAPHDist)

end

"""
QQQQ -doc to do
"""
function update_params_2to1(dist::MAPHDist)

end

"""
QQQQ
"""
mutable struct MAPHSufficientStats
    B::Vector{Float64} #initial starts
    Z::Vector{Float64} #time spent
    N::Matrix{Float64} #transitions
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