include("MAPH.jl")
include("moment_function.jl")

"""
$(TYPEDFIELDS)
The struct representing the CTMC with an underlying maph distribution
$(TYPEDEF)
"""

struct ContinuousTimeMarkovChain
    "the transient states"
    transient_states:: Vector{Int}

    "the absorbing states"
    absorbing_states:: Vector{Int}

    "transient probabilities"
    P_λ ::Matrix{Float64}

    "absorbing probabilities"
    P_μ:: Matrix{Float64}

    "jumping chain"
    P_jump :: Matrix{Float64}

    "diagonal rate"
    Λ:: Vector{Float64}

end

function MAPH_TO_CTMC(maph::MAPHDist)
    m, n  = model_size(maph)
    transient_states = 1:m
    absorbing_states = (m+1):(m+n)
    P_λ = I - inv(Diagonal(maph.T))*maph.T
    P_μ = maph.R
    P1 = hcat(P_λ, P_μ)
    P2 = hcat(zeros(n,m), Matrix(Diagonal(ones(n))))
    P_jump =[P1; P2]
    Λ = vcat(-diag(maph.T), zeros(n))
    return ContinuousTimeMarkovChain(transient_states, absorbing_states, P_λ, P_μ, P_jump, Λ)
end



"""
$(TYPEDFIELDS)
The struct which record time and absorbing state of a CTMC process 
$(TYPEDEF)
"""
struct SingleObservation
    "absorbing state"
    a_state::Int

    "absorbing time"
    a_time::Float64

    "sojourn_times"
    sojourn_times::Vector{Float64}

    "all_states"
    all_states::Vector{Int}
end

function rand(maph::MAPHDist; full_trace = true)
    m, n = model_size(maph)
    all_states = 1:(m+n)
    ctmc = MAPH_TO_CTMC(maph) 

    states = Int[]
    sojourn_times = Float64[]
    
    state = sample(ctmc.transient_states, weights(maph.α))
    t = 0.0
    @info "now generating random paths"
    while state ∈ ctmc.transient_states
        sojourn_time = rand(Exponential(1/ctmc.Λ[state]))
        if full_trace
            push!(states, state)
            push!(sojourn_times, sojourn_time)
        end
        t += sojourn_time
        state = sample(all_states, weights(ctmc.P_jump[state, :]))
    end
    push!(states, state)
    
    return SingleObservation(state, t, sojourn_times, states)
end


"""
$(TYPEDFIELDS)
The mutable struct which representing the statistics to fit the MAPH distribtuion 
$(TYPEDEF)
"""
mutable struct MAPHSufficientStats
    "initial starts"
    B::Vector{Float64} 

    "time spent"
    Z::Vector{Float64} 

    "transitions between transient states"
    M::Matrix{Float64} 

    "transitions between transient to abosrbing states"
    N::Matrix{Float64} 

    MAPHSufficientStats(B::Vector{Float64}, Z::Vector{Float64}, M::Matrix{Float64},  N::Matrix{Float64}) = new(B,Z,M,N)
    function MAPHSufficientStats(maph::MAPHDist) 
        m, n = model_size(maph)
        @assert (m > 0) && (n > 0) " need at least one transient and absorbing phases"
        new(zeros(m), zeros(m), zeros(m,m), zeros(m,n))
    end
end

+(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B+ss2.B, ss1.Z+ss2.Z, ss1.M+ ss2.M, ss1.N+ss2.N)
/(ss::MAPHSufficientStats,n::Real) = MAPHSufficientStats(ss.B/n, ss.Z/n, ss.M/n, ss.N/n)
/(ss1::MAPHSufficientStats,ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B ./ ss2.B, ss1.Z ./  ss2.Z, ss1.M ./  ss2.M,  ss1.N ./  ss2.N)
-(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B-ss2.B, ss1.Z-ss2.Z, ss1.M - ss2.M, ss1.N - ss2.N)
    




