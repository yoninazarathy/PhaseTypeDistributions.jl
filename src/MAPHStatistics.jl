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
    a::Int

    "absorbing time"
    y::Float64

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
$(TYPEDEF)
The struct which representing the statistics to fit the MAPH distribtuion 
$(TYPEDFIELDS)
"""
struct MAPHSufficientStats
    "initial starts"
    B::Vector{Float64} 

    "time spent"
    Z::Vector{Float64} 

    "transitions between transient states"
    M::Matrix{Float64} 

    "transitions between transient to abosrbing states"
    N::Matrix{Float64} 

end



+(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B+ss2.B, ss1.Z+ss2.Z, ss1.M+ ss2.M, ss1.N+ss2.N)
/(ss::MAPHSufficientStats,n::Real) = MAPHSufficientStats(ss.B/n, ss.Z/n, ss.M/n, ss.N/n)
/(ss1::MAPHSufficientStats,ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B ./ ss2.B, ss1.Z ./  ss2.Z, ss1.M ./  ss2.M,  ss1.N ./  ss2.N)
-(ss1::MAPHSufficientStats, ss2::MAPHSufficientStats) = MAPHSufficientStats(ss1.B-ss2.B, ss1.Z-ss2.Z, ss1.M - ss2.M, ss1.N - ss2.N)
*(n::Real, ss::MAPHSufficientStats) = MAPHSufficientStats(ss.B *n, ss.Z*n, ss.M*n, ss.N*n)

function very_crude_c_solver(y::Float64, i::Int, j::Int, k::Int, maph::MAPHDist)
    quadgk(u -> (maph.α * exp(maph.T*u))[i] * (exp(maph.T*(y-u))*maph.D[:,k])[j] , 0, y, rtol=1e-8) |> first
end

"""
$(METHODLIST)
Compute the expected value of the sufficient stats
"""
function compute_sufficient_stats(observation::SingleObservation, 
                            maph::MAPHDist; 
                            c_solver = very_crude_c_solver)

    m, n = model_size(maph)
    
    a(y::Float64) = maph.α * exp(maph.T*y)
    b(y::Float64, k::Int) = exp(maph.T*y) * maph.D[:,k]
    c(y::Float64, i::Int, j::Int, k::Int) = very_crude_c_solver(y, i, j, k, maph)

    EB(y::Float64, i::Int, k::Int) = maph.α[i] * b(y, k)[i] / reduce(vcat, (maph.α * b(y, k)))
    EZ(y::Float64, i::Int, k::Int) = c(y, i, i, k) / reduce(vcat, (maph.α * b(y,k)))
    ENT(y::Float64, i::Int, j::Int, k::Int) = i != j ? maph.T[i,j] .* c(y, i, j, k) / reduce(vcat, (maph.α * b(y,k))) : 0
    ENA(y::Float64, i::Int, j::Int, k::Int) = j == k ? a(y)[i] * maph.D[i,k] / reduce(vcat, (maph.α * b(y,k))) : 0 

    B = map(i -> EB(observation.y, i, observation.a - m), 1:m)
    Z = map(i -> EZ(observation.y, i, observation.a - m), 1:m)

    M = reduce(hcat, map(i ->  map(j -> ENT(observation.y, i, j, observation.a - m), 1:m), 1:m))
    N = reduce(hcat, map(j -> map(i -> ENA(observation.y, i, j, observation.a -m) , 1:m ), 1:n))

    return MAPHSufficientStats(B, Z, M ,N)
end

# function weighted_maph_sufficient_stats(stats::MAPHSufficientStats, absorbing_prob::Vector{Float64})
#     B = stats.B
#     Z = stats.Z
#     M = 


# end



function compute_expected_stats(all_obs::Vector{SingleObservation}, maph::MAPHDist; c_solver = very_crude_c_solver, lower_quantile = 0.01, upper_quantile = 0.99)
    m, n = model_size(maph)
    all_absorbing_states = map(obs -> obs.a, all_obs)
    all_absorbing_times = map(obs -> obs.y, all_obs)
    min_time = quantile(all_absorbing_times, lower_quantile)
    max_time = quantile(all_absorbing_times, upper_quantile)

    #filter out those ones with abnormal time
    filtered_obs = filter(obs -> (min_time < obs.y < max_time),  all_obs)

    unique_absorbing_state = sort(unique(all_absorbing_states))
    #sort the obs data to different absorb states
    data_sorted_by_states = map(k-> filter(obs -> obs.a == k, filtered_obs), unique_absorbing_state)
    
    absorption_counts = countmap(all_absorbing_states)
    absorbing_prob = map(i -> absorption_counts[i] / length(all_absorbing_states), unique_absorbing_state)
    expected_stats_by_different_states = map(data_by_absorb_state -> mean(compute_sufficient_stats.(data_by_absorb_state, Ref(maph))), data_sorted_by_states)
    weighted_stats = map(k -> absorbing_prob[k] * expected_stats_by_different_states[k], eachindex(absorbing_prob))
    

    return sum(weighted_stats)


end


