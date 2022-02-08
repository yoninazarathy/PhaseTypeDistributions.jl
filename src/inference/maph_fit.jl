"""
Create a named tuple observation.
"""
observation_from_full_traj(times::Vector{Float64}, states::Vector{Int64}) = (y = sum(times),a = last(states))

"""
QQQQ
"""
function sufficient_stat_from_trajectory(d::MAPHDist, sojourn_times::Array{Float64}, states::Array{Int})::MAPHSufficientStats
    p,q = model_size(d)
    transient_states = (q+1):(q+p)

    ss =  MAPHSufficientStats(d)


    for s = 1:p
        if states[1] == transient_states[s]
            ss.B[s] +=1
        end
    end

    for i = 1:(length(states)-1)
        ss.Z[states[i]-q] += sojourn_times[i]
        ss.N[states[i]-q, states[i+1]] += 1
    end

    return ss
end

"""
QQQQ
"""
function stats_to_dist(maph::MAPHDist,ss::MAPHSufficientStats)::MAPHDist
    p,q = model_size(maph)
    α = ss.B'
    T = ss.N[:,(q+1):(q+p)]./ss.Z
    T0 = ss.N[:,1:q]./ss.Z
    for i = 1:p
        temp = T[i,i]
        T[i,i]= -sum(T[i,:])-sum(T0[i,:])+temp
    end
    return MAPHDist(α,T,T0)
end


"""
Fits ... QQQQ
"""
function fit_maph(times::Vector{Float64}, absorbing_states::Vector{Int}, p::Int)
    @assert length(times) == length(absorbing_states) "Vector of times and absorbing states mismatch in length"
    unique_absorbing_states = unique(absorbing_states)
    q = length(unique_absorbing_states)
    @assert minimum(unique_absorbing_states) == 1 && maximum(unique_absorbing_states) == q  "Mismatch with absorbing states"
    n = length(times)
    print("Fitting MAPH with p=$p hidden states, q=$q absorbing states, and n=$n observations")

    dist = MAPHDist(p, compute_descriptive_stats(times, absorbing_states)...)

    #QQQQ - here have the fitting loop (here is the beef!)
    #QQQQ - Zhihao continue during the week....
    
    return dist
end

"""
Computes some descriptive statistics of a data sequence
"""
function compute_descriptive_stats(times::Vector{Float64}, abosrbing_states::Vector{Int})
    #QQQQ - good to have a test function for this function (in tests)
    @assert length(times) == length(abosrbing_states) "Vector of times and absorbing states mismatch in length"
    unique_absorbing_states = unique(abosrbing_states)
    q = length(unique_absorbing_states)
    @assert minimum(unique_absorbing_states) == 1 && maximum(unique_absorbing_states) == q  "Mismatch with absorbing states"
    n = length(times)

    probs = [mean(abosrbing_states .== i) for i in 1:q]
    means = [mean(times[abosrbing_states .== i]) for i in 1:q]
    scvs = [var(times[abosrbing_states .== i]) for i in 1:q] ./ (means .^ 2)

    @show probs

    return probs, means, scvs
end


"""
Fits ... QQQQ
"""
function fit!(maph::MAPHDist,data::MAPHObsData)::MAPHDist
    #EM Loop

    p,q = model_size(maph)


    for k in 1:100
        ss = sufficient_stats(data[1],maph)
        for i in 2:10^2
            ss = ss+sufficient_stats(data[i],maph)
        end
        ss = ss/10^2
                maph = stats_to_dist(maph,ss)
    end

    return maph
end
