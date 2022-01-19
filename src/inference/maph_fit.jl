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
