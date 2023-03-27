"""
Create a named tuple observation.
"""
observation_from_full_traj(times::Vector{Float64}, states::Vector{Int64}) = (y = sum(times),a = last(states))

"""
QQQQ
"""
function sufficient_stat_from_trajectory(d::MAPHDist, sojourn_times::Array{Float64}, states::Array{Int})::MAPHSufficientStats
    m, n = model_size(d)
    ss =  MAPHSufficientStats(d)

    #states are indexed as "absorbing first" so subtract number of absorbing for correct index (__ - n)

    #Set B
    ss.B[states[1] - n] = 1  

    #Set Z - Loop over transient states
    for i = 1:(length(states) - 1)
        ss.Z[states[i] - n] += sojourn_times[i]
    end

    #Set M - Loop over transitions between transient states e.g. with 4 states, (Trans, Trans, Trans, Abs), we will have 2 transitions
    #if length(states) = 2 , then (trans, abs) and never loop
    for i = 1:(length(states) - 2)
        ss.M[states[i] - n, states[i+1] - n] += 1
    end

    #Set N
    ss.N[states[end - 1] - n, states[end]] = 1

    return ss
end

"""
`ss` is the average expected sufficient statistics. 
"""
function maximum_likelihood_estimate_second_parameter(ss::MAPHSufficientStats, maph::MAPHDist)
    m, n = model_size(maph)
    α̂ = max.(ss.B, 0)
    N_total = sum(ss.N, dims = 2) + sum(ss.M, dims=2)
    q̂ = [N_total[i]/ss.Z[i] for i = 1:m]
    absorb = absorption_probs(maph)
    R̂ = [absorb[k]*ss.B[i]/ss.B[i] for i =1:m, k = 1:n]
    P̂ = [[ss.M[i,j]/N_total[i] for i = 1:m, j = 1:m] for k = 1:n]
    
    @show "running34324"
    @show "hi"
    maph.α = α̂'
    maph.T = Diagonal(q̂)
    maph.R = R̂
    maph.P = P̂
    update_params_2to1!(maph)
    return maph
end



#Maybe rename this function to "maximization step..." or similar
function maximum_likelihood_estimate(p::Int, q::Int, ss::MAPHSufficientStats)

    #QQQQ pickup here next time.... (from eq 16)

    α_next = max.(ss.B,0)
    t_next = max.(ss.N[:,1:q] ./ ss.Z,0)

    t_next[isnan.(t_next)] .= 0

    T_next = zeros(p,p) 
    T0_next = zeros(p,q)

    for i = 1:p
        T_next[i,:] = max.(ss.N[i,(q+1):(q+p)]./ss.Z[i],0)
        T_next[i,isnan.(T_next[i,:])].=0
        T_next[i,i] = -(sum(t_next[i,:]) + sum(T_next[i,:])) #-(t_next[i]+sum(T_next[i,:]))
        T0_next[i,:] = max.(ss.N[i,1:q]./ss.Z[i],0)
        T0_next[i,isnan.(T0_next[i,:])].=0
    end
    #################

    T = ss.N[:,(q+1):(q+p)]./ss.Z
    T0 = ss.N[:,1:q]./ss.Z
    for i = 1:p
        T[i,i] += -sum(T[i,:])-sum(T0[i,:])
    end


    return α_next, T_next, T0_next
end


"""
Fits ... QQQQ
"""
function fit_maph(  times::Vector{Float64}, 
                    absorbing_states::Vector{Int}, 
                    p::Int; 
                    max_iter::Int = 100, 
                    ω::Float64 = 10)::MAPHDist

    length(times) != length(absorbing_states) && error("Vector of times and absorbing states mismatch in length")
    unique_absorbing_states = unique(absorbing_states)
    q = length(unique_absorbing_states)
    (minimum(unique_absorbing_states) != 1 || maximum(unique_absorbing_states) != q) && error("Mismatch with absorbing states")
    n = length(times)
    print("Fitting MAPH with p=$p hidden states, q=$q absorbing states, and n=$n observations")

    #compute descriptive stats
    ds = compute_descriptive_stats(times, absorbing_states)
    
    #construct a moment based initialization
    dist = MAPHDist(p, ds...; ω) #Moment based 
    
    iter = 1
    while iter < max_iter

        continue; #return to top QQQQ

        #QQQQ - refactor next couple of lines...

        ## Expectation steps
        computed_ss = MAPHSufficientStats[]
        for i = 1:n
            obs = (y = times[i], a = absorbing_states[i])
            ss_i = sufficient_stats(obs, dist)
            push!(computed_ss, ss_i)
        end

        ## Maximisation steps
        mean_ss = mean(computed_ss)

        α_next, T_next, T0_next = maximum_likelihood_estimate(p, q, mean_ss) 
        #remove numerical instabilities 

        #@show T_next, T0_next

        α_next = max.(α_next, 0) #QQQQ check why - maybe goes slightly negative
        α_next /= sum(α_next)
        dist = MAPHDist(α_next', T_next, T0_next)

        iter +=1
    end
 
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

