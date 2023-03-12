"""
Returns an array of data filtered according to absorbing state being the index of the array.
"""
function absorb_filter_data(data, maph::MAPHDist)
    p, q = model_size(maph)
    filter_data = []

    for i = 1:q
        temp_data = filter(data) do obs first(obs).a == i end
        push!(filter_data, temp_data)
    end    

    #QQQQ this whole function can be just this one line ??
    # [filter((obs)->first(obs).a == i, data) for i in 1:q]
    return filter_data
end


"""
Returns an array of data filtered according to the sojourn time being less than the index of the time array.
"""
function time_filter_data(data, n::Int64)

    times = ((d)->first(d).y).(data)
    min_time = quantile(times, 0.1)
    max_time = quantile(times, 0.9)
    @show min_time, max_time
    time_vec = Array(LinRange(min_time, max_time, n)) #QQQQ - cleanup later

    filter_data = []

    for i = 1:(n-1)
        temp_data = filter(data) do obs first(obs).y ≥ time_vec[i] && first(obs).y < time_vec[i+1] end 
        push!(filter_data,temp_data)
        @show length(temp_data)
    end    

    return filter_data, time_vec
end

"""
Test the sufficient stats with full trace data
"""
function full_trace_sufficient_stats_test(;N=10^4)
    Λ₄, λ45, λ54, Λ₅ = 5, 2, 7, 10
    μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]

    maph = MAPHDist([0.5,0.5]',T_example, T0_example)

    p,q = model_size(maph)
    test_stats = MAPHSufficientStats(maph)

    for _ in 1:N
        times, states = rand(maph, full_trace = true) 
        ss = sufficient_stat_from_trajectory(maph, times, states)
        test_stats.N += ss.N
        test_stats.Z += ss.Z
        test_stats.B += ss.B
    end

    test_stats.B = test_stats.B/sum(test_stats.B)
    test_stats.Z = test_stats.Z/sum(test_stats.Z)
    test_stats.N = test_stats.N/sum(test_stats.N)

    test_α = test_stats.B

    @show test_stats.Z

    computered_intensity = (test_stats.N.*(1 ./ test_stats.Z))

    @show computered_intensity

    # @show(stats.N,test_stats.N)
    return true
end


"""
QQQQ
"""
function sufficient_stats_test(; sim_runs = 10)

    #Set heuristically the number of times steps to be the sqrt of the number of runs
    timesteps = round(Int,sqrt(sim_runs))

    Λ₄, λ45, λ54, Λ₅ = 5., 2., 7., 10.
    μ41, μ42, μ43, μ51, μ52, μ53 = 1., 1., 1., 1., 1., 1. 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
    initial_dist = [0.5,0.5]

    maph = MAPHDist(initial_dist', T_example, T0_example)

    p,q = model_size(maph)
   
    data = []
    full_trace =[]
    
    @showprogress "Simulating data" for i in 1:sim_runs
        times, states = rand(maph, full_trace = true) 
        push!(full_trace, (times,states))
        push!(data, (observation_from_full_traj(times,states),i))
    end

    absorb = absorb_filter_data(data, maph)
    absorb_data = vcat(absorb[1],absorb[2],absorb[3])
    
    #QQQQ this is for now just looking at trajectory observed in first state
    time_bin, time_vec = time_filter_data(absorb_data, timesteps)

    Z_errors = Float64[]
    N_errors = Float64[]

    #loop over all bins
    # println("\n start initialization...")
    # @showprogress "Checking sufficient stats" 


    for i in 1:length(time_bin)
        ss_i = MAPHSufficientStats[]
        for trace in full_trace[last.(time_bin[i])]
            ss = sufficient_stat_from_trajectory(maph, trace[1], trace[2])
            push!(ss_i, ss)
        end
        
        if !isempty(ss_i)
            mean_observed_ss = mean(ss_i)
            time_slice = first(data[last(last.(time_bin[i]))])[1]
            obs = first(data[last(last.(time_bin[i]))])
            @show obs
            computed_ss = sufficient_stats(obs, maph)
            # @show computed_ss.N, mean_observed_ss.N

            errs_N = (mean_observed_ss.N - computed_ss.N)#./ computed_ss
            push!(N_errors, norm(mean_observed_ss.N - computed_ss.N)/time_slice)

            # computed_ss.B ≈ initial_dist || return false

            # @show computed_ss.Z, mean_observed_ss.Z
            push!(Z_errors, norm(mean_observed_ss.Z - computed_ss.Z)/time_slice)

            # @show obs
            # sufficient_stats()
        end
    end
    
    return true
    return N_errors, Z_errors, time_vec

    # ab1 = filter(x->x.a==1,first.(data))
    # ab2 = filter(x->x.a==2,first.(data))
    # ab3 = filter(x->x.a==3,first.(data))
    
    # probs = [length(ab1),length(ab2),length(ab3)]./length(data)
    # means = [mean(first.(ab1)),mean(first.(ab2)),mean(first.(ab3))]
    # scvs = [scv(first.(ab1)),scv(first.(ab2)),scv(first.(ab3))]

    # dist = MAPHDist(10,probs,means,scvs)
    # println("finish initialization")

    # println("starting simulations")
    # fit!(dist,first.(data))
    # println("\nfinished simulations")
    
    return true

end








"""
QQQQ - An example with a deterministic path....
"""
function sufficient_stats_test3()
    
    ϵ = 0.00000001

    T_example = [-(1.0+4ϵ) 1.0-4ϵ ϵ 
                ϵ -(1.0 + 4ϵ) 1.0-4ϵ
                ϵ ϵ -(1 + 4ϵ) ]

    T0_example = [ϵ ϵ ϵ; 
                  ϵ ϵ ϵ;
                  1-4ϵ ϵ ϵ]

    
    T_example2 = [-(1.0+3ϵ) 1.0 ϵ
                 ϵ -(1.0+3ϵ) 1.0
                 ϵ ϵ -(1+3ϵ)]
    T0_example2 = [ϵ ϵ;
                   ϵ ϵ;
                   1-3ϵ ϵ]   

    display(sum(T0_example2[:,1]))
    maph = MAPHDist([0.5,0.5, 0.0]',T_example, T0_example)

    obs = (y=100.0, a=2)
    sufficient_stats(obs, maph)
end

"""
QQQQ - this is probably temporary.
"""
function analyze_ss_with_plots()

    Λ₄, λ45, λ54, Λ₅ = 5, 2, 7, 10
    μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
    initial_dist = [0.5,0.5]

    maph = MAPHDist(initial_dist', T_example, T0_example)

    SingleObs = NamedTuple{(:y, :a), Tuple{Float64, Int64}}

    for y in 0.5:0.5
        obs = (y=y, a = 1)
        ss = sufficient_stats(obs,maph)
        ss.N[1,1] #this means the expected number of jumps from the first transient state  
                        #to the first absorbing state (a=1) 
        # @show ss.N
        return ss.N
    end

end