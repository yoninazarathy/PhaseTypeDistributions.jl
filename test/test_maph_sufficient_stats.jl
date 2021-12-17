"""
Returns an array of data filtered according to absorbing state being the index of the array.
"""
function absorb_filter_data(data, maph::MAPHDist)
    p, q = model_size(maph)
    filter_data = []

    for i = 1:q
        temp_data = filter(data) do obs first(obs).a == i end
        push!(filter_data,temp_data)
    end    
    
    return filter_data
end


"""
QQQQ
"""
function time_filter_data(data, n::Int64)

    max_time = maximum(((d)->first(d).y).(data) )#[data[i].y for i = 1:length(data)])
    time_vec = Array(LinRange(0, max_time, n)) #QQQQ - cleanup later

    filter_data = []

    for i = 1:(n-1)
        temp_data = filter(data) do obs first(obs).y ≥ time_vec[i] && first(obs).y < time_vec[i+1] end 
        push!(filter_data,temp_data)
    end    

    return(filter_data)
end

"""
QQQQ
"""
function sufficient_stats_test1()
    Λ₄, λ45, λ54, Λ₅ = 5, 2, 7, 10
    μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]

    maph = MAPHDist([0.5,0.5]',T_example, T0_example)
    stats = MAPHSufficientStats(maph)
    # @show maph
    # @show model_size(maph)

    data = [(y=2.3,a=2),(y=5.32,a=1),(y=15.32,a=2)]
    # update_sufficient_stats(maph, data,stats)

    # @show stats

    # QQQQ update_sufficient_stats(maph,data,stats)

    test_stats = MAPHSufficientStats(maph)

    for _ in 1:10^5
        times, states = rand(maph, full_trace = true) 
        ss = sufficient_stat_from_trajectory(maph, times, states)
        test_stats.N += ss.N
        test_stats.Z += ss.Z
        test_stats.B += ss.B
    end

    test_stats.B = test_stats.B/sum(test_stats.B)
    test_stats.Z = test_stats.Z/sum(test_stats.Z)
    test_stats.N = test_stats.N/sum(test_stats.N)

    # @show(stats.N,test_stats.N)
    @show test_stats
end


"""
QQQQ
"""
function sufficient_stats_test2()
    Λ₄, λ45, λ54, Λ₅ = 5, 2, 7, 10
    μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]

    maph = MAPHDist([0.5,0.5]', T_example, T0_example)
   
    data = []
    full_trace =[]
    
    println("starting generating data")
    for i in 1:2*10^4
        times, states = rand(maph, full_trace = true) 
        push!(full_trace, (times,states))
        push!(data, (observation_from_full_traj(times,states),i))
        i % 10^5 == 0 && print(".")
    end
    println("\nfinished generating data")

    absorb = absorb_filter_data(data, maph)
    time_bin = time_filter_data(absorb[1], 1000)

    #loop over all bins
    println("\n start initialization...")
    for i in 1:length(time_bin)
        ss_i = MAPHSufficientStats[]
        for trace in full_trace[last.(time_bin[i])]
            ss = sufficient_stat_from_trajectory(maph, trace[1], trace[2])
            push!(ss_i, ss)
        end
        
        if !isempty(ss_i)
            mean_observed_ss = mean(ss_i)
            obs = first(data[last(last.(time_bin[i]))])
            computed_ss = sufficient_stats(obs, maph)
            errs_N = (mean_observed_ss.N - computed_ss.N) ./ computed_ss.N #./ computed_ss
            # @show obs
            # sufficient_stats()
        end
    end
    
    ab1 = filter(x->x.a==1,first.(data))
    ab2 = filter(x->x.a==2,first.(data))
    ab3 = filter(x->x.a==3,first.(data))
    
    probs = [length(ab1),length(ab2),length(ab3)]./length(data)
    means = [mean(first.(ab1)),mean(first.(ab2)),mean(first.(ab3))]
    scvs = [scv(first.(ab1)),scv(first.(ab2)),scv(first.(ab3))]

    dist = MAPHDist(10,probs,means,scvs)
    println("finish initialization")

    println("starting simulations")
    fit!(dist,first.(data))
    println("\nfinished simulations")
    
    return dist

    # @show first(data)
    # # @show m[1]
    # @show maximum(data.y)
    # @show maximum[data[i].y for i = 1:length(data)]
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
