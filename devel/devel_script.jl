using Pkg
Pkg.activate((@__DIR__) * "/..")
using Revise
# include("../src/PhaseTypeDistributions.jl")
#import .PhaseTypeDistributions
using PhaseTypeDistributions
using ProgressMeter

#This is now a scratch pad for development

###############
# Develop MLE #
###############
if true

    probs, means, scvs = [0.458, 0.2722, 0.2698], [0.2706401195766266, 0.27505749846070254, 0.2872995493009737], [1.2, 1.028862135296857, 1.5]
    total_mean = probs'means
    total_scv = probs'scvs
    # @show total_mean, total_scv
    d = MAPHDist(10, probs, means, scvs; ω = 100.);
    L = 100
    obs = [rand(d) for _ in 1:L]
    # @show mean(first.(obs))
    # @show scv(first.(obs))
    # @show [mean(last.(obs).==k) for k in 1:3]
 
    ss = [sufficient_stats(ob, d) for ob in obs]
    @show mean([ss1.Z for ss1 in ss])

    

end



if false
    Λ₄, λ45, λ54, Λ₅ = 10., 5., 2., 5.
    μ41, μ42, μ43, μ51, μ52, μ53 = 3., 1., 1., 1., 1., 1. 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
    initial_dist = [0.5,0.5]

    maph = MAPHDist(initial_dist', T_example, T0_example)
    init_maph = deepcopy(maph)

    m,n = model_size(maph)
    
    data = []
    sufficient_stats_data = []

    sim_runs = 10_000

    # times, states = rand(maph, full_trace = true)
    # @show (times, states)

    @showprogress "Simulating data (and computing sufficient states on data)" for i in 1:sim_runs
        times, states = rand(maph, full_trace = true) 
        #push!(sufficient_stats_data, sufficient_stat_from_trajectory(maph, times, states) )
        push!(data, (observation_from_full_traj(times, states),i))
    end 


    times_data = [data[i][1].y for i in 1:length(data)]
    states_data = [data[i][1].a for i in 1:length(data)]
    
    display(compute_descriptive_stats(times_data, states_data))

    # filtered_data = absorb_filter_data(data,maph)

    # compute_descriptive_stats

    # # computed_stats = []

    #filtered_data = absorb_filter_data(data,maph)

    # total_stats = []

    # for data in filtered_data
    #     computed_stats = []
    #     @showprogress "Estimating mean sufficient stats on data"  for i in 1:length(data)
    #         obs = first(data[i])
    #         computed_ss = sufficient_stats(obs, maph)
    #         push!(computed_stats, computed_ss)
    #     end
    #     push!(total_stats, computed_stats)
    # end
    
    # @show reduce(vcat,total_stats)[1]
    
    
    # @showprogress "Estimating mean sufficient stats on data" for i in 1:length(data)
        #     obs = first(data[i])
    #     computed_ss = sufficient_stats(obs, maph)
    #     push!(computed_stats, computed_ss)
    # end
    
    # @show absorb_filter_data(data,maph)[1]
    
    # average_expected_ss = mean(computed_stats)
    # @show average_expected_ss
    # mle = maximum_likelihood_estimate_second_parameter(total_stats, maph)
    
    # @show mle
    # @show mean(mle), mean(init_maph)
    # @show scv(mle), scv(init_maph)
end