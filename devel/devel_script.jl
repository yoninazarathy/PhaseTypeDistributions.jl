using Pkg
Pkg.activate((@__DIR__) * "/..")

include("../src/PhaseTypeDistributions.jl")

using .PhaseTypeDistributions

using ProgressMeter

#This is now a scratch pad for development.

Λ₄, λ45, λ54, Λ₅ = 5., 2., 2., 5.
μ41, μ42, μ43, μ51, μ52, μ53 = 1., 1., 1., 1., 1., 1. 
T_example = [-Λ₄ λ45; λ54 -Λ₅]
T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
initial_dist = [0.5,0.5]

maph = MAPHDist(initial_dist', T_example, T0_example)
init_maph = deepcopy(maph)

m,n = model_size(maph)
   
data = []
sufficient_stats_data = []

sim_runs = 1000


@showprogress "Simulating data (and computing sufficient states on data)" for i in 1:sim_runs
    times, states = rand(maph, full_trace = true) 
    push!(sufficient_stats_data, sufficient_stat_from_trajectory(maph, times, states) )
    push!(data, (observation_from_full_traj(times, states),i))
end

computed_stats = []


@showprogress "Estimating mean sufficient stats on data" for i in 1:length(data)
    obs = first(data[i])
    computed_ss = sufficient_stats(obs, maph)
    push!(computed_stats, computed_ss)
end
average_expected_ss = mean(computed_stats)
@show average_expected_ss
mle = maximum_likelihood_estimate_second_parameter(average_expected_ss, maph)
@show mean(mle), mean(init_maph)


# using Plots
#N_e, Z_e, time_vec = sufficient_stats_test(sim_runs = 10^4)

# p1 = scatter(time_vec, N_e, xlabel = "Absorbtion time", label="N errors")
# p2 = scatter(time_vec, Z_e, xlabel = "Absorbtion time", label="Z errors")
# plot(p1,p2)

# Ntest = analyze_ss_with_plots()

############################
## Playground area here....``

# fit_maph([1.2,3.2,5.2,2.3,2.1,10.8],[1,2,1,1,2,1],8)
