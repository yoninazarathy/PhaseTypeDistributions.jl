using Pkg
Pkg.activate((@__DIR__) * "/..")

include("../src/PhaseTypeDistributions.jl")

using .PhaseTypeDistributions

#This is now a scratch pad for development.

# Λ₄, λ45, λ54, Λ₅ = 5., 2., 2., 5.
# μ41, μ42, μ43, μ51, μ52, μ53 = 1., 1., 1., 1., 1., 1. 
# T_example = [-Λ₄ λ45; λ54 -Λ₅]
# T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
# initial_dist = [0.5,0.5]

# maph = MAPHDist(initial_dist', T_example, T0_example)
# @show mean(maph)

# using Plots
#N_e, Z_e, time_vec = sufficient_stats_test(sim_runs = 10^4)

# p1 = scatter(time_vec, N_e, xlabel = "Absorbtion time", label="N errors")
# p2 = scatter(time_vec, Z_e, xlabel = "Absorbtion time", label="Z errors")
# plot(p1,p2)

# Ntest = analyze_ss_with_plots()

############################
## Playground area here....``

# fit_maph([1.2,3.2,5.2,2.3,2.1,10.8],[1,2,1,1,2,1],8)
