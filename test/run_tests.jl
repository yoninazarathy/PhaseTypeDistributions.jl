#This is a temporary script while developing
using Test
using Pkg
Pkg.activate((@__DIR__) * "/..")

include("../src/PhaseTypeDistributions.jl")
include("test_maph_fit.jl")
include("test_maph_init.jl")
include("test_maph_sufficient_stats.jl")
include("test_structured_ph.jl")
include("test_maph_basic.jl")

using Random
Random.seed!(2)

# @test hypoexp_test()
# @test hyperexp_test()
# @test maph_moments_and_rand_test(;N=1000000)
@test test_fit_example1(;sim_runs=10^4)
#@test full_trace_sufficient_stats_test()
#@test sufficient_stats_test()

# using Plots
#N_e, Z_e, time_vec = sufficient_stats_test(sim_runs = 10^4)

# p1 = scatter(time_vec, N_e, xlabel = "Absorbtion time", label="N errors")
# p2 = scatter(time_vec, Z_e, xlabel = "Absorbtion time", label="Z errors")
# plot(p1,p2)

# Ntest = analyze_ss_with_plots()

############################
## Playground area here....``

# fit_maph([1.2,3.2,5.2,2.3,2.1,10.8],[1,2,1,1,2,1],8)
# @test test_fit_example1()
# test_init()
