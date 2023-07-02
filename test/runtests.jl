using PhaseTypeDistributions
using Test
using ProgressMeter

include("../src/PhaseTypeDistributions.jl")
include("test_maph_fit.jl")
include("test_maph_init.jl")
include("test_maph_sufficient_stats.jl")
include("test_structured_ph.jl")
include("test_maph_basic.jl")

using Random
Random.seed!(2)


# @test maph_type_construction_tests()
# @test hypoexp_test()
# @test hyperexp_test()
# @test maph_moments_and_rand_test(;N=1000000)
# @test test_fit_example1(;sim_runs=10^2)
# @test test_maph_init();
# @test test_maph_perturbation()
# @test full_trace_sufficient_stats_test()
#@test sufficient_stats_test()
