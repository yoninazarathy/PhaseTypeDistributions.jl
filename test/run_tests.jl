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

# @test hypoexp_test()
# @test hyperexp_test()
# @test maph_moments_and_rand_test()
# @test maph_init_test()
@test full_trace_sufficient_stats_test()
############################
## Playground area here....``


# test_init()