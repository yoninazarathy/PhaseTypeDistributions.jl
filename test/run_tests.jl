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

<<<<<<< HEAD
@test hypoexp_test()
@test hyperexp_test()
@test maph_moments_and_rand_test()
@test init_test()
@test test_full_trace_sufficient_stats()
=======
@test test_hypoexp()
@test test_hyperexp()
@test test_maph_moments_and_rand()
@test test_maph_init()

>>>>>>> d71935519d27b6b8e7403df7fab131c832b45a11
############################
## Playground area here....``


# test_init()

