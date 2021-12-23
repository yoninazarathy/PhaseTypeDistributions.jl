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

@test test_hypoexp()
@test test_hyperexp()
@test test_maph_moments_and_rand()
@test test_maph_init()

############################
## Playground area here....``

# test_init()

