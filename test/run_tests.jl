#This is a temporary script while developing
using Test
using Pkg
Pkg.activate((@__DIR__) * "/..")
include("../src/PhaseTypeDistributions.jl")

include("test_maph_fit.jl")
include("test_maph_init.jl")
include("test_maph_sufficient_stats.jl")
include("test_structured_ph.jl")

@test hypoexp_test()
@test hyperexp_test()

############################
## Playground area here....

# test_init()

