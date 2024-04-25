# using PhaseTypeDistributions

using Test
using ProgressMeter
using Random
using MAPHDistributions
using Plots
using GR

Random.seed!(2)

# include("../src/MAPH.jl")
# include(".../src/moment_function")
# include(".../")

Λ₄, λ45, λ54, Λ₅ = 15.0, 5.0, 7.0, 16.0
μ41, μ42, μ43, μ51, μ52, μ53 = 4.0, 3.0, 3.0, 1.0, 7.0, 1.0 
T = [-Λ₄ λ45; λ54 -Λ₅]
D = [μ41 μ42 μ43; μ51 μ52 μ53]
α = [0.5 0.5]

maph = MAPH_constructor(α, T, D )
maph2 = deepcopy(maph)

ys = collect(0:0.001:2)

# pdf1 = sub_distribution(maph2, 1, ys)

pdfs = sub_distribution.(Ref(maph2), collect(1:3), Ref(ys))

all_obs = map(n -> rand(maph), 1:100)


maph2 = maph_initialization(all_obs, 8)

Maximization_step!(all_obs, maph)

# @show get_emperical_absorb_prob(all_obs)

# for k = 1:20
#     Maximization_step!(all_obs, maph)
# end

# pdfs2 = sub_distribution.(Ref(maph), collect(1:3), Ref(ys))


# plot(ys, pdfs[1],  label = "original maph absorbed in state $(1+3)")


# gr()
# plot(ys, pdfs[1], label = "original maph absorbed in state 3")
# for i = 2:3
#     display(plot!(ys, pdfs[i],  label = "original maph absorbed in state $(i+2)"))
# end

# for i = 1:3
#     display(plot!(ys, pdfs2[i],  label = "fitted maph absorbed in state $(i+2)", ls = :dot, lw = 2))
# end
# include("test_maph_fit.jl")
# include("test_maph_init.jl")
# include("test_maph_sufficient_stats.jl")
# include("test_structured_ph.jl")
# include("test_maph_basic.jl")




# @test maph_type_construction_tests()
# @test hypoexp_test()
# @test hyperexp_test()
# @test maph_moments_and_rand_test(;N=1000000)
# @test test_fit_example1(;sim_runs=10^2)
# @test test_maph_init();
# @test test_maph_perturbation()
# @test full_trace_sufficient_stats_test()
#@test sufficient_stats_test()

