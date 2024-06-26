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

# pdfs = sub_distribution.(Ref(maph2), collect(1:3), Ref(ys))

all_obs = map(n -> rand(maph), 1:200)
maph_fit = maph_initialization(all_obs, 2; ω = 30.0,  θ = 5.0)
EM_fit!(all_obs, maph_fit)



# @show get_emperical_absorb_prob(all_obs)

# for k = 1:20
#     Maximization_step!(all_obs, maph)
# end

# pdfs2 = sub_distribution.(Ref(maph), collect(1:3), Ref(ys))

