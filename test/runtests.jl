# using PhaseTypeDistributions

using Test
using ProgressMeter
using Random
using MAPHDistributions
using Plots
using GR
using Distributions
using LinearAlgebra

Random.seed!(2)

# include("../src/MAPH.jl")
# include(".../src/moment_function")
# include(".../")

# Λ₄, λ45, λ54, Λ₅ = 15.0, 5.0, 7.0, 16.0
# μ41, μ42, μ43, μ51, μ52, μ53 = 4.0, 3.0, 3.0, 1.0, 7.0, 1.0 
# T = [-Λ₄ λ45; λ54 -Λ₅]
# D = [μ41 μ42 μ43; μ51 μ52 μ53]
# α = [0.5 0.5]

# m = 50
# n = 5
# P = rand(m, m)
# P = P - Diagonal(P)
# P = P ./ sum(P, dims = 2) .* rand(m)
# R = rand(m, n)
# R = R ./ sum(R, dims = 2)
# q = rand(Exponential(1), m)

# maph = MAPH_constructor(α, q, R, P)

maph = maph_random_parameters(4 , 3)

@show maph.D
# maph2 = deepcopy(maph)


# maph = MAPH_constructor(α, T, D)

# is_valid_R_P(maph.R, maph.P)


# maph = maph_random_parameters(3, 5)

# maph2 = deepcopy(maph)

# ys = collect(0:0.001:2)

# # pdf1 = sub_distribution(maph2, 1, ys)

# # pdfs = sub_distribution.(Ref(maph2), collect(1:3), Ref(ys))

# all_obs = map(n -> rand(maph), 1:50)

# # stats = compute_sufficient_stats(all_obs[1], maph)

# data_length, stats, s_stats = E_step!(all_obs, maph)

# M_step!(data_length, stats, s_stats, maph)




# M_step!(data_length, stats, sum_stats_per_state, maph)
# compute_sorted_stats(all_obs, maph)
# EM_fit!(all_obs, maph, 1)
# maph_fit = maph_initialization(all_obs, 10; ω = 30.0,  θ = 5.0)
# # # EM_fit!(all_obs, maph_fit, 1



# @show get_emperical_absorb_prob(all_obs)

# for k = 1:20
#     Maximization_step!(all_obs, maph)
# end

# pdfs2 = sub_distribution.(Ref(maph), collect(1:3), Ref(ys))

