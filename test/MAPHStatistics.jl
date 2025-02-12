Λ₄, λ45, λ54, Λ₅ = 15.0, 5.0, 7.0, 16.0
μ41, μ42, μ43, μ51, μ52, μ53 = 4.0, 3.0, 3.0, 1.0, 7.0, 1.0 
T = [-Λ₄ λ45; λ54 -Λ₅]
D = [μ41 μ42 μ43; μ51 μ52 μ53]
α = [0.5, 0.5]

maph = MAPHDist(α, T, D)
single_ob = rand(maph)
@show compute_sufficient_stats(single_ob, maph)

all_obs = [rand(maph) for _ in 1:10]
@show MAPHDistributions.compute_sorted_stats( all_obs, maph)

# a(maph::MAPHDist, y::Real) = maph.α' * exp(maph.T*y)
# b(maph::MAPHDist, y::Real, k::Int) = exp(maph.T*y) * maph.D[:,k]
# c(maph::MAPHDist, y::Real, i::Int, j::Int, k::Int) = very_crude_c_solver(y, i, j, k, maph)

# @show a(maph, 0.5)

# @show b(maph, 0.5, 1)

# @show c(maph, 0.5, 1,1 ,1 )