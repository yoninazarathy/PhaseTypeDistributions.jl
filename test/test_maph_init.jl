"""
QQQQ
"""
function test_maph_init()
    desired_absorbptions = [0.3, 0.2, 0.4, 0.08, 0.02]
    desired_means = [2.0, 3.1, 2.3, 4.5, 0.2]
    desired_scvs = [1.3, 0.7, 2.5, 2.4, 0.5]

    for p in 30:30
        dist  = MAPHDist(p, desired_absorbptions, desired_means, desired_scvs)
        # @show p
        # @show mean(dist), desired_means)
        # @show absorption_probs(dist)
    end
    return true
end

function maph_type_construction_tests()
    Λ₄, λ45, λ54, Λ₅ = 15.0, 5.0, 7.0, 16.0
    μ41, μ42, μ43, μ51, μ52, μ53 = 4.0, 3.0, 3.0, 1.0, 7.0, 1.0 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
    initial_dist = [0.5, 0.5]
    dist = MAPHDist(initial_dist', T_example, T0_example)
    # dist2 = MAPHDist(initial_dist', -diag(T_example), dist.R, dist.P)

    @show dist.P

    # @show dist.R
    # @show dist.R*ones(3)
    # @show initial_dist'*dist.R
    # @show absorption_probs(dist)
    # @show dist.T0
    # @show dist2.T0
    # @show mean(dist)
    # @show mean(dist2)
    return true
end