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
    # eps = 0.0000001
    Λ₄, λ45, λ54, Λ₅ = 15.0, 5.0, 7.0, 16.0
    # Λ₄, λ45, λ54, Λ₅ = 15.0+3eps, eps, eps, 16.0+3eps
    μ41, μ42, μ43, μ51, μ52, μ53 = 4.0, 3.0, 3.0, 1.0, 7.0, 1.0 
    # μ41, μ42, μ43, μ51, μ52, μ53 = 15.0, eps, eps, eps, eps, 16.0 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
    initial_dist = [0.5, 0.5]
    dist = MAPHDist(initial_dist', T_example, T0_example)

    # println("T:")
    # display(dist.T)

    # println("T0:")
    # display(dist.T0)

    # println("R:")
    # display(dist.R)

    # println("P:")
    # for (k,p) in enumerate(dist.P)
    #     println("k = $k:")
    #     display(p)
    # end

    dist2 = MAPHDist(initial_dist', -diag(T_example), dist.R, dist.P)

    # @show dist
    # @show dist2

    # @show dist.P

    # @show dist.R
    # @show dist.R*ones(3)
    # @show initial_dist'*dist.R
    # @show absorption_probs(dist)
    # @show dist.T0
    # # @show dist2.T0
    # @show mean(dist)
    # @show mean(dist2)

    # @show scv(dist)
    # @show scv(dist2)

    return (mean(dist) ≈ mean(dist)) && (scv(dist) ≈ scv(dist))
end