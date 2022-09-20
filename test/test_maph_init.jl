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
        @show mean(dist), desired_means)
        # @show absorption_probs(dist)
    end
    return true
end
