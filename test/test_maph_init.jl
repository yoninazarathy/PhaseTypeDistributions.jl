"""
QQQQ
"""
function test_maph_init()
    desired_absorbptions = [0.3, 0.2, 0.4, 0.08, 0.02]
    desired_means = [2.0, 3.1, 2.3, 4.5, 0.2]
    desired_scvs = [1.3, 1.7, 2.5, 0.4, 0.05]

    p = 10
    dist  = MAPHDist(p, desired_absorbptions, desired_means, desired_scvs)


    # for p in 3:5
    #     dist  = MAPHDist(p, desired_absorbptions, desired_means, desired_scvs)
    #     @show p
    #     @show mean(dist)
    #     @show absorption_probs(dist)
    # end
    # data = [rand(dist) for _ in 1:10^4]
    # pi_est = [count((x)->x.a == i, data)/length(data) for i in 1:5]
    # μ_est = [mean(first.(filter((x)->x.a == i, data))) for i in 1:5] 

  


    # for p in [2, 5, 10, 30]
    #     dist = MAPHDist(p, desired_absorbptions, desired_means, desired_scvs)
    #     @show dist 
    #     sum(absorption_probs(dist) - desired_absorbptions) ≈ 0 && return false
    #     #QQQQ more testing
    # end
    # # @show model_size(out)
    # # data = [rand(out) for _ in 1:10^5]
    # # pi_est = [count((x)->x.a == i, data)/length(data) for i in 1:5]
    # # μ_est = [mean(first.(filter((x)->x.a == i, data))) for i in 1:5] 
    # (absorption_probs(out) ≈ desired_absorbptions') || return false
    

    # Λ₄, λ45, λ54, Λ₅ = 20, 17, 7, 10
    # μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    # T_example = [-Λ₄ λ45; λ54 -Λ₅]
    # T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
    # initial_dist = [0.5,0.5]
    # maph = MAPHDist(initial_dist', T_example, T0_example)
    # @show absorption_probs(maph)
    # @show mean(maph)

    
    # mean(out) ≈ desired_means || return false
    # @show pi_est
    # @show μ_est
    # @show pi_est .* μ_est
    return true
end
