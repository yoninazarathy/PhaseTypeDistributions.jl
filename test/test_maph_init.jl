<<<<<<< HEAD
function init_test()
    p=20
    probs = [0.3,0.5,0.1,0.08,0.02]
    means = [2.0,3.1,2.3,4.5,0.2]
    scvs = [1.1,2.3,2.0,0.33,1.5]
    out = MAPHDist(p,probs,means,scvs)
    data = [rand(out) for _ in 1:10^5]
    π_est = [count((x)->x.a == i, data)/length(data) for i in 1:5]
    μ_est = [mean(first.(filter((x)->x.a == i, data))) for i in 1:5] 

    norm(probs.-π_est)<0.1 || return false
    norm(means.-π_est .* μ_est)<0.1 || return false

=======
"""
QQQQ
"""
function test_maph_init()
    desired_absorbptions = [0.3,0.5,0.1,0.08,0.02]
    desired_means = [2.0,3.1,2.3,4.5,0.2]
    out = MAPHDist(20,desired_absorbptions,desired_means,[1.1,2.3,2.0,0.33,1.5])
    # @show model_size(out)
    # data = [rand(out) for _ in 1:10^5]
    # pi_est = [count((x)->x.a == i, data)/length(data) for i in 1:5]
    # μ_est = [mean(first.(filter((x)->x.a == i, data))) for i in 1:5] 

    absorption_probs(out) ≈ desired_absorbptions' || return false
    
    # QQQQ  - continue here...

    
    # mean(out) ≈ desired_means || return false
    # @show pi_est
    # @show μ_est
    # @show pi_est .* μ_est
>>>>>>> d71935519d27b6b8e7403df7fab131c832b45a11
    return true
end
