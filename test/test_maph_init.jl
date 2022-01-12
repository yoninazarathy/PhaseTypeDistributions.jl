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

    return true
end
