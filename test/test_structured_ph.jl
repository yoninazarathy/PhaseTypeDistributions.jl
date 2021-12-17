function hypoexp_test()
    mean_desired = 5.7
    scv_desired = 0.23
    d = hypo_exp_init(mean_desired, scv_desired)
    (mean_desired ≈ mean(d)) || return false 
    (scv_desired ≈ scv(d)) || return false #QQQQ
    return true
end

function hyperexp_test()
    mean_desired = 1.0
    scv_desired = 1.0
    d = hyper_exp_init(mean_desired, scv_desired)
    @show mean(d)
    # (mean_desired ≈ mean(d)) || return false 
    # (scv_desired ≈ scv(d)) || return false #QQQQ
    return true
end