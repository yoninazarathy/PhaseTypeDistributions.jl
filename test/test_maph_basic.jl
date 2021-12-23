function test_maph_moments_and_rand()
    #QQQQ create an MAPH Distributions
    Λ₄, λ45, λ54, Λ₅ = 5, 2, 7, 10
    μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]

    maph = MAPHDist([0.3,0.7]', T_example, T0_example)

    data = [rand(maph) for _ in 1:10^5]


    μ_est = [mean(first.(filter((x)->x.a == i, data))) for i in 1:3] 

    σ_est = [var(first.(filter((x)->x.a == i, data))) for i in 1:3]

    scv_est = [scv(first.(filter((x)->x.a == i, data))) for i in 1:3]

    π_est = [count((x)->x.a == i, data)/length(data) for i in 1:3]


    first_moment = mean(maph)

    variance = var(maph)

    norm(first_moment.-μ_est)<0.01 || return false

    norm(variance.-σ_est)<0.01 || return false

    norm(scv(maph).-scv_est)<0.1 || return false

    norm(absorption_probs(maph).-π_est)<0.01 || return false

    #QQQQ generate 10^6 rands... and see that estimated mean and scv match the mean and scv
    return true
end