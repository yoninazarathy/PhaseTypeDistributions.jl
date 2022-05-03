function maph_moments_and_rand_test(;N=1000)
    #QQQQ create an MAPH Distributions
    Λ₄, λ45, λ54, Λ₅ = 7.0, 2.0, 7.0, 14.0
    μ41, μ42, μ43, μ51, μ52, μ53 = 2.0, 1.0, 2.0, 1.0, 4.0, 2.0 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]

    maph = MAPHDist([0.3,0.7]', T_example, T0_example)

    data = [rand(maph) for _ in 1:N]
    π_est = [count((x)->x.a == i, data)/length(data) for i in 1:3]

    μ_est = [mean(first.(filter((x)->x.a == i, data))) for i in 1:3]

    

    
    σ2_est = [var(first.(filter((x)->x.a == i, data))) for i in 1:3]

    @show mean(maph), μ_est
    @show var(maph),σ2_est


    scv_est = [scv(first.(filter((x)->x.a == i, data))) for i in 1:3]

    @show scv(maph), scv_est
    


    # first_moment = mean(maph)
    # variance = var(maph)

    # (norm(first_moment .- μ_est) < 10/sqrt(N)) || (return false)
    # (norm(variance .- σ2_est) < 10/sqrt(N)) || (return false)
    # (norm(scv(maph) .- scv_est) < 20/sqrt(N)) || (return false)
    # (norm(absorption_probs(maph) .- π_est) < 5/sqrt(N)) || (return false)

    return true
end