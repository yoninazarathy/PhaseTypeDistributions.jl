"""
QQQQ
"""
function test_example4()
    Λ₄, λ45, λ54, Λ₅ = 5, 2, 7, 10
    μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]

    maph = MAPHDist([0.5,0.5]', T_example, T0_example)

    println("starting simulations")
    data = [rand(maph) for _ in 1:10^2]
    
    p,q = model_size(maph)
    est_maph = MAPHInit(p,q)
    # est_maph = MAPHDist(;model_size(maph)...)
    @assert model_size(est_maph) == model_size(maph)

    #fit! gets an MAPHDist object for two possible reasons:
        #Reason #1 (always) - to know p and q.
        #Reason #2 (sometimes) - to have a starting guess. QQQQ - later give it a flag to say 
    fit!(est_maph,data)

    println("\nfinished simulations")
end
