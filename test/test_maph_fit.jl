"""
QQQQ
"""
function test_fit_example1(;sim_runs=10^2)
    Λ₄, λ45, λ54, Λ₅ = 5, 2, 7, 10
    μ41, μ42, μ43, μ51, μ52, μ53 = 1, 1, 1, 1, 1, 1 
    T_example = [-Λ₄ λ45; λ54 -Λ₅]
    T0_example = [μ41 μ42 μ43; μ51 μ52 μ53]
    initial_dist = [0.5,0.5]

    maph = MAPHDist(initial_dist', T_example, T0_example)
   
    absorbtion_probs_of_true_dist = absorption_probs(maph)
    mean_of_true_dist = mean(maph)
    scvs_of_true_dist = scv(maph)


    times = zeros(sim_runs)
    absorbing_states = zeros(Int,sim_runs)
    
    @showprogress "Simulating data" for i in 1:sim_runs
        time, state = rand(maph, full_trace = false) 
        times[i] = time
        absorbing_states[i] = state
    end

    @show absorbing_states


    estimated_dist = fit_maph(times, absorbing_states, 2)

    absorbtion_probs_of_estimated_dist = absorption_probs(estimated_dist)
    mean_of_estimated_dist = mean(estimated_dist)
    scvs_of_estimated_dist = scv(estimated_dist)

    @show absorbtion_probs_of_true_dist
    @show absorbtion_probs_of_estimated_dist

    return true
end


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
