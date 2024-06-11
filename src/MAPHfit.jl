function fit!(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    m, n = model_size(maph)

    # E-step
    @info "start fitting E-step, computing sufficient stats"
    stats  = compute_expected_stats(all_obs, maph)

    #QQQQ fix mean issue for absorbing probs...

    #M -step
    s_stats = sum(stats)
    @info "M-step, start estimating the paramters"
    α =  mean(stats).B
    q = map(i -> (sum(s_stats.N[i, :]) + sum(s_stats.M[i, :]) ) / s_stats.Z[i] , 1:m)

    ρ = replace(reduce(hcat, map(k -> stats[k].B ./ s_stats.B, 1:n)), NaN => 0)
    @assert sum(α) ≈ 1.0
    @assert sum(α .* ρ) ≈ 1.0
    @assert all(q .≥ 0)  "we have $q"
 
    P = map(k -> stats[k].M ./ sum(stats[k].N), 1:length(stats) )

    maph.α = reshape(α, (1,m)) 
    update!(maph, q, ρ, P)    
end