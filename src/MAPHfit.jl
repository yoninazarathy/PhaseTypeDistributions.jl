function E_step!(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    stats = compute_expected_stats(all_obs, maph)
    s_stats = sum(stats)
    return stats, s_stats
end

function M_step!(stats::Vector{MAPHSufficientStats}, s_stats::MAPHSufficientStats, maph::MAPHDist)
    m, n = model_size(maph)
    new_α =  mean(stats).B
    new_q = map(i -> (sum(s_stats.N[i, :]) + sum(s_stats.M[i, :]) ) / s_stats.Z[i] , 1:m)
    new_ρ = replace(reduce(hcat, map(k -> stats[k].B ./ s_stats.B, 1:n)), NaN => 0)
    @assert sum(new_α) ≈ 1.0
    @assert sum(new_α .* new_ρ) ≈ 1.0
    @assert all(new_q .≥ 0)  "we have $new_q"
    new_P = map(k -> stats[k].M ./ sum(stats[k].N), 1:length(stats) )

    maph.α = reshape(new_α, (1,m)) 
    update!(maph, new_q, new_ρ, new_P)    
end







function EM_fit!(all_obs::Vector{SingleObservation}, maph::MAPHDist, iterations::Int = 100)
    #E-step
    for _ ∈ 1:iterations
        stats, s_stats = E_step!(all_obs, maph)
        M_step!(stats,s_stats, maph)
    end
end

   