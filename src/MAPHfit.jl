function E_step!(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    sorted_stats = compute_sorted_stats(all_obs, maph)

    sum_stats_per_state = map(ss -> sum(ss), sorted_stats)
    total_stats = sum(sum_stats_per_state)

    # @show keys(absorb_prob)
    return sum_stats_per_state, total_stats
end

function M_step!(num_obs::Int, stats::Vector{MAPHSufficientStats}, sum_stats::MAPHSufficientStats, maph::MAPHDist)
    m, n = model_size(maph)

    new_α =  sum_stats.B ./ num_obs

    new_q = map(i -> (sum_stats.N[i] + sum(sum_stats.M[i, :]) ) , 1:m) ./ sum_stats.Z
    new_ρ = replace(reduce(hcat, map(k -> sum(stats[k]).B ./ sum_stats.B, 1:n)), NaN => 0)
    @assert sum(new_α) ≈ 1.0
    @assert sum(new_α .* new_ρ) ≈ 1.0
    @assert all(new_q .≥ 0)  "we have $new_q"
    new_P = map(k -> sum(stats[k]).M ./ sum(stats[k]).N, 1:n)

    maph.α = reshape(new_α, (1,m)) 
    update!(maph, new_q, new_ρ, new_P)    
end







function EM_fit!(all_obs::Vector{SingleObservation}, maph::MAPHDist, iterations::Int = 3)
    #E-step
    num_obs = length(all_obs)
    for _ ∈ 1:iterations
        stats, s_stats = E_step!(all_obs, maph)
        M_step!(num_obs, stats, s_stats, maph)
    end
end

   