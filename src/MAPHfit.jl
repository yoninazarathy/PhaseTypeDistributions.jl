function E_step!(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    @time data_length, sorted_stats = compute_sorted_stats(all_obs, maph)

    sum_stats_per_state = map(ss -> sum(ss), sorted_stats)
    total_stats = sum(sum_stats_per_state)

    return data_length, sum_stats_per_state, total_stats
end

function M_step!(data_length::Int, stats::Vector{MAPHSufficientStats}, sum_stats::MAPHSufficientStats, maph::MAPHDist)
    m, n = model_size(maph)
    number_times_exit = [sum_stats.N[i] + sum(sum_stats.M[i,:]) for i ∈ 1:m]
    number_times_exit_per_state = [[stats[k].N[i] + sum(stats[k].M[i, :]) for i ∈ 1:m] for k ∈ 1:n]
    new_α =  sum_stats.B ./data_length
    new_q = [number_times_exit[i] / sum_stats.Z[i] for i ∈ 1:m]
    @show sum_stats.Z
    @show number_times_exit
    new_ρ = [stats[k].B[i] / sum_stats.B[i] for i = 1:m, k = 1:n]
    new_P = [[stats[k].M[i,j] / number_times_exit_per_state[k][i] for i ∈ 1:m, j ∈ 1:m ] for k ∈ 1:n]

    @assert sum(new_α) ≈ 1.0
    @assert sum(new_α .* new_ρ) ≈ 1.0
    @assert all(new_q .≥ 0)  "we have $new_q"
    maph.α = reshape(new_α, (1,m)) 
    update!(maph, new_q, new_ρ, new_P)    

end







function EM_fit!(all_obs::Vector{SingleObservation}, maph::MAPHDist, iterations::Int = 20)
    for _ ∈ 1:iterations
        @time data_length, stats, s_stats = E_step!(all_obs, maph)
        @time M_step!(data_length, stats, s_stats, maph)
    end
end

   