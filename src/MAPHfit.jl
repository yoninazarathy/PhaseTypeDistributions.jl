include("MAPH.jl")
include("MAPHStatistics.jl")


function Maximization_step!(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    m, n = model_size(maph)

    stats  = compute_expected_stats(all_obs, maph)
    @show stats

    α =  mean(stats).B
    @show [sum(stats).N[i, :] for i in 1:m]
    @show [sum(stats).Z[i] for i in 1:m]
    q = map(i -> sum(sum(stats).N[i, :]) / sum(stats).Z[i] , 1:m)
    ρ = replace(reduce(hcat, map(k -> stats[k].B ./ sum(stats).B, 1:n)), NaN => 0)
    @show sum(α .* ρ)
    @info "hello"
    @assert sum(α) ≈ 1.0
    @assert sum(α .* ρ) ≈ 1.0
    @show sum(α .* ρ)
    # if !all(q .> 0)
    @show q
    # end
    @assert all(q .> 0)

    P = map(k -> stats[k].M ./ sum(stats[k].N), 1:length(stats) )

    maph.α = reshape(α, (1,m)) 
    update!(maph, q, ρ, P)    
end