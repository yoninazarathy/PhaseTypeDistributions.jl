include("MAPH.jl")
include("MAPHStatistics.jl")


function Maximization_step!(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    m, n = model_size(maph)

    stats  = compute_expected_stats(all_obs, maph)

    α =  mean(stats).B
    q = map(i -> sum(sum(stats).N[i, :]) / sum(stats).Z[i] , 1:m)
    ρ = reduce(hcat, map(k -> stats[k].B ./ sum(stats).B, 1:n))

    @assert sum(α) ≈ 1.0
    @assert sum(α .* ρ) ≈ 1.0
    @assert all(q .> 0)

    P = map(k -> stats[k].M ./ sum(stats[k].N), 1:length(stats) )
    maph.α = reshape(α, (1,m)) 
    
    # maph.α = α
    update!(maph, q, ρ, P)    
end