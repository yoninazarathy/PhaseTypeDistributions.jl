include("MAPH.jl")
include("MAPHStatistics.jl")


function Maximization_step(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    m, n = model_size(maph)

    stats  = compute_expected_stats(all_obs, maph)
    mean(stats).B
    α =  mean(stats).B
    q = map(i -> sum(mean(stats).N[i, :]) / mean(stats).Z[i] , 1:m)
    ρ = map(k -> stats[k].B ./ mean(stats).B, 1:n)
    P = sum(stats).M ./ sum(sum(stats).N, dims = 2)

    ρ
    # maph = MAPH_constructor(α, q, ρ, P)
    


     



end