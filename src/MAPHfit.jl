include("MAPH.jl")
include("MAPHStatistics.jl")


function Maximization_step(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    m, n = model_size(maph)

    stats  = compute_expected_stats(all_obs, maph)
    mean(stats).B
    α =  mean(stats).B
    q = map(i -> sum(sum(stats).N[i, :]) / sum(stats).Z[i] , 1:m)
    ρ = reduce(hcat, map(k -> stats[k].B ./ sum(stats).B, 1:n))
    P = sum(stats).M ./ sum(sum(stats).N, dims = 2)
    
    @show α
    @show q
    @show ρ
    @show P

    


    

end