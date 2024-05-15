include("MAPH.jl")
include("MAPHStatistics.jl")


function fit!(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    m, n = model_size(maph)

    # E-step
    @info "start fitting E-step, computing sufficient stats"
    stats  = compute_expected_stats(all_obs, maph)

    #find the error! .. the first states in initilization is unconnected to the other states.. we need to fix

    #M -step
    @info "M-step, start estimating the paramters"
    α =  mean(stats).B
    q = map(i -> (sum(sum(stats).N[i, :]) + sum(sum(stats).M[i, :]) ) / sum(stats).Z[i] , 1:m)

    ρ = replace(reduce(hcat, map(k -> stats[k].B ./ sum(stats).B, 1:n)), NaN => 0)
    @assert sum(α) ≈ 1.0
    @assert sum(α .* ρ) ≈ 1.0
    @assert all(q .> 0)
 
    P = map(k -> stats[k].M ./ sum(stats[k].N), 1:length(stats) )

    maph.α = reshape(α, (1,m)) 
    update!(maph, q, ρ, P)    
end