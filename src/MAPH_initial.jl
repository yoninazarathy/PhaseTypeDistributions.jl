include("MAPH.jl")
include("MAPHStatistics.jl")


function maph_initialization(all_obs::Vector{SingleObservation}, m::Int, n::Int)

    @show sortperm(get_emperical_absorb_prob(all_obs), rev = false)
    @show get_emperical_absorb_prob(all_obs)

    
end
