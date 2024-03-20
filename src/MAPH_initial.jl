include("MAPH.jl")
include("MAPHStatistics.jl")


function maph_initialization(all_obs::Vector{SingleObservation}, m::Int, ω::Real = 100.0)
    @assert m ≥ 1 "cannot have empty phase"

    prob_per_state = get_emperical_absorb_prob(all_obs)
    mean_per_state, scv_per_state= get_emperical_statistics(all_obs, ω)
    prob_per_state = sort(prob_per_state; byvalue = true, rev = true)

    phases_required_per_state = Dict()

    for key ∈ keys(scv_per_state)
        phases_required_per_state[key] = scv_per_state[key] > 1.0 ? 2.0 : ceil(1/scv_per_state[key])        
    end

    m_required = sum(collect(values(phases_required_per_state)))
    num_phases = 0


    if m_required ≤ m -1
        @info "we have enough phases to start the maph_initialization"
        num_phases = m
    elseif m-1 < phases_required_per_state[collect(keys(prob_per_state))[1]]
        @info "will initialize with an exponential distribution"
        num_phases = 0
    elseif (phases_required_per_state[collect(keys(prob_per_state))[1]] ≤ m-1) &&  (m-1 < m_required)
        k = 1
        phases = 0
        while phases ≤ m-1
            phases_required = phases_required_per_state[collect(keys(prob_per_state))[k]]
            phases += phases_required
            if phases > m-1
                break
            end
            k += 1 
        end
        @info "will only use part of the absorbing states to initialize, states used are absorbing state $(Int.(collect(keys(prob_per_state))[1:(k-1)]))"
        num_phases = sum(map(i -> phases_required_per_state[collect(keys(prob_per_state))[i]], 1:(k-1)))
    end
    



    
end
