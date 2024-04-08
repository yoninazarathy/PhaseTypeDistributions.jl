include("MAPH.jl")
include("MAPHStatistics.jl")


function maph_initialization(all_obs::Vector{SingleObservation}, m::Int, ω::Real = 100.0)
    @assert m ≥ 1 "cannot have empty phase"

    prob_per_state = get_emperical_absorb_prob(all_obs)
    mean_per_state, scv_per_state= get_emperical_statistics(all_obs, ω)
    prob_per_state = sort(prob_per_state; byvalue = true, rev = true)

    phases_required_per_state = Dict()

    for key ∈ keys(prob_per_state)
        phases_required_per_state[key] = scv_per_state[key] > 1.0 ? 2.0 : ceil(1/scv_per_state[key])        
    end

    m_required = sum(collect(values(phases_required_per_state)))
    @info "to use the statistics from each absorbing states we need at least $(m_required + 1)"
    num_phases = 0
    num_phase_type_distributions = 0
    
    @info "each of the absorbing state requires num of phases as $(phases_required_per_state)"

    if m_required ≤ m - 1
        @info "we have enough phases to start the maph_initialization, states used are absorbing state $(Int.(collect(keys(prob_per_state))))"
        num_phases = m
        num_phase_type_distributions = length(keys(prob_per_state))
    elseif m-1 < phases_required_per_state[collect(keys(prob_per_state))[1]]
        @info "will initialize with an exponential distribution"
        num_phases = 0
        maph = (α = ones(1, m) ./ m ; T = diagm(-ω .* ones(m)) ; D = ω .* ones(m, length(prob_per_state)) ./ (length(prob_per_state)) ; MAPH_constructor(α, T, D))
        return maph
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
        num_phase_type_distributions = k-1
    end
    @info "we now start to initialze with one exponential distribution with parameter ω and PH distributions initialized by the absorbing states $(Int.(collect(keys(prob_per_state))[1:num_phase_type_distributions]))"
    

    ph_distributions = Dict()
    for key ∈ collect(keys(prob_per_state))[1:num_phase_type_distributions]
        if scv_per_state[key] > 1.0
            @info "we now construct a hyper exponential distribution for absorbing state $(key), with mean $(mean_per_state[key]) and scv $(scv_per_state[key])"
            ph_distributions[key] = hyper_exp_dist(mean_per_state[key], scv_per_state[key])
        else
            @info "we now construct a hypo exponential distribution for absorbing state $(key), with mean $(mean_per_state[key]) and scv $(scv_per_state[key])"
            ph_distributions[key] = hypo_exp_dist(mean_per_state[key], scv_per_state[key])
        end
    end
    expo_phases = m - Int(num_phases)

    α = hcat(ones(1, expo_phases) ./ (expo_phases), zeros(1, Int(num_phases)))

    T_expo = diagm(-ω .* ones(expo_phases))
    D_expo = ω .* ones(expo_phases, length(prob_per_state)) ./ (length(prob_per_state))
    @show T_expo
    @show D_expo

    T = T_expo
    D = D_expo
    for key ∈ reverse(collect(keys(mean_per_state)))
        if key ∈ keys(ph_distributions)
           T = cat(T, ph_distributions[key].T, dims = (1,2))
           D =  sum(ph_distributions[4].T, dims = 1)'
           @show D
        end
    end

    @show  sum(ph_distributions[4].T, dims = 1)
    @show D
    # @show reverse(collect(keys(mean_per_state)))


    # @show cat(T_expo, collect(values(ph_distributions))[1].T, dims = (1,2))
        


end