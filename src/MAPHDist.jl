"""

Create an MAPHDist of dimensions pxq where q is the length of `probs`, `means`, and `scvs` and p is specified.

This tries to do a best "moment fit" for the probability of absorbitions, means, and scvs

"""
function MAPHDist(p::Int, probs::Vector{Float64}, means::Vector{Float64}, scvs::Vector{Float64}; ω::Float64 = 50.0)
    @show "Constructor based on moments"

    q = length(probs)
    length(means) != q && error("Dimension mismatch")
    length(scvs) != q && error("Dimension mismatch")
   
    π_order = sortperm(probs,rev = true)


    π = probs[π_order]

    @show π
    
    sorted_means = means[π_order]
    sorted_scvs = scvs[π_order]

    @show sorted_means


    sorted_means = matched_mean.(sorted_means, ω)
    sorted_scvs = matched_scv.(sorted_scvs,sorted_means, ω)

    num_phases = [sorted_scvs[i] ≥ 1.0 ? 2 : ceil(Int, 1/sorted_scvs[i]) for i in 1:q]
    required_phases = sum(num_phases)

    if required_phases ≤ p-1
        K = q
    elseif num_phases[1] ≤ p-1 < required_phases
        K = findlast((x)-> sum(x)≤p-1 ,  [sum(num_phases[1:k]) for k=1:q])
    elseif p-1 < num_phases[1] 
        K = 0
    end


    K ∉ (1:q) && error("Not enough phases to start")

    p2 = sum(num_phases[1:K])
    p1 = p - p2

    used_absorb = π_order[1:K]

    @show K, num_phases, p2, p1, used_absorb, π_order, π[used_absorb]
    
    used_dist = []
    for ua in used_absorb
        if sorted_scvs[ua] ≥ 1
            push!(used_dist, hyper_exp_init(sorted_means[ua], sorted_scvs[ua]))
        else
            push!(used_dist, hypo_exp_init(sorted_means[ua], sorted_scvs[ua]))
        end
    end


    α_test, T_test, T0_test = cat_dist(π[1:K],used_dist,ω,Int(p1),Int(p2),q)



    maph = MAPHDist(α_test,T_test,T0_test)

    display(maph.T)
    display(maph.T0)
    display(maph.α)

    return maph

end




function q_matrix(d::MAPHDist)::Matrix{Float64}
    p, q = model_size(d)
    return [zeros(q,p) zeros(q,q) ; d.T0 d.T]
end

function p_matrix(d::MAPHDist)::Matrix{Float64}
    p, q = model_size(d)
    PT = I-inv(Diagonal(d.T))*d.T
    PT0 = -inv(Diagonal(d.T))*d.T0
    return [I zeros(q,p); PT0 PT]
end

function rand(d::MAPHDist; full_trace = false)
    p, q = model_size(d)
    transient_states = (q+1):(q+p)
    all_states = 1:(q+p)
    Pjump = p_matrix(d)
    Λ = vcat(zeros(q), -diag(d.T))

    if full_trace
        states = Int[]
        sojourn_times = Float64[]
    end

    state = sample(transient_states,weights(d.α)) #initial state
    t = 0.0
    while state ∈ transient_states
        sojourn_time = rand(Exponential(1/Λ[state]))
        if full_trace
            push!(states,state)
            push!(sojourn_times,sojourn_time)
        end
        t += sojourn_time
        state = sample(all_states,weights(Pjump[state,:]))
    end

    if full_trace
        push!(states,state)
        return (sojourn_times, states)
    else
        return (y = t, a = state)
    end
end


"""
This function permutes the output (abosrbing) states of an MAPH Distribtuion.
"""
# function permute(d::MAPHDist,perm::Vector{Int})
#     return MAPHDist(d.α, d.T, d.T0[:,perm])
# end