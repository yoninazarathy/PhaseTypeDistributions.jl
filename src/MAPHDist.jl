"""

Create an MAPHDist of dimensions pxq where q is the length of `probs`, `means`, and `scvs` and p is specified.

This tries to do a best "moment fit" for the probability of absorbitions, means, and scvs

"""
function MAPHDist(p::Int, probs::Vector{Float64}, means::Vector{Float64}, scvs::Vector{Float64})
    q = length(probs)
    length(means) != q && error("Dimension mismatch")
    length(scvs) != q && error("Dimension mismatch")
   
    πhat_order = sortperm(probs)
    sorted_scvs = scvs[πhat_order]
    sorted_means = means[πhat_order]

    num_phases = [sorted_scvs[i] ≥ 1 ? 2 : ceil(1/sorted_scvs[i]) for i in 1:q]

    required_phases = sum(num_phases)
    K=0
    K = findfirst((x)->x ≤ p, [sum(num_phases[k:end]) for k=1:q])
    K ∉ (1:q) && error("Not enough phases to start")


    dist = []
    for k = 1:q
        if sorted_scvs[k]≥1
            push!(dist,hyper_exp_init(sorted_means[k],sorted_scvs[k]))
        end
        if sorted_scvs[k]<1
            push!(dist,hypo_exp_init(sorted_means[k],sorted_scvs[k]))
        end
    end

    α = zeros(p)'
    T = ones(p,p).*eps()
    T0 = ones(p,q).*eps() 


    reverse_order = sortperm(πhat_order)
    
    reversed_dist = dist[reverse_order]


    if p ≥ required_phases
        α_temp, T_temp, T0_temp = merge_dist(probs,reversed_dist)
        m,n = model_size(MAPHDist(α_temp, T_temp, T0_temp))
        α[1:m] = α_temp
        T[1:m,1:m] = T_temp
        T0[1:m,1:n] = T0_temp
    else
        idx_drop = πhat_order[1:K-1]
        probs = probs[eachindex(probs) .∉ Ref(idx_drop)]
        reversed_dist = reversed_dist[eachindex(reversed_dist) .∉ Ref(idx_drop)]
        α_temp, T_temp, T0_temp = merge_dist(probs,reversed_dist)
        m,n = model_size(MAPHDist(α_temp, T_temp, T0_temp))
        α[1:m] =α_temp
        T[1:m,1:m] = T_temp
        T0[1:m,1:n]=T0_temp
    end

    maph = MAPHDist(α,T,T0)

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
