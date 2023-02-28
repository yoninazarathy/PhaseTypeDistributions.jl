"""

Create an MAPHDist of dimensions pxq where q is the length of `probs`, `means`, and `scvs` and p is specified.

This tries to do a best "moment fit" for the probability of absorbitions, means, and scvs

"""
function MAPHDist(  p::Int, 
                    probs::Vector{Float64}, 
                    means::Vector{Float64}, 
                    scvs::Vector{Float64}; 
                    ω::Float64 = 50.0) #QQQQ - we'll calibrate it sojourn_time

    q = length(probs)
    length(means) != q && error("Dimension mismatch")
    length(scvs) != q && error("Dimension mismatch")
    !(sum(probs) ≈ 1.0) && error("Probabilities need to sum up to 1.0")

    #Sort according to importance based on absorbitions
    π_order = sortperm(probs, rev = true)
    π = probs[π_order]
    means = means[π_order]
    scvs = scvs[π_order]

    #Adjust means and scvs according to ω
    matched_scv(scv::Float64, mean::Float64, ω::Float64) = (scv*(ω*mean+1)^2-1)/(ω*mean)^2
    matched_mean(mean::Float64,ω::Float64) = mean - 1/ω
    scvs = matched_scv.(scvs,means, ω)
    means = matched_mean.(means, ω)

    num_phases = [scv ≥ 1.0 ? 2 : ceil(Int, 1/scv) for scv in scvs]
    required_phases = sum(num_phases)

    if required_phases ≤ p-1
        K = q
    elseif num_phases[1] ≤ p-1 < required_phases
        K = findlast((x)-> sum(x)≤p-1 , accumulate(+,num_phases))
    else #Case of p-1 < num_phases[1]
        K = 0
    end

    K ∉ (1:q) && error("Not enough phases to start") #K = 0 QQQQ - handle this case later
    p2 = sum(num_phases[1:K])
    p1 = p - p2

    used_dist = [scvs[i] ≥ 1 ? hyper_exp_init(means[i], scvs[i]) : hypo_exp_init(means[i], scvs[i]) for i in 1:K]
    α, T, T0 = cat_dist(π[1:K], used_dist, ω, p1, p2, q)

    #permute T0 back to the given order
    T0 = T0[:,sortperm(π_order)]

    return MAPHDist(α, T, T0)
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