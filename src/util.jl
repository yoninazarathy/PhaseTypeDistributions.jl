"""
Compute the SCV of data.
"""
function scv(data::Vector{Float64})
    return (std(data)/mean(data))^2
end 

"""
Merge PH distributions into a MAPH distribution 
"""
function merge_dist(probs::Vector{Float64},dist::Vector{Any})
    merged_α = probs.*[reshape(dist[i].α, model_size(dist[i])) for i = 1:length(dist)]
    # merged_α = (probs ./ sum(probs)) .* [reshape(dist[i].α, model_size(dist[i])) for i = 1:length(dist)]

    α_temp = reduce(vcat,merged_α)'
    merged_T = [dist[i].T for i = 1:length(dist)]
    # merged_T = (probs ./ sum(probs)) .* [dist[i].T for i = 1:length(dist)]
    T_temp = cat(merged_T..., dims = (1,2))
    @show dist
    # merged_T0 = -(probs./sum(probs)) .* [sum(dist[i].T,dims=2) for i = 1:length(dist)]
    merged_T0 = -1.0.*[sum(dist[i].T,dims=2) for i = 1:length(dist)]
    @show merged_T0
    T0_temp = cat(merged_T0...,dims=(1,2))
    display(T0_temp)
    display(T_temp)

    return α_temp, T_temp, T0_temp
end



function vertical_merge_matrix(M::Matrix{Float64},p::Int64)
    m, n = size(M)

    new_M = zeros((p,n))

    for i=1:(p-1)
        new_M[i,:] = M[i,:]
    end

    new_M[p,:] = sum(M[p:m,:],dims=1)

    return new_M
end

function horizontal_merge_matrix(M::Matrix{Float64},p::Int64)
    m, n = size(M)

    new_M = zeros((m,p))

    for i=1:(p-1)
        new_M[:,i] = M[:,i]
    end

    new_M[:,p] = sum(M[:,p:m],dims=2)

    return new_M
end