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
<<<<<<< HEAD
    merged_α = (probs./sum(probs)).*[reshape(dist[i].α,length(dist[i].α)) for i = 1:length(dist)]
    α_temp = reduce(vcat,merged_α)'
    merged_T = (probs./sum(probs)).*[dist[i].T for i = 1:length(dist)]

    T_temp = cat(merged_T...,dims = (1,2))

    merged_T0 = -(probs./sum(probs)).*[sum(dist[i].T,dims=2) for i = 1:length(dist)]

=======
    merged_α = (probs ./ sum(probs)) .* [reshape(dist[i].α, model_size(dist[i])) for i = 1:length(dist)]
    α_temp = reduce(vcat,merged_α)'
    merged_T = (probs ./ sum(probs)) .* [dist[i].T for i = 1:length(dist)]
    T_temp = cat(merged_T..., dims = (1,2))
    merged_T0 = -(probs./sum(probs)) .* [sum(dist[i].T,dims=2) for i = 1:length(dist)]
>>>>>>> d71935519d27b6b8e7403df7fab131c832b45a11
    T0_temp = cat(merged_T0...,dims=(1,2))

    return α_temp, T_temp, T0_temp
end

