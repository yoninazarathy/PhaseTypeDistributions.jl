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


function cat_dist(probs::Vector{Float64},dist::Vector{Any},ω::Float64,p1::Int64,p2::Int64,q::Int64)
    merged_T = [dist[i].T for i = 1:length(dist)]
    TS2S2 = cat(merged_T..., dims = (1,2))

    merged_α = [ω*probs[i]/sum(probs).*dist[i].α for i = 1:length(dist)]

    merged_T0 = [-1.0.*sum(dist[i].T,dims = 2) for i =1:length(dist)]

    @show merged_T0

    S2T0 = cat(merged_T0...,dims = (1,2))

    v = cat(merged_α..., dims = (2,2))
    
    TS1S2 = repeat(v,p1)

    TS1S2 = cat(merged_α..., dims = (2,2))


    Imatrix = Matrix(-ω.*I(Int(p1)))

    M = cat(Imatrix,zeros(p2,p1),dims = (1,1))

    N = cat(TS1S2,TS2S2,dims = (1,1))

    T = cat(M,N, dims = (2,2))

    T0 = cat(zeros(p1,length(dist)),S2T0,dims = (1,1))

    if length(dist) < q
        @show 1
        T0 = cat(T0,zeros(p1+p2,q-length(dist)),dims = (2,2))
    end


    α = zeros(1,p1+p2)
    
    α[1:p1] .= 1/p1

    α = vec(α)'


    return  α,T,T0

end





function absorption_vector_create(T1::Matrix{Float64},prob::Vector{Float64})
    p = length(T)
    q = length(probs)
    T0 = zeros((p,q))

    for i = 1:p
        T0[i,:] = T1[i,1].*probs
    end

end









