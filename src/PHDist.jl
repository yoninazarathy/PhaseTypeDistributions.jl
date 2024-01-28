"""
$(TYPEDFIELDS)
The mutable struct which representing the MAPH distribution in different parameterizations
"""
mutable struct PHDist
    "m × 1 Initial probability distribution across the phases. Elements should sum to 1"
    α::Union{Matrix{Float64}, Vector{Float64}}

    "m x m transition rate matrix for the MAPH states"    
    T::Matrix{Float64}
    PHDist(α::Union{Matrix{Float64}, Vector{Float64}}, T::Matrix{Float64}) = new(α, T)
end

function ph_constructor(α::Union{Matrix{Float64}, Vector{Float64}}, T::Matrix{Float64})
    @assert size(T,1) == size(T, 2) "T must be a square matrix"
    @assert size(α, 2) == size(T, 1) "The length of α must be equal to the number of rows in T"
    return PHDist(α, T)
end








