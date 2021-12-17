using LinearAlgebra, QuadGK, StatsBase, Distributions, Statistics
import Base: rand, +, /, -
import Distributions: cdf, ccdf

include("types.jl")
include("util.jl")
include("MAPHDist.jl")
include("inference/maph_fit.jl")
include("inference/maph_compute_sufficient_stats.jl")