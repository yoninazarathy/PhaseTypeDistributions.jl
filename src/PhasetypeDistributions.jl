using Revise

module PhaseTypeDistributions

using LinearAlgebra, QuadGK, StatsBase, Distributions, Statistics, ProgressMeter, Dates
import Base: rand, +, /, -
import Distributions: cdf, ccdf, mean, var

include("types.jl")
include("util.jl")
include("MAPHDist.jl")
include("inference/maph_fit.jl")
include("structured_ph.jl")
include("moments.jl")
include("inference/maph_compute_sufficient_stats.jl")

export  MAPHDist, 
        mean, 
        model_size, 
        sufficient_stat_from_trajectory, 
        observation_from_full_traj, 
        sufficient_stats, 
        absorb_filter_data, 
        time_filter_data, 
        maximum_likelihood_estimate_second_parameter #QQQQ add all other exports


end #end module