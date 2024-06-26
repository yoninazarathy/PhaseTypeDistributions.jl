module MAPHDistributions

using LinearAlgebra, QuadGK, StatsBase, Distributions, Statistics, ProgressMeter, Dates, StaticArrays, DocStringExtensions
import Base: rand, +, /, -, *
import Distributions: cdf, ccdf, mean, var

include("MAPH.jl")
export MAPHDist, MAPH_constructor, update!
include("PHDist.jl")
export get_absorbing_vector
# include("inference/maph_fit.jl")
include("moment_function.jl")
export absorption_probs, kth_moment, mean, var, scv, sub_distribution, sub_pdf, mgf, model_size
include("MAPHStatistics.jl")
export MAPHSufficientStats, +, -, /, *, ContinuousTimeMarkovChain, MAPH_TO_CTMC, SingleObservation, rand,  compute_sufficient_stats, compute_expected_stats, stats_filter, get_emperical_absorb_prob
include("MAPHfit.jl")
export EM_fit!
include("MAPH_initial.jl")
export maph_initialization


# include("types.jl")
# include("util.jl")

#include("inference/maph_compute_sufficient_stats.jl")

# export  MAPHDist, 
#         mean,scv,var,
#         model_size, 
#         sufficient_stat_from_trajectory, 
#         observation_from_full_traj, 
#         sufficient_stats, 
#         absorb_filter_data, 
#         time_filter_data, 
#         maximum_likelihood_estimate_second_parameter,
#         compute_descriptive_stats
#         randn#QQQQ add all other exports


end #end module

