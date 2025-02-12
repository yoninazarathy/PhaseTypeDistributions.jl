module MAPHDistributions

using LinearAlgebra, QuadGK, StatsBase, Distributions, Statistics, Dates, StaticArrays, DocStringExtensions
using JuMP, HiGHS, MathOptInterface
import Base: rand, +, /, -, *
import Distributions: cdf, ccdf, mean, var

include("MAPHDist.jl")
export MAPHDist

include("PHDist.jl")
export get_absorbing_vector
# include("inference/maph_fit.jl")
include("moment_function.jl")
export absorption_probs, kth_moment, mean, var, scv, sub_distribution, sub_pdf, mgf, model_size
include("MAPHStatistics.jl")
export MAPHSufficientStats, +, -, /, *, ContinuousTimeMarkovChain, MAPH_TO_CTMC, SingleObservation, rand,  compute_sufficient_stats, compute_sorted_stats, stats_filter, get_emperical_absorb_prob, very_crude_c_solver
include("MAPHfit.jl")
export EM_fit, E_step, M_step
include("MAPH_initial.jl")
export maph_initialization

const MOI = MathOptInterface


end #end module

