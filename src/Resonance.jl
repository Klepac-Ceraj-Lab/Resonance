module Resonance

# Data loading     
export datafiles,
       scratchfiles,
       analysisfiles,
       outputfiles,
       figurefiles,
       tablefiles,
       inputfiles,
       Metadata,
       ReadCounts,
       TaxonomicProfiles,
       UnirefProfiles,
       KOProfiles,
       ECProfiles,
       Neuroimaging,
       BrainVolume,
       hemisphere,
       hashemisphere

# Plotting
export loadings,
       varexplained,
       mdsaxis,
       permanovas,
       plot_permanovas,
       plot_permanovas!,
       mantel,
       plot_mantel,
       plot_mantel!,
       plot_pcoa!,
       commonname,
       comm2wide,       
       runlms

# Machine Learning
export  # 0. structs/types
        Prediction,
        Classification,
        Regression,
        ProbeData,
        # 1. preprocessing functions
        filter_prevalences,
        # 2. training functions
        partitionvec,
        probe_prod_randomforest,
        # postprocessing functions
        CustomRangeNormalizer,
        compute_custom_scale

# Percentiles
export  AgeBracketPercentiles,
        GrowthCurve,
        compute_age_bracket,
        compute_growth_curve,
        get_cogscore_percentile,
        get_brain_percentile,
        plot_multiple_growthcurves!,
        plot_all_results!

using Reexport
using ReTest

using Airtable
using Arrow
using AlgebraOfGraphics
using CairoMakie
using CategoricalArrays
using CodecZlib
using ColorSchemes
using Combinatorics
using Dates
using Dictionaries
using GLM
using LinearAlgebra
using MultipleTesting
using MultivariateStats
using PERMANOVA
using ProgressLogging
using Random
using Setup # dev package at `./Setup`
using SparseArrays
using Statistics
using Tables
using ThreadsX
using MLJ
using DecisionTree
using CubicSplines

@reexport using BiobakeryUtils
@reexport using DataFrames
@reexport using Chain
@reexport using CSV
@reexport using XLSX

const transform = DataFrames.transform

include("airtable.jl")
include("rawdata.jl")
include("microbiomejl_interface.jl")
include("files.jl")
include("genefamilies.jl")
include("plotting.jl")
include("omnibus.jl")
include("prediction.jl")
include("data_loading.jl")
include("percentiles.jl")
include("lms.jl")

end
