module Resonance

# Data loading     
export datafiles,
       scratchfiles,
       analysisfiles,
       outputfiles,
       figurefiles,
       modelfiles,
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
       format_species_labels,
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
        compute_custom_scale,
        calculate_fitness,
        weighted_hpimportances,
        compute_joined_importances,
        plot_comparativedemo_importance_barplots!,
        attribute_colors,
        plot_comparative_lmvsrf_scatterplots!,
        plot_taxon_deepdive!,
        singlemodel_importances_suppltable,
        plot_importances_pareto!

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
using Random
using Setup # dev package at `./Setup`
using SparseArrays
using Statistics
using Tables
using ThreadsX
using MLJ
using DecisionTree

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
include("lms.jl")
include("mbx.jl")

end
