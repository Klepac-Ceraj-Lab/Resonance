module Resonance

# Data loading     
export datafiles,
       scratchfiles,
       analysisfiles,
       outputfiles,
       figurefiles,
       Metadata,
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
       commonname

# Machine Learning
export  # 0. structs/types
        Prediction,
        Classification,
        Regression,
        UnivariateRandomForestClassifier,
        UnivariateRandomForestRegressor,
        UnivariatePredictorEnsemble,
        # 1. preprocessing functions
        dropmissing,
        dropnan,
        myxor,
        filter_age_bracket,
        build_metadata_prediction_df,
        prepare_future_prediction_df,
        meanclass,
        compute_tietjenmoore,
        test_tietjenmoore,
        univariate_tietjenmoore,
        try_outliers,
        # 2. training functions
        train_randomforest,
        # postprocessing functions
        predict_proba,
        predict,
        report_merits,
        singlemodel_singlesplit_importance,
        singlemodel_allsplits_importances,
        singlemodel_summary_importances,
        multimodel_individual_summaryimportances,
        multimodel_aggregate_summaryimportances,
        singlemodel_binarytopn_importances,
        multimodel_individual_binarytopns,
        multimodel_aggregate_binarytopns,
        descript_inputs,
        singlemodel_summary_prevalences,
        singlemodel_summary_abundances,
        singlemodel_summary_sdevs,
        multimodel_individual_prevalences,
        multimodel_individual_abundances,
        biserial_correlation,
        feature_split_correlation_analysis,
        singlemodel_singlesplit_correlations,
        singlemodel_allsplits_correlations,
        singlemodel_allsplits_correlationsummary,
        multimodel_individual_correlations,
        # plotting functions
        build_confusion_matrix,
        average_confusion_matrix,
        confmatrix2barplot,
        singlemodel_merit_barplot!,
        singlemodel_merit_scatterplot!,
        singlemodel_avgimportance_barplot!,
        multimodel_avgimportance_barplot!,
        singlemodel_logistic_regression!,
        multimodel_logistic_regression!

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

using Arrow
using AlgebraOfGraphics
using CairoMakie
using CategoricalArrays
using ColorSchemes
using Combinatorics
using Dates
using LinearAlgebra
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
using CubicSplines

@reexport using BiobakeryUtils
@reexport using DataFrames
@reexport using Chain
@reexport using CSV
@reexport using XLSX

const transform = DataFrames.transform

include("files.jl")
include("wrangle.jl")
include("genefamilies.jl")
include("plotting.jl")
include("omnibus.jl")
include("prediction.jl")
include("data_loading.jl")
include("percentiles.jl")

end
