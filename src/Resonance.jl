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
       MetabolicProfiles,
       Neuroimaging,
       BrainVolume,
       hemisphere,
       hashemisphere

       # Data wrangling
export airtable_metadata,
       brain_ingest,
       findprevstool,
       count_set,
       upset_dots!,
       load_metabolites, # should replace with load(MetabolicProfiles())
       pull_row,
       countmap,
       codebreastfeeding!,
       stp_overlap,
       comm_overlap

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
export  # structs/types
        Prediction,
        Classification,
        Regression,
        ResonancePredictor,
        ResonanceUnivariatePredictor,
        ResonanceMultivariatePredictor,
        UnivariateRandomForestClassifier,
        UnivariateRandomForestRegressor,
        UnivariatePredictorEnsemble,
        # preprocessing functions
        dropmissing,
        dropnan,
        nonna_mean,
        nonmissing_mean,
        nonna_nonmissing_mean,
        myxor,
        tryparsecol,
        check_longdata_metaduplicates!,
        unstack_techreplicates,
        build_metadata_prediction_df,
        univariate_tietjenmoore,
        test_tietjenmoore,
        try_outliers,
        # training functions
        train_randomforest,
        # postprocessing functions
        report_merits,
        get_singlemodel_singlesplit_importance,
        get_singlemodel_allsplits_importances,
        get_singlemodel_summary_importances,
        get_singlemodel_binarytopn_importances,
        get_multimodel_individual_summaryimportances,
        get_multimodel_individual_binarytopns,
        get_multimodel_aggregate_summaryimportances,
        get_multimodel_aggregate_binarytopns,
        build_confusion_matrix,
        average_confusion_matrix,
        confmatrix2barplot,
        regression_bestprediction,
        # plotting functions
        singlemodel_merit_barplot!,
        singlemodel_merit_scatterplot!,
        singlemodel_avgimportance_barplot!,
        multimodel_avgimportance_barplot!

using Reexport
using ReTest

using Arrow
using Airtable
using AlgebraOfGraphics
using CairoMakie
using CategoricalArrays
using CodecZlib
using ColorSchemes
using Combinatorics
using Dates
using Dictionaries
using LinearAlgebra
using MultivariateStats
using PERMANOVA
using Random
using SparseArrays
using Statistics
using Tables
using ThreadsX
using MLJ
using GLM
using DecisionTree

@reexport using BiobakeryUtils
@reexport using DataFrames
@reexport using Chain
@reexport using CSV
@reexport using XLSX

const transform = DataFrames.transform

include("files.jl")
include("airtable.jl")
include("wrangle.jl")
include("metabolites.jl")
include("genefamilies.jl")
include("plotting.jl")
include("kneaddata.jl")
include("omnibus.jl")
include("prediction.jl")
include("data_loading.jl")

end
