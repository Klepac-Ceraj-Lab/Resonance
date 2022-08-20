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
       startup

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
export build_future_df,
       check_longdata_metaduplicates!,
       microbiome_predictors,
       tryparsecol,
       univariate_tietjenmoore,
       test_tietjenmoore,
       try_outliers,
       unstack_techreplicates,
       build_metadata_prediction_df

using Reexport
using ReTest

using Arrow
using Airtable
using AlgebraOfGraphics
using CairoMakie
using CategoricalArrays
using CodecZlib
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
include("startup.jl")
include("kneaddata.jl")
include("omnibus.jl")
include("prediction.jl")
include("data_loading.jl")

end
