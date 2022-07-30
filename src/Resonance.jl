module Resonance

export datafiles,
       scratchfiles,
       analysisfiles,
       airtable_metadata,
       brain_ingest,
       findprevstool,
       count_set,
       upset_dots!,
       load_metabolites,
       pull_row,
       loadings,
       varexplained,
       mdsaxis,
       permanovas,
       plot_permanovas,
       mantel,
       commonname,
       countmap,
       codebreastfeeding!,
       stp_overlap,
       comm_overlap,
       startup,
       build_future_df,
       check_longdata_metaduplicates!,
       microbiome_predictors,
       tryparsecol,
       univariate_tietjenmoore,
       test_tietjenmoore,
       try_outliers,
       unstack_techreplicates

using Reexport
using Airtable
using CairoMakie
using AlgebraOfGraphics
using MultivariateStats
using Dictionaries
using CodecZlib
using Statistics
using Arrow
using Tables
using SparseArrays
using CodecZlib
using FilePaths
using CategoricalArrays
using PERMANOVA
using LinearAlgebra
using Random
using ThreadsX

@reexport using BiobakeryUtils
@reexport using DataFrames
@reexport using CSV
@reexport using XLSX

const transform = DataFrames.transform

include("files.jl")
include("airtable.jl")
include("wrangle.jl")
include("metabolites.jl")
include("genefamilies.jl")
include("plotting.jl")
include("brain.jl")
include("startup.jl")
include("kneaddata.jl")
include("mantel.jl")
include("prediction.jl")

end
