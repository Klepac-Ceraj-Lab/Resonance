module Resonance

export airtable_metadata,
       brain_ingest,
       findprevstool,
       count_set,
       upset_dots!,
       load_metabolites,
       pull_row,
       loadings,
       varexplained,
       commonname,
       countmap,
       codebreastfeeding!,
       stp_overlap

using Reexport
using Airtable
using CairoMakie
using AlgebraOfGraphics
using MultivariateStats
using Dictionaries
using CodecZlib
using Statistics

@reexport using BiobakeryUtils
@reexport using DataFrames
@reexport using CSV
@reexport using XLSX
@reexport using DataFramesMeta

include("airtable.jl")
include("wrangle.jl")
include("metabolites.jl")
include("genefamilies.jl")
include("plotting.jl")
include("brain.jl")

end
