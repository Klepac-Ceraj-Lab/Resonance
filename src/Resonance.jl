module Resonance

export airtable_metadata,
       brain_ingest,
       count_set,
       upset_dots!,
       load_metabolites,
       pull_row,
       loadings,
       varexplained,
       commonname

using Reexport
using Airtable
using CairoMakie
using AlgebraOfGraphics
using MultivariateStats
using Dictionaries
using CodecZlib

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
