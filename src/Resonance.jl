module Resonance

export airtable_metadata,
       brain_ingest,
       count_set,
       upset_dots!,
       load_metabolites,
       pull_row

using Reexport
using Airtable
using CairoMakie
using AlgebraOfGraphics

@reexport using BiobakeryUtils
@reexport using DataFrames
@reexport using CSV
@reexport using XLSX
@reexport using Chain
@reexport using DataFramesMeta

include("airtable.jl")
include("wrangle.jl")
include("metabolites.jl")
include("plotting.jl")
include("brain.jl")

end
