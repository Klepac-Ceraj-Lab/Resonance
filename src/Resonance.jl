module Resonance

export airtable_metadata,
       brain_ingest,
       count_set,
       upset_dots!

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
include("plotting.jl")

end
