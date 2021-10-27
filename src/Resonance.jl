module Resonance

export airtable_metadata,
       brain_ingest,
       count_set,
       upset_dots!

using Reexport
using Airtable
using CairoMakie

@reexport using Microbiome
@reexport using DataFrames
@reexport using CSV
@reexport using XLSX

include("airtable.jl")
include("wrangle.jl")
include("plotting.jl")

end
