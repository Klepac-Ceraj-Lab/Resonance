module Resonance

export airtable_metadata,
       brain_ingest

using Reexport
using Airtable

@reexport using Microbiome
@reexport using DataFrames
@reexport using CSV
@reexport using XLSX

include("airtable.jl")
include("wrangle.jl")

end
