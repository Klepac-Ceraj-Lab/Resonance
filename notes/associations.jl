using Resonance
using DataFramesMeta
using CairoMakie

samplemeta = airtable_metadata()
allmeta = CSV.read("data/wrangled.csv", DataFrame)
metabolites = CSV.read("data/metabolites.csv", DataFrame)
species = CSV.read("data/species", DataFrame)
unirefs = load_genefamilies()


