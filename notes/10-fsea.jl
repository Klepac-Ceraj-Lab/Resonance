using Resonance
using ProgressLogging
using CategoricalArrays
using CairoMakie
using AlgebraOfGraphics
using Microbiome.Distances
using Microbiome.MultivariateStats
using MultipleTesting
using AlgebraOfGraphics

omni, etoh, tps, complete_brain, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()
unique!(tps, ["subject", "timepoint"])

set!(species, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(genes, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(metabolites, leftjoin(select(etoh, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
