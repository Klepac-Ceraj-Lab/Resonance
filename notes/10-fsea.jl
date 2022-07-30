using Resonance
using ProgressLogging
using CategoricalArrays
using CairoMakie
using AlgebraOfGraphics
using Microbiome.Distances
using Microbiome.MultivariateStats
using HypothesisTests
using MultipleTesting
using AlgebraOfGraphics
using Statistics
using DataFrames.InvertedIndices
using DataFrames.PrettyTables
using ThreadsX

omni, etoh, tps, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()
unique!(tps, ["subject", "timepoint"])

set!(species, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(genes, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(metabolites, leftjoin(select(etoh, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))

neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(genes)))

#-

gsgenes = genes[:, .!ismissing.(get(genes, :cogScore))]

cscor = cor(abundances(gsgenes), get(gsgenes, :cogScore), dims=2)

fsdf = DataFrame(
    map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        fig = Resonance.plot_fsea(acs, cs; label = gs)
        save("figures/fsea_$(replace(gs, ' '=>"-")).png", fig)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write(datafiles("fsea_all.csv"), fsdf)

pretty_table(first(fsdf, 10); backend = Val(:latex))
#- 

nodupes_samples = unique(DataFrame(metadata(gsgenes)), :subject).sample

gsnodupes = gsgenes[:, map(s-> name(s) ∈ nodupes_samples, samples(gsgenes))]

cscor = cor(abundances(gsnodupes), get(gsnodupes, :cogScore), dims=2)

fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write(datafiles("fsea_nodupe.csv"))

#-

u6_samples = unique(subset(DataFrame(metadata(gsgenes)), :ageMonths => ByRow(<(6))), :subject).sample

gsu6 = gsgenes[:, map(s-> name(s) ∈ u6_samples, samples(gsgenes))]

cscor = cor(abundances(gsu6), get(gsu6, :cogScore), dims=2)

fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write(datafiles("fsea_u6.csv"))


#-

o12_samples = unique(subset(DataFrame(metadata(gsgenes)), :ageMonths => ByRow(>(12))), :subject).sample

gso12 = gsgenes[:, map(s-> name(s) ∈ o12_samples, samples(gsgenes))]

cscor = cor(abundances(gso12), get(gso12, :cogScore), dims=2)

fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write(datafiles("fsea_o12.csv"))


#-

brgenes = genes[:, .!ismissing.(get(genes, Symbol("Left-Thalamus")))]

brfsea = DataFrame()

for region in Resonance.brainmeta
    @warn region
    cscor = cor(abundances(brgenes), get(brgenes, Symbol(region)), dims=2)

    fsdf = DataFrame(
        map(collect(keys(neuroactive))) do gs
            @info gs
            ixs = neuroactive[gs]
            isempty(ixs) && return (; region, geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            cs = filter(!isnan, cscor[ixs])
            isempty(cs) && return (; region, geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

            acs = filter(!isnan, cscor[Not(ixs)])
            mwu = MannWhitneyUTest(cs, acs)

            return (; region, geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
        end
    )
    subset!(fsdf, :pvalue=> ByRow(!isnan))
    append!(brfsea, fsdf)
end


brfsea.qvalue = adjust(brfsea.pvalue, BenjaminiHochberg())
sort!(brfsea, :qvalue)
CSV.write(datafiles("fsea_brain_all.csv"), brfsea)

pretty_table(first(brfsea, 10); backend = Val(:latex))

#- 

