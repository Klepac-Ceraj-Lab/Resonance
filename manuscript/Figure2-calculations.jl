using Resonance
using CairoMakie
using Statistics
using HypothesisTests
using MultipleTesting
using KernelDensity
using MultivariateStats
using Distributions
using CategoricalArrays
using ThreadsX
using ColorSchemes
using GLM
using Dates

## Data Loading

mdata = Resonance.load(Metadata())

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
unirefs = unirefs[:, map(s-> !ismissing(s) && s < 120, get(unirefs, :ageMonths))]
unistrat =  filter(f-> hastaxon(f) || name(f) == "UNMAPPED", unirefs)
unirefs = filter(!hastaxon, unirefs) # don't use species stratification for summaries

metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata=mdata)
metabolites = metabolites[:, [!ismissing(a) && a < 14 for a in get(metabolites, :ageMonths)]]
isdefined(Main, :metdm) || (metdm = braycurtis(metabolites))
metpco = fit(MDS, metdm; distances=true)

brain = Resonance.load(Neuroimaging(); timepoint_metadata=mdata)

brain_roi = [
    "right-lateral-occipital",
    "left-lateral-occipital",
    "right-inferior-parietal",
    "left-inferior-parietal",
    "right-middle-temporal",
    "left-middle-temporal",
    "right-cerebellum-white-matter",
    "left-cerebellum-white-matter",
    "right-thalamus-proper",
    "left-thalamus-proper",
]
brainmeta = let
    brainsub = brain[:, startswith.(samplenames(brain), "FG")]
    df = DataFrame(:sample => samplenames(brainsub), (Symbol(reg) => vec(abundances(brainsub[reg, :])) for reg in brain_roi)...)
end

set!(unirefs, brainmeta)

## Calculate correlations

unimdata = DataFrame(metadata(unirefs))
allages = unique(subset(unimdata, :cogScorePercentile => ByRow(!ismissing)), :subject)
u6 = unique(subset(unimdata, :ageMonths => ByRow(<(6)), :cogScorePercentile => ByRow(!ismissing)), :subject)
o18 = unique(subset(unimdata, :ageMonths => ByRow(>(18)), :cogScorePercentile => ByRow(!ismissing)), :subject)

allcomm = let keepuni = vec(prevalence(unirefs[:, allages.sample]) .> 0)
    unirefs[keepuni, allages.sample]
end

u6comm = let keepuni = vec(prevalence(unirefs[:, u6.sample]) .> 0)
    unirefs[keepuni, u6.sample]
end

o18comm = let keepuni = vec(prevalence(unirefs[:, o18.sample]) .> 0)
    unirefs[keepuni, o18.sample]
end

allcors = vec(cor(get(allcomm, :cogScorePercentile), abundances(allcomm), dims=2))

u6cors = vec(cor(get(u6comm, :cogScorePercentile), abundances(u6comm), dims=2))
o18cors = vec(cor(get(o18comm, :cogScorePercentile), abundances(o18comm), dims=2))

allcors_age = vec(cor(get(allcomm, :ageMonths), abundances(allcomm), dims=2))
u6cors_age = vec(cor(get(u6comm, :ageMonths), abundances(u6comm), dims=2))
o18cors_age = vec(cor(get(o18comm, :ageMonths), abundances(o18comm), dims=2))

all_neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(allcomm)))
all_neuroactive_full = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(allcomm)); consolidate=false)
u6_neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(u6comm)))
u6_neuroactive_full = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(u6comm)); consolidate=false)
o18_neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(o18comm)))
o18_neuroactive_full = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(o18comm)); consolidate=false)

## Run FSEA

### All ages

let 
    tmp = DataFrame(ThreadsX.map(collect(keys(all_neuroactive))) do gs
        ixs = all_neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, allcors[ixs])
        isempty(cs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, allcors[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScorePercentile", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(all_neuroactive))) do gs
        ixs = all_neuroactive[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, allcors_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, allcors_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outputfiles("fsea_all_consolidated.csv"), tmp)
end


let 
    tmp = DataFrame(ThreadsX.map(collect(keys(all_neuroactive_full))) do gs
        ixs = all_neuroactive_full[gs]
        isempty(ixs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, allcors[ixs])
        isempty(cs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, allcors[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScorePercentile", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(all_neuroactive_full))) do gs
        ixs = all_neuroactive_full[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, allcors_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, allcors_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outputfiles("fsea_all.csv"), tmp)

end

### Kids Under 6 months

let 
    tmp = DataFrame(ThreadsX.map(collect(keys(u6_neuroactive))) do gs
        ixs = u6_neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, u6cors[ixs])
        isempty(cs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, u6cors[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScorePercentile", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(u6_neuroactive))) do gs
        ixs = u6_neuroactive[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, u6cors_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, u6cors_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outputfiles("fsea_u6_consolidated.csv"), tmp)
end



let 
    tmp = DataFrame(ThreadsX.map(collect(keys(u6_neuroactive_full))) do gs
        ixs = u6_neuroactive_full[gs]
        isempty(ixs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, u6cors[ixs])
        isempty(cs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, u6cors[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScorePercentile", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(u6_neuroactive_full))) do gs
        ixs = u6_neuroactive_full[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, u6cors_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, u6cors_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outputfiles("fsea_u6.csv"), tmp)
end

### Kids over 18 months

let 
    tmp = DataFrame(ThreadsX.map(collect(keys(o18_neuroactive))) do gs
        ixs = o18_neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, o18cors[ixs])
        isempty(cs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, o18cors[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScorePercentile", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(o18_neuroactive))) do gs
        ixs = o18_neuroactive[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, o18cors_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, o18cors_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outputfiles("fsea_o18_consolidated.csv"), tmp)
end


let 
    tmp = DataFrame(ThreadsX.map(collect(keys(o18_neuroactive_full))) do gs
        ixs = o18_neuroactive_full[gs]
        isempty(ixs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, o18cors[ixs])
        isempty(cs) && return (; cortest = "cogScorePercentile", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, o18cors[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScorePercentile", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(o18_neuroactive_full))) do gs
        ixs = o18_neuroactive_full[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, o18cors_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, o18cors_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outputfiles("fsea_o18.csv"), tmp)
end
