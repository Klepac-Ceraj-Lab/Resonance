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
using JLD2

## Data Loading

mdata = Resonance.load(Metadata())

unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = mdata) # this can take a bit
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
allages = unique(subset(unimdata, :cogScore => ByRow(!ismissing)), :subject)
u6 = unique(subset(unimdata, :ageMonths => ByRow(<(6)), :cogScore => ByRow(!ismissing)), :subject)
o18 = unique(subset(unimdata, :ageMonths => ByRow(>(18)), :cogScore => ByRow(!ismissing)), :subject)

unirefs_00to120 = let filt = get(unirefs, :filter_00to120)
    keepuni = vec(prevalence(unirefs[:, filt]) .> 0)
    unirefs[keepuni, filt]
end

unirefs_00to06 = let filt = get(unirefs, :filter_00to06)
    keepuni = vec(prevalence(unirefs[:, filt]) .> 0)
    unirefs[keepuni, filt]
end

unirefs_18to120 = let filt = get(unirefs, :filter_18to120)
    keepuni = vec(prevalence(unirefs[:, filt]) .> 0)
    unirefs[keepuni, filt]
end


cors_00to120 = vec(cor(get(unirefs_00to120, :cogScore), abundances(unirefs_00to120), dims=2))

cors_00to06 = vec(cor(get(unirefs_00to06, :cogScore), abundances(unirefs_00to06), dims=2))
cors_18to120 = vec(cor(get(unirefs_18to120, :cogScore), abundances(unirefs_18to120), dims=2))

cors_00to120_age = vec(cor(get(unirefs_00to120, :ageMonths), abundances(unirefs_00to120), dims=2))
cors_00to06_age = vec(cor(get(unirefs_00to06, :ageMonths), abundances(unirefs_00to06), dims=2))
cors_18to120_age = vec(cor(get(unirefs_18to120, :ageMonths), abundances(unirefs_18to120), dims=2))

neuroactive_00to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs_00to120)))
neuroactive_full_00to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs_00to120)); consolidate=false)
neuroactive_00to06 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs_00to06)))
neuroactive_full_00to06 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs_00to06)); consolidate=false)
neuroactive_18to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs_18to120)))
neuroactive_full_18to120 = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs_18to120)); consolidate=false)

isdir(scratchfiles("figure2")) || mkpath(scratchfiles("figure2"))

jldsave(scratchfiles("figure2", "figure2_data.jld2");
    cors_00to120, cors_00to06, cors_18to120, 
    cors_00to120_age, cors_00to06_age, cors_18to120_age, 
    neuroactive_00to120, neuroactive_00to06, neuroactive_18to120, 
    neuroactive_full_00to120, neuroactive_full_00to06, neuroactive_full_18to120,
    unimdata
)

## Run FSEA

### All ages

let 
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive_00to120))) do gs
        ixs = neuroactive_00to120[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to120[ixs])
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to120[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(neuroactive_00to120))) do gs
        ixs = neuroactive_00to120[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to120_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to120_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(scratchfiles("figure2", "fsea_consolidated_00to120.csv"), tmp)
end


let 
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive_full_00to120))) do gs
        ixs = neuroactive_full_00to120[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to120[ixs])
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to120[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(neuroactive_full_00to120))) do gs
        ixs = neuroactive_full_00to120[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to120_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to120_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(scratchfiles("figure2", "fsea_all.csv"), tmp)

end

### Kids Under 6 months

let 
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive_00to06))) do gs
        ixs = neuroactive_00to06[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to06[ixs])
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to06[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(neuroactive_00to06))) do gs
        ixs = neuroactive_00to06[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to06_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to06_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(scratchfiles("figure2", "fsea_consolidated_00to06.csv"), tmp)
end



let 
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive_full_00to06))) do gs
        ixs = neuroactive_full_00to06[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to06[ixs])
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to06[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(neuroactive_full_00to06))) do gs
        ixs = neuroactive_full_00to06[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_00to06_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_00to06_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(scratchfiles("figure2", "fsea_u6.csv"), tmp)
end

### Kids over 18 months

let 
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive_18to120))) do gs
        ixs = neuroactive_18to120[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_18to120[ixs])
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_18to120[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(neuroactive_18to120))) do gs
        ixs = neuroactive_18to120[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_18to120_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_18to120_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(scratchfiles("figure2", "fsea_consolidated_18to120.csv"), tmp)
end


let 
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive_full_18to120))) do gs
        ixs = neuroactive_full_18to120[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_18to120[ixs])
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_18to120[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    tmp2 = DataFrame(ThreadsX.map(collect(keys(neuroactive_full_18to120))) do gs
        ixs = neuroactive_full_18to120[gs]
        isempty(ixs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cors_18to120_age[ixs])
        isempty(cs) && return (; cortest = "age",  geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cors_18to120_age[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "age", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    append!(tmp, tmp2)
    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = MultipleTesting.adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(scratchfiles("figure2", "fsea_o18.csv"), tmp)
end
