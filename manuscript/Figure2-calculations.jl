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
using MultipleTesting

const adjust = MultipleTesting.adjust
## Data Loading

mdata = Resonance.load(Metadata())

species = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)

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

specdf = comm2wide(species)
specdf.quartile = categorical(let
    l, u = quantile(specdf.cogScore, [0.25, 0.75])
    map(s-> s < l ? "lower" : s > u ? "upper" : "middle", specdf.cogScore)
end; levels=["lower", "middle", "upper"], ordered = true)

non_spec_cols = [
    "sample", "subject", "timepoint", "ageMonths", "sex", "race", "education", "date",
    "cogScore", "sample_base", "read_depth", "filter_00to120", "filter_00to06", "filter_18to120", "quartile"
    ]

stotals = map(i-> sum(specdf[i, Not(non_spec_cols)]), 1:nrow(specdf))
for spc in names(specdf, Not(non_spec_cols))
    specdf[!, spc] ./= stotals
end

ecsdf = comm2wide(ecs)
ecsdf.quartile = specdf.quartile

## GLMs

let 
    indf = specdf[specdf.filter_18to120, :]
    outfile = tablefiles("lms_species_18to120.csv")
    lmresults = DataFrame()

    for spc in names(indf, Not(non_spec_cols))
        @info spc
            
        over0 = indf[!, spc] .> 0
        sum(over0) / size(indf, 1) > 0.20 || continue
        # ab = collect(indf[!, spc] .+ (minimum(indf[over0, spc])) / 2) # add half-minimum non-zerovalue

        df = indf[:, ["ageMonths", "cogScore", "quartile", "read_depth", "education"]]
        df.bug = asin.(sqrt.(indf[!, spc]))

        mod = lm(@formula(bug ~ cogScore + ageMonths + read_depth + education), df; dropcollinear=false)
        ct = DataFrame(coeftable(mod))
        ct.species .= spc
        ct.kind .= "cogScore"
        append!(lmresults, ct)

        mod = lm(@formula(bug ~ quartile + ageMonths + read_depth + education), df; dropcollinear=false)
        ct = DataFrame(coeftable(mod))
        subset!(ct, :Name=> ByRow(q-> !contains(q, "middle")))
        droplevels!(df.quartile)
        ct.species .= spc
        ct.kind .= "quartile"
        append!(lmresults, ct)
    end

    select!(lmresults, Cols(:species, :Name, :))
    rename!(lmresults, "Pr(>|t|)"=>"pvalue");

    @chain lmresults begin
        subset!(:Name => ByRow(x->
            !any(y-> contains(x, y), 
                ("(Intercept)", "ageMonths", "read_depth", "education")
                )
            )
        )

        groupby(:kind)
        transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
        sort!(:qvalue)
    end

    CSV.write(outfile, lmresults)
    lmresults[:, Not(["Lower 95%", "Upper 95%"])]
end

#- presence/absence

let
    indf = specdf[specdf.filter_18to120, :]
    outfile = tablefiles("lms_species_18to120_pa.csv")
    lmresults = DataFrame()

    for spc in names(indf, Not(non_spec_cols))
        @info spc
            
        0.15 < prevalence(indf[!, spc]) < 0.9 || continue
        df = indf[:, ["ageMonths", "cogScore", "quartile", "read_depth", "education"]]
        df.bug = Int.(indf[!, spc] .> 0)

        
        mod = glm(@formula(bug ~ cogScore + ageMonths + read_depth + education), df, Binomial(), ProbitLink(); dropcollinear=false)
        ct = DataFrame(coeftable(mod))
        ct.species .= spc
        ct.kind .= "cogScore"
        append!(lmresults, ct)

    end

    select!(lmresults, Cols(:species, :Name, :))
    rename!(lmresults, "Pr(>|z|)"=>"pvalue");

    @chain lmresults begin
        subset!(:Name => ByRow(x->
            !any(y-> contains(x, y), 
                ("(Intercept)", "ageMonths", "read_depth", "education")
                )
            )
        )

        groupby(:kind)
        transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
        sort!(:qvalue)
    end

    CSV.write(outfile, lmresults)
    lmresults[:, Not(["Lower 95%", "Upper 95%", "kind"])]
end


#-

# ## FSEA

## Calculate correlations

unimdata = DataFrame(metadata(unirefs))

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

## Run LMs



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
    CSV.write(scratchfiles("figure2", "fsea_18to120.csv"), tmp)
end
