using Resonance
using CairoMakie
using Statistics
using HypothesisTests
using MultipleTesting
# using KernelDensity
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

isdir(tablefiles("figure2")) || mkpath(tablefiles("figure2"))

## Data Loading

mdata = Resonance.load(Metadata())
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 

seqs.sample = [s for s in seqs.sample]
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
species = filter(f-> taxrank(f) == :species, taxa)
relativeabundance!(species)
unirefs = Resonance.load(UnirefProfiles(); timepoint_metadata = seqs) # this can take a bit
ecs = Resonance.load(ECProfiles(); timepoint_metadata = seqs)
kos = Resonance.load(KOProfiles(); timepoint_metadata = seqs)

# metabolites = Resonance.load(MetabolicProfiles(); timepoint_metadata = mdata)
brain = Resonance.load(Neuroimaging(), timepoint_metadata = seqs, samplefield="sample")


# @assert all(samplenames(species) .== samplenames(unirefs))
@assert all(samplenames(species) .== samplenames(ecs) .== samplenames(kos) .== samplenames(unirefs))

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
    "right-insula",
    "left-insula",
    "Gray-matter",
    "White-matter"
]
brainmeta = let
    brainsub = brain[:, startswith.(samplenames(brain), "SEQ")]
    df = DataFrame(:sample => samplenames(brainsub), (Symbol(reg) => vec(abundances(brainsub[reg, :])) for reg in brain_roi)...)
end

set!(unirefs, brainmeta)

specdf = comm2wide(species)

specdf.quartile = categorical(let
l, u = quantile(skipmissing(specdf.cogScore), [0.25, 0.75])
map(s-> ismissing(s) ? missing : s < l ? "lower" : s > u ? "upper" : "middle", specdf.cogScore)
end; levels=["lower", "middle", "upper"], ordered = true)

non_spec_cols = Cols(
    "sample", "subject", "omni", "timepoint", "ageMonths", "sex", "race", "education", # @Hugemiler removed `"seqid", `
    "cogScore", "read_depth", "quartile",
    r"filter_", r"Mullen", "edfloat"
    )
    
braindf = leftjoin(brainmeta, specdf; on="sample")

ecsdf = comm2wide(ecs)
ecsdf.quartile = specdf.quartile

## GLMs

runlms(specdf[specdf.filter_00to120, :], tablefiles("figure2", "lms_species_00to120.csv"), names(specdf, Not(non_spec_cols)))
runlms(specdf[specdf.filter_18to120, :], tablefiles("figure2", "lms_species_18to120.csv"), names(specdf, Not(non_spec_cols)))
runlms(specdf[specdf.filter_00to06, :], tablefiles("figure2", "lms_species_00to06.csv"), names(specdf, Not(non_spec_cols)))

## Calculate correlations

## Create the unirefs objects

unimdata = DataFrame(get(unirefs))

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
relativeabundance!(unirefs_18to120)
relativeabundance!(unirefs_00to06)
relativeabundance!(unirefs_00to120)

## Run LMs

let 
    indf = DataFrame(get(unirefs_18to120))[:, ["ageMonths", "cogScore", "read_depth", "education"]]
    outfile = tablefiles("figure2", "lms_unirefs_18to120.csv")
    lmresults = DataFrame(ThreadsX.map(features(unirefs_18to120)) do feat
        ab = vec(abundances(unirefs_18to120[feat, :]))
        # if (prevalence(ab) > 0.2)
        #     return( NamedTuple{(Symbol("feature"),Symbol("Name"),Symbol("Coef."),Symbol("Std. Error"),Symbol("t"),Symbol("Pr(>|t|)"),Symbol("Lower 95%"),Symbol("Upper 95%"),Symbol("Cor"),Symbol("kind")),Tuple{String, String, Float64, Float64, Float64, Float64, Float64, Float64, Float64, String}}([string(feat), "cogScore", NaN, NaN, NaN, NaN, NaN, NaN, NaN, "unirefs"]) )
        # end
        # ab = collect(indf[!, feat] .+ (minimum(indf[over0, feat])) / 2) # add half-minimum non-zerovalue
        df = indf[!, :]
        df.bug = asin.(sqrt.(ab))

        mod = lm(@formula(bug ~ cogScore + ageMonths + read_depth + education), df; dropcollinear=false)
        ct = DataFrames.Tables.rowtable(coeftable(mod))[2]
        @assert ct.Name == "cogScore"
        return (; ct..., feature = string(feat), Cor = cor(df.bug, df.cogScore), kind="unirefs") # TODO: ADD THE CORRELATION (df.big and df.cogScore) to this Tuple
    end)

    select!(lmresults, Cols(:feature, :Name, :))
    rename!(lmresults, "Pr(>|t|)"=>"pvalue");

    CSV.write(outfile, lmresults)
    lmresults[:, Not(["Lower 95%", "Upper 95%"])]
end

let 
    indf = DataFrame(get(unirefs_00to120))[:, ["ageMonths", "cogScore", "read_depth", "education"]]
    outfile = tablefiles("figure2", "lms_unirefs_00to120.csv")
    lmresults = DataFrame(ThreadsX.map(features(unirefs_00to120)) do feat
        ab = vec(abundances(unirefs_00to120[feat, :]))
        if (prevalence(ab) > 0.2)
            return( NamedTuple{(Symbol("feature"),Symbol("Name"),Symbol("Coef."),Symbol("Std. Error"),Symbol("t"),Symbol("Pr(>|t|)"),Symbol("Lower 95%"),Symbol("Upper 95%"),Symbol("kind")),Tuple{String, String, Float64, Float64, Float64, Float64, Float64, Float64, String}}([string(feat), "cogScore", NaN, NaN, NaN, NaN, NaN, NaN, "unirefs"]) )
        end
        # ab = collect(indf[!, feat] .+ (minimum(indf[over0, feat])) / 2) # add half-minimum non-zerovalue
        df = indf[!, :]
        df.bug = asin.(sqrt.(ab))

        mod = lm(@formula(bug ~ cogScore + ageMonths + read_depth + education), df; dropcollinear=false)
        ct = DataFrames.Tables.rowtable(coeftable(mod))[2]
        @assert ct.Name == "cogScore"
        return (; ct..., feature = string(feat), Cor = cor(df.bug, df.cogScore), kind="unirefs") # TODO: ADD THE CORRELATION (df.big and df.cogScore) to this Tuple
    end)

    select!(lmresults, Cols(:feature, :Name, :))
    rename!(lmresults, "Pr(>|t|)"=>"pvalue");

    CSV.write(outfile, lmresults)
    lmresults[:, Not(["Lower 95%", "Upper 95%"])]
end

let 
    indf = DataFrame(get(unirefs_00to06))[:, ["ageMonths", "cogScore", "read_depth", "education"]]
    outfile = tablefiles("figure2", "lms_unirefs_00to06.csv")
    lmresults = DataFrame(ThreadsX.map(features(unirefs_00to06)) do feat
        ab = vec(abundances(unirefs_00to06[feat, :]))
        if (prevalence(ab) > 0.2)
            return( NamedTuple{(Symbol("feature"),Symbol("Name"),Symbol("Coef."),Symbol("Std. Error"),Symbol("t"),Symbol("Pr(>|t|)"),Symbol("Lower 95%"),Symbol("Upper 95%"),Symbol("kind")),Tuple{String, String, Float64, Float64, Float64, Float64, Float64, Float64, String}}([string(feat), "cogScore", NaN, NaN, NaN, NaN, NaN, NaN, "unirefs"]) )
        end
        # ab = collect(indf[!, feat] .+ (minimum(indf[over0, feat])) / 2) # add half-minimum non-zerovalue
        df = indf[!, :]
        df.bug = asin.(sqrt.(ab))

        mod = lm(@formula(bug ~ cogScore + ageMonths + read_depth + education), df; dropcollinear=false)
        # @show DataFrames.Tables.rowtable(coeftable(mod))
        ct = DataFrames.Tables.rowtable(coeftable(mod))[2]
        @assert ct.Name == "cogScore"

        return (; ct..., feature = string(feat), Cor = cor(df.bug, df.cogScore), kind="unirefs") # TODO: ADD THE CORRELATION (df.big and df.cogScore) to this Tuple
    end)

    select!(lmresults, Cols(:feature, :Name, :))
    rename!(lmresults, "Pr(>|t|)"=>"pvalue");

    CSV.write(outfile, lmresults)
    lmresults[:, Not(["Lower 95%", "Upper 95%"])]
end

#-

## Run FSEA

### All ages

let infile = tablefiles("figure2", "lms_unirefs_00to120.csv")
    outfile = tablefiles("figure2", "fsea_consolidated_00to120.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature)) # I needed ENV["DATA_FILES"] = "/brewster/kevin/scratch/raw_data" to run this
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end

let infile = tablefiles("figure2", "lms_unirefs_00to120.csv")
    outfile = tablefiles("figure2", "fsea_all_00to120.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature); consolidate=false)
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end

#-

### Kids Under 6 months

u6cors = vec(cor(get(unirefs_00to06, :cogScore), abundances(unirefs_00to06), dims=2))

let infile = tablefiles("figure2", "lms_unirefs_00to06.csv")
    outfile = tablefiles("figure2", "fsea_consolidated_00to06.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = u6cors #Ts = df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature))
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end

let infile = tablefiles("figure2", "lms_unirefs_00to06.csv")
    outfile = tablefiles("figure2", "fsea_all_00to06.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature); consolidate=false)
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end
### Kids over 18 months

o18cors = vec(cor(get(unirefs_18to120, :cogScore), abundances(unirefs_18to120), dims=2))

let infile = tablefiles("figure2", "lms_unirefs_18to120.csv")
    outfile = tablefiles("figure2", "fsea_consolidated_18to120.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = df.Cor # df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature))
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end

let infile = tablefiles("figure2", "lms_unirefs_18to120.csv")
    outfile = tablefiles("figure2", "fsea_all_18to120.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = df.Cor # df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature); consolidate=false)
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        es = Resonance.enrichment_score(cs, acs)

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end
