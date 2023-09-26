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
    "sample", "subject", "seqid", "omni", "timepoint", "ageMonths", "sex", "race", "education",
    "cogScore", "read_depth", "quartile",
    r"filter_", r"Mullen", "edfloat"
    )
    
braindf = leftjoin(brainmeta, specdf; on="sample")

ecsdf = comm2wide(ecs)
ecsdf.quartile = specdf.quartile

## GLMs

runlms(species[:, get(species, :filter_00to120)], tablefiles("figure2", "lms_species_00to120.csv"))
runlms(species[:, get(species, :filter_18to120)], tablefiles("figure2", "lms_species_18to120.csv"))
runlms(species[:, get(species, :filter_00to06)], tablefiles("figure2", "lms_species_00to06.csv"))

let 
    indf = unique(sort(braindf, [:subject, :timepoint]; rev=true), :subject)
    outfile = tablefiles("figure2", "lms_species_brain.csv")
    lmresults = DataFrame()

    for roi in brain_roi
        @info roi
        for spc in names(indf, names(specdf, r"^[A-Z][a-z]+_[a-z0-9_]+"))
                
            over0 = indf[!, spc] .> 0
            sum(over0) / size(indf, 1) > 0.20 || continue
            # ab = collect(indf[!, spc] .+ (minimum(indf[over0, spc])) / 2) # add half-minimum non-zerovalue

            df = indf[:, ["ageMonths", "read_depth", "education"]]
            df.bug = asin.(sqrt.(indf[!, spc]))
            df.brain = indf[!, roi]

            mod = lm(@formula(bug ~ brain + ageMonths + read_depth + education), df; dropcollinear=false)
            ct = DataFrame(coeftable(mod))
            ct.species .= spc
            ct.kind .= roi
            append!(lmresults, ct)

        end
    end

    select!(lmresults, Cols(:species, :kind, :))
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

runlms(species[:, get(species, :filter_18to120)], tablefiles("figure2", "lms_species_18to120_pa.csv"); model_kind=:logistic, prevalence_threshold=0.15)

#-

# ## FSEA

## Calculate correlations

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

relativeabundance!(unirefs_00to120)
relativeabundance!(unirefs_18to120)
relativeabundance!(unirefs_00to06)




## Run LMs

runlms(unirefs_00to120, tablefiles("figure2", "lms_unirefs_00to120.csv"); prevalence_threshold = 0.1, model_kind=:logistic)
runlms(unirefs_18to120, tablefiles("figure2", "lms_unirefs_00to120.csv"); prevalence_threshold = 0.1, model_kind=:logistic)
runlms(unirefs_00to06, tablefiles("figure2", "lms_unirefs_00to120.csv"); prevalence_threshold = 0.1, model_kind=:logistic)

#-

## Run FSEA

### All ages

let infile = tablefiles("figure2", "lms_unirefs_00to120.csv")
    outfile = tablefiles("figure2", "fsea_consolidated_00to120.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    stats = df.stat
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature))
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, enrichment = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, stats[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, enrichment = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, stats[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        # mwu = MannWhitneyUTest(cs, acs)
        
        (pseudop, es) = Resonance.fsea_permute([cs; acs], 1:length(cs))

        return (; cortest = "cogScore", geneset = gs, enrichment = es,  pvalue=pseudop)
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end

let infile = tablefiles("figure2", "lms_unirefs_00to120.csv")
    outfile = tablefiles("figure2", "fsea_all_00to120.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    stats = df.stat
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature); consolidate=false)
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, enrichment = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, stats[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, enrichment = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, stats[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        # mwu = MannWhitneyUTest(cs, acs)
        
        (pseudop, es) = Resonance.fsea_permute([cs; acs], 1:length(cs))

        return (; cortest = "cogScore", geneset = gs, enrichment = es,  pvalue=pseudop)
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end

#-

### Kids Under 6 months

let infile = tablefiles("figure2", "lms_unirefs_00to06.csv")
    outfile = tablefiles("figure2", "fsea_consolidated_00to06.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature))
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        # es = Resonance.enrichment_score(cs, acs) # needs updating - will break

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
        # es = Resonance.enrichment_score(cs, acs) # needs updating - will break

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end
### Kids over 18 months

let infile = tablefiles("figure2", "lms_unirefs_18to120.csv")
    outfile = tablefiles("figure2", "fsea_consolidated_18to120.csv")
    df = subset(CSV.read(infile, DataFrame), "pvalue"=> ByRow(!isnan))
    Ts = df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature))
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        # es = Resonance.enrichment_score(cs, acs) # needs updating - will break

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
    Ts = df.t
    neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), df.feature); consolidate=false)
    tmp = DataFrame(ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(x-> !isnan(x) && x < 7, Ts[ixs]) # Some *very* small coeficients have spuriously large T-stats
        isempty(cs) && return (; cortest = "cogScore", geneset = gs, U = NaN, median = NaN, enrichment = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(x-> !isnan(x) && x < 7, Ts[Not(ixs)]) # Some *very* small coeficients have spuriously large T-stats
        mwu = MannWhitneyUTest(cs, acs)
        # es = Resonance.enrichment_score(cs, acs) # needs updating - will break

        return (; cortest = "cogScore", geneset = gs, U = mwu.U, median = mwu.median, enrichment = es, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end)

    subset!(tmp, :pvalue=> ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, :qvalue)
    CSV.write(outfile, tmp)
end
