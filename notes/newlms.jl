using Resonance
using Chain
using VKCComputing


base = LocalBase()
newtaxa = CSV.read("input/ECHO_metabolomics_preselectedinputs.csv", DataFrame)
updated = CSV.read("input/fmp_alltp.csv", DataFrame)
transform!(updated, AsTable(["ageMonths", "assessmentAgeDays", "assessmentAgeMonths", "scanAgeDays", "scanAgeMonths"]) => ByRow(row-> begin
        am, aad, aam, sad, sam = row
        calcage = coalesce(aad / 365 * 12 + aam, sad / 365 * 12 + sam)
        if !ismissing(am)
            !ismissing(calcage) && am != calcage && @info am
        end

        return coalesce(am, calcage)
        end) => "ageMonths"
)

function get_batches(base, bsp_id)
    bsp = get(base["Biospecimens"], bsp_id, nothing)
    isnothing(bsp) && throw(KeyError("Biospecimen id $bsp_id doesn't exist in database"))

    batches = get(bsp, Symbol("sequencing_batch (from seqprep)"), nothing)
    isnothing(batches) && return String[]
    
    return unique([b[:uid] for b in base[batches]])
end

@assert all(bs-> length(bs) == 1, get_batches.(Ref(base), String.(newtaxa.omni)))

newtaxa.batch = first.(get_batches.(Ref(base), String.(newtaxa.omni)))

b19 = subset(newtaxa, "batch"=> ByRow(==("mgx025")))
b19 = leftjoin(select(b19, Not(["ageMonths", "cogScore"])), select(updated, "subject", "timepoint", "ageMonths", "cogScore"); on=["subject", "timepoint"])

subset!(b19, AsTable(["ageMonths", "cogScore"])=> ByRow(row-> !any(ismissing, row)))

#-

using GLM

lms = DataFrame()

for sp in [
    "Erysipelatoclostridium_ramosum",
    "Eggerthella_lenta",
    "Escherichia_coli",
    "Bifidobacterium_longum",
    "Veillonella_parvula",
    "Bifidobacterium_breve",
    "Klebsiella_variicola",
    "Enterococcus_faecalis",
    "Streptococcus_salivarius",
    "Klebsiella_pneumoniae",
    "Bifidobacterium_bifidum",
    "Flavonifractor_plautii",
    "Ruminococcus_gnavus",
    "Veillonella_atypica",
    "Klebsiella_quasipneumoniae",
    "Gordonibacter_pamelaeae",
    "Clostridioides_difficile",
    "Clostridium_innocuum",
    "Streptococcus_parasanguinis",
]
    indf = select(b19, "ageMonths", "cogScore", "education", sp)

    
    over0 = indf[!, sp] .> 0
    prev = sum(over0) / size(indf, 1)
    ab = collect(indf[!, sp] .+ (minimum(indf[over0, sp])) / 2) # add half-minimum non-zerovalue

    df = select(indf, "ageMonths", "cogScore", "education"; copycols=false)
    df.bug = log.(ab)

    mod = lm(@formula(bug ~ cogScore + ageMonths + education), df; dropcollinear=false)
    ct = DataFrame(coeftable(mod))
    ct.species .= sp
    ct.kind .= "cogScore"
    ct.prevalence .= prev
    append!(lms, subset!(ct, "Name"=> ByRow(==("cogScore"))))
end

#-

using MixedModels
using CategoricalArrays
using MultipleTesting

medf = subset(newtaxa, "ageMonths"=> ByRow(a-> !ismissing(a) && a > 18), "cogScore"=> ByRow(!ismissing))
medf.subject = categorical(medf.subject)
length(unique(medf.subject))

nsamples = @chain medf begin
    groupby("subject")
    DataFrames.combine("timepoint"=> length => "nsamples")
end

count(>(1), nsamples.nsamples)

memods = DataFrame()

# [
#     "Erysipelatoclostridium_ramosum",
#     "Eggerthella_lenta",
#     "Escherichia_coli",
#     "Bifidobacterium_longum",
#     "Veillonella_parvula",
#     "Bifidobacterium_breve",
#     "Klebsiella_variicola",
#     "Enterococcus_faecalis",
#     "Streptococcus_salivarius",
#     "Klebsiella_pneumoniae",
#     "Bifidobacterium_bifidum",
#     "Flavonifractor_plautii",
#     "Ruminococcus_gnavus",
#     "Veillonella_atypica",
#     "Klebsiella_quasipneumoniae",
#     "Gordonibacter_pamelaeae",
#     "Clostridioides_difficile",
#     "Clostridium_innocuum",
#     "Streptococcus_parasanguinis",
# ]

for sp in names(medf, r"[A-Z][a-z]+_[a-z]+")
    indf = select(medf, "ageMonths", "cogScore", "education", "subject", sp)

    
    over0 = indf[!, sp] .> 0
    prev = sum(over0) / size(indf, 1)
    prev < 0.10 && continue
    ab = asin.(sqrt.(indf[!, sp] ./ sum(indf[!, sp])))

    df = select(indf, "ageMonths", "cogScore", "education", "subject"; copycols=false)
    df.bug = ab

    mod = fit(MixedModel, @formula(bug ~ cogScore + ageMonths + education + (1|subject)), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= sp
    ct.kind .= "cogScore"
    ct.prevalence .= prev
    append!(memods, subset!(ct, "Name"=> ByRow(==("cogScore"))))
end

memods."qvalue" = adjust(memods."Pr(>|z|)", BenjaminiHochberg())
sort!(memods.qvalue)

#-

using MixedModels
using CategoricalArrays
using MultipleTesting
using CairoMakie


msel_subs = names(updated, r"^Mullen::.+T$")
wppsi_subs = names(updated, r"^Wppsi_IV::.+Composite$")
wisc_subs = names(updated, r"^Wisc_V::.+Composite$")

mseldf = @chain updated begin
    subset("ageMonths"=> ByRow(!ismissing), AsTable(msel_subs)=> ByRow(row-> !any(ismissing, row)))
    select("ageMonths", "education", "subject", "timepoint", msel_subs...)
    leftjoin!(unique(select(newtaxa, Not("ageMonths", "education")), ["subject", "timepoint"]); on = ["subject", "timepoint"])
    subset!("Bifidobacterium_longum"=> ByRow(!ismissing))
    sort("timepoint")
    unique(["subject"])
    transform("subject"=> categorical=> "subject")
end

hist(mseldf.ageMonths)

#-

mulms = DataFrame()

# for sp in [
#         "Erysipelatoclostridium_ramosum",
#         "Eggerthella_lenta",
#         "Escherichia_coli",
#         "Bifidobacterium_longum",
#         "Veillonella_parvula",
#         "Bifidobacterium_breve",
#         "Klebsiella_variicola",
#         "Enterococcus_faecalis",
#         "Streptococcus_salivarius",
#         "Klebsiella_pneumoniae",
#         "Bifidobacterium_bifidum",
#         "Flavonifractor_plautii",
#         "Ruminococcus_gnavus",
#         "Veillonella_atypica",
#         "Klebsiella_quasipneumoniae",
#         "Gordonibacter_pamelaeae",
#         "Clostridioides_difficile",
#         "Clostridium_innocuum",
#         "Streptococcus_parasanguinis",
#     ]
for sp in names(mseldf, r"^[A-Z][a-z]+_\w+")
    indf = select(mseldf, "ageMonths", "education", "subject", sp)

    
    over0 = indf[!, sp] .> 0
    prev = sum(over0) / size(indf, 1)
    prev < 0.10 && continue
    ab = asin.(sqrt.(indf[!, sp] ./ sum(indf[!, sp])))
    
    for scale in msel_subs

        df = select(indf, "ageMonths", "education", "subject"; copycols=false)
        df.scale = mseldf[!, scale]
        count(!ismissing, df.scale) < 20 && continue
        df.bug = ab

        mod = fit(MixedModel, @formula(bug ~ scale + ageMonths + education + (1|subject)), df)
        ct = DataFrame(coeftable(mod))
        ct.species .= sp
        ct.kind .= scale
        ct.prevalence .= prev
        append!(mulms, subset!(ct, "Name"=> ByRow(==("scale"))))
    end
end

@chain mulms begin
    rename!("Pr(>|z|)"=>"pvalue")
    groupby("kind")
    transform!("pvalue"=> (pval -> adjust(collect(pval), BenjaminiHochberg()))=> "qvalue")
    sort!("qvalue")
end


#-
wiscdf = @chain updated begin
    subset("ageMonths"=> ByRow(!ismissing), AsTable(wisc_subs)=> ByRow(row-> any(!ismissing, row)))
    select("ageMonths", "education", "subject", "timepoint", wisc_subs...)
    leftjoin!(unique(select(newtaxa, Not("ageMonths", "education")), ["subject", "timepoint"]); on = ["subject", "timepoint"])
    subset!("Bifidobacterium_longum"=> ByRow(!ismissing))
    sort("timepoint")
    unique(["subject"])
    transform("subject"=> categorical=> "subject")
end

wslms = DataFrame()

for sp in [
        "Erysipelatoclostridium_ramosum",
        "Eggerthella_lenta",
        "Escherichia_coli",
        "Bifidobacterium_longum",
        "Veillonella_parvula",
        "Bifidobacterium_breve",
        "Klebsiella_variicola",
        "Enterococcus_faecalis",
        "Streptococcus_salivarius",
        "Klebsiella_pneumoniae",
        "Bifidobacterium_bifidum",
        "Flavonifractor_plautii",
        "Ruminococcus_gnavus",
        "Veillonella_atypica",
        "Klebsiella_quasipneumoniae",
        "Gordonibacter_pamelaeae",
        "Clostridioides_difficile",
        "Clostridium_innocuum",
        "Streptococcus_parasanguinis",
    ]
# for sp in names(wiscdf, r"^[A-Z][a-z]+_[a-z0-9]+(_[a-z0-9]+)?$")
    @info sp
    indf = select(wiscdf, "ageMonths", "education", "subject", sp)
    
    over0 = indf[!, sp] .> 0
    prev = sum(over0) / size(indf, 1)
    prev < 0.10 && continue
    ab = asin.(sqrt.(indf[!, sp] ./ sum(indf[!, sp])))
    
    for scale in wisc_subs

        df = select(indf, "ageMonths", "education", "subject"; copycols=false)
        df.scale = wiscdf[!, scale]
        count(!ismissing, df.scale) < 20 && continue
        df.bug = ab

        mod = fit(MixedModel, @formula(bug ~ scale + ageMonths + education + (1|subject)), df)
        ct = DataFrame(coeftable(mod))
        ct.species .= sp
        ct.kind .= scale
        ct.prevalence .= prev
        append!(wslms, subset!(ct, "Name"=> ByRow(==("scale"))))
    end
end

@chain wslms begin
    rename!("Pr(>|z|)"=>"pvalue")
    groupby("kind")
    transform!("pvalue"=> (pval -> adjust(collect(pval), BenjaminiHochberg()))=> "qvalue")
    sort!("qvalue")
end


