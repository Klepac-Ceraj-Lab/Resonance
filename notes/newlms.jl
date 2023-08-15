using Resonance
using Chain
using VKCComputing


base = LocalBase()
newtaxa = CSV.read("input/ECHO_metabolomics_preselectedinputs.csv", DataFrame)
updated = CSV.read("input/fmp_alltp.csv", DataFrame)

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
