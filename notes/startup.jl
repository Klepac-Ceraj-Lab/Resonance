using Resonance
using AlgebraOfGraphics

omni = CSV.read("data/wrangled/omnisamples.csv", DataFrame)
etoh = CSV.read("data/wrangled/etohsamples.csv", DataFrame)

##

tps  = CSV.read("data/wrangled/timepoints.csv", DataFrame)
tps."Left-Thalamus" = map(eachrow(tps)) do row
    (t, p) = (row."Left-Thalamus", row."Left-Thalamus-Proper")
    all(ismissing, (t,p)) && return missing
    return max(coalesce(t, 0), coalesce(p, 0))
end

tps."Right-Thalamus" = map(eachrow(tps)) do row
    (t, p) = (row."Right-Thalamus", row."Right-Thalamus-Proper")
    all(ismissing, (t,p)) && return missing
    return max(coalesce(t, 0), coalesce(p, 0))
end

mainmeta = [
    "ageMonths",
    "age0to3mo",
    "age3to6mo",
    "age6to12mo",
    "age12moplus",
    "mother_HHS_Education",
    "simple_race",
    "cogScore",
    "has_segmentation",
    "ECHOTPCoded",
]

brainmeta = [
            #  "CortexVol",
            #  "CorticalWhiteMatterVol",
            #  "SubCortGrayVol",
            #  "TotalGrayVol",
            #  "BrainSegVol-to-eTIV",
            #  "CerebralWhiteMatterVol",
            #  "lhCorticalWhiteMatterVol",
            #  "lhCerebralWhiteMatterVol",
            #  "lhCortexVol",
             "Left-Thalamus",
             "Left-Lateral-Ventricle",
             "Left-Cerebellum-White-Matter",
             "Left-Cerebellum-Cortex",
             "Left-Caudate",
             "Left-Putamen",
             "Left-Pallidum",
             "Left-Hippocampus",
             "Left-Amygdala",
             "Left-Accumbens-area",
             "Left-VentralDC",
             "Left-choroid-plexus",
            #  "rhCorticalWhiteMatterVol",
            #  "rhCerebralWhiteMatterVol",
            #  "rhCortexVol",
             "Right-Thalamus",
             "Right-Lateral-Ventricle",
             "Right-Cerebellum-White-Matter",
             "Right-Cerebellum-Cortex",
             "Right-Caudate",
             "Right-Putamen",
             "Right-Pallidum",
             "Right-Hippocampus",
             "Right-Amygdala",
             "Right-Accumbens-area",
             "Right-VentralDC",
             "Right-choroid-plexus",
             "Brain-Stem",
             "CSF"
]

for m in brainmeta
    tps[!, m] .= tps[!, m] ./ map(x-> ismissing(x) || x == 0 ? missing : x, tps."EstimatedTotalIntraCranialVol")
end

select!(tps, ["subject", "timepoint", mainmeta..., brainmeta...])
rename!(tps, Dict(k=> replace(k, "-"=>"_") for k in brainmeta))
foreach(i-> (brainmeta[i] = replace(brainmeta[i], "-"=>"_")), eachindex(brainmeta))


complete_brain = completecases(tps[:, brainmeta])

##

metabolites = CSV.read("data/wrangled/metabolites.csv", DataFrame)
ms = [Resonance.Metabolite(row[:uid], row[:Metabolite], row[:MZ], row[:RT]) for row in eachrow(metabolites)]
metabolites = CommunityProfile(Matrix(metabolites[!, 9:end]), ms, MicrobiomeSample.(names(metabolites)[9:end]))
set!(metabolites, leftjoin(etoh, tps[!, ["subject",  "timepoint", mainmeta...]], on=[:subject, :timepoint], makeunique=true))

species = CSV.read("data/wrangled/species.csv", DataFrame)
species = CommunityProfile(Matrix(species[!, 2:end]), taxon.(species[!, 1]), MicrobiomeSample.(names(species)[2:end]))
set!(species, leftjoin(omni, tps[!, ["subject", "timepoint", mainmeta...]], on=[:subject, :timepoint], makeunique=true))
species = species[:, map(!ismissing, get(species, :subject)) .& map(!ismissing, get(species, :timepoint))]