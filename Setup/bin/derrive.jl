using Resonance
using Chain
using CairoMakie
using AlgebraOfGraphics

loaddata() = quote
mdata = Resonance.load_raw_metadata()
taxa = Resonance.load_raw_metaphlan()
ecs = Resonance.load_raw_humann(; kind = "ecs", names=true, stratified=false)
kos = Resonance.load_raw_humann(; kind = "kos", names=true, stratified=false)
# unirefs = Resonance.load_raw_humann(; kind = "genefamilies", names=false, stratified=false)
@assert samplenames(kos) == samplenames(ecs) == samplenames(taxa)
mdata.isMom = map(tp-> !ismissing(tp) && tp != "-8" && contains(tp, r"^pre"i), mdata.ECHOTPCoded)
mdata.isKid = map(tp-> !ismissing(tp) && tp != "-8" && !contains(tp, r"^pre"i), mdata.ECHOTPCoded)
end
#- Samples




# #- mgx / mbx


# @info "samples with techreps:"
# @info @chain samplenames(taxa) begin
#     DataFrame(; sample=_)
#     transform("sample"=> ByRow(s-> replace(s, r"_S\d+"=>"")) => "sample_base")
#     groupby("sample_base")
#     combine("sample"=> length=> "n_reps")
#     subset("n_reps"=> ByRow(>(1)))
# end

# @info """
# ## Raw Data

# - N rows in metadata (subject=> timepoint pairs): $(nrow(mdata))
# - N sequenced samples with taxaonomic profiles: $(nsamples(taxa))
# - N sequenced samples with functional profiles: $(nsamples(ecs))
# - N unique child stool samples: $(
#     nrow(subset(mdata, "omni"=> ByRow(!ismissing), "isKid" => identity)))
# - N unique maternal stool samples: $(
#     nrow(subset(mdata, "omni"=> ByRow(!ismissing), "isMom" => identity)))
# - N child subjects with at least 1 stool sample: $(
#     nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity), "subject"))
# )
#   - N child subjects with stool sample < 6 months: $(
#     nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
#                           "ageMonths" => ByRow(a-> !ismissing(a) && a < 6)), "subject"))

#   )
#   - N child subjects with stool sample < 8 months: $(
#     nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
#                           "ageMonths" => ByRow(a-> !ismissing(a) && a < 8)), "subject"))

#   )
#   - N child subjects with stool sample < 12 months: $(
#     nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
#                           "ageMonths" => ByRow(a-> !ismissing(a) && a < 12)), "subject"))

#   )
#   - N child subjects with stool sample > 18 months: $(
#     nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
#                           "ageMonths" => ByRow(a-> !ismissing(a) && a > 18)), "subject"))

#   )
# - N maternal subjects with at least 1 stool sample: $(
#     nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isMom"=> identity), "subject"))
# )

# """

# #- Important Metadata (kids)

# keep_fields = ["subject", "timepoint", "ageMonths", "education", "race", "cogScore", "sex"
#                "isKid", "isMom", "omni", "omni_collectionDate", "etoh"]
# reduced = @chain mdata begin
#     select(keep_fields)
#     rename("omni"=> "sample_base", "omni_collectionDate"=>"omni_collectionDate")
#     subset("sample_base"=> ByRow(!ismissing))
# end
    
# taxmd = DataFrame(metadata(taxa))
# leftjoin!(taxmd, select(reduced, "sample_base", "cogScore"), on=[:sample_base])
# set!(taxa, taxmd)

# taxmatch = taxa[:, .!ismissing.(get(taxa, :subject))]
# isKid = get(taxmatch, :isKid)
# isMom = get(taxmatch, :isMom)

# @info """
# ## For timpoints with sequenced stool sample

# - Total (N): $(nsamples(taxmatch))
# - Kids (N): $(count(isKid))
#     - Has age: $(count(!ismissing, get(taxmatch, :ageMonths)[isKid]))
#     - Has education: $(count(!ismissing, get(taxmatch, :education)[isKid]))
#     - Has race: $(count(r-> !ismissing(r) && r != "Unknown", get(taxmatch, :race)[isKid]))
#     - Has cogScore: $(count(r-> !ismissing(r) && r != "Unknown", get(taxmatch, :cogScore)[isKid]))
#     - Has all: $(count(map(r-> !ismissing(r) && r != "Unknown", get(taxmatch, :race)[isKid]) .& 
#                        map(!ismissing, get(taxmatch, :ageMonths)[isKid]) .&
#                        map(!ismissing, get(taxmatch, :education)[isKid]) .&
#                        map(!ismissing, get(taxmatch, :cogScore)[isKid])
#                 ))
# - Moms (N): $(count(isMom))
#     - Has education: $(count(!ismissing, get(taxmatch, :education)[isMom]))
#     - Has race: $(count(r-> !ismissing(r) && r != "Unknown", get(taxmatch, :race)[isMom]))
#     - Has all: $(count(map(r-> !ismissing(r) && r != "Unknown", get(taxmatch, :race)[isMom]) .& 
#                        map(!ismissing, get(taxmatch, :education)[isMom])
#                 ))
# """

# @info """
# ## Kids with all demographics (cogScore, education, race, sex)

# #- Exports

# @info "Subsetting Table"

# kids = subset(taxmatch, "isKid"=> identity)
# kids0120 = subset(kids, "ageMonths"=> ByRow(a-> !ismissing(a) && a < 120))
