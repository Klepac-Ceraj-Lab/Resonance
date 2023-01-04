using Resonance
using Chain
using CairoMakie
using AlgebraOfGraphics

mdata = Resonance.load_raw_metadata()
taxa = Resonance.load_raw_metaphlan()
ecs = Resonance.load_raw_humann(; kind = "ecs", names=true, stratified=false)
kos = Resonance.load_raw_humann(; kind = "kos", names=true, stratified=false)
unirefs = Resonance.load_raw_humann(; kind = "genefamilies", names=false, stratified=false)
@assert samplenames(kos) == samplenames(ecs) == samplenames(taxa) == samplenames(unirefs)
mdata.isMom = map(tp-> !ismissing(tp) && tp != "-8" && contains(tp, r"^pre"i), mdata.ECHOTPCoded)
mdata.isKid = map(tp-> !ismissing(tp) && tp != "-8" && !contains(tp, r"^pre"i), mdata.ECHOTPCoded)

knead = Resonance.load(ReadCounts())
for p in [ecs, kos, unirefs]
    set!(p, select(knead, "sample_uid"=>"sample", AsTable(["final pair1", "final pair2"])=> ByRow(row-> row[1]+row[2]) =>"read_depth"))
end

#- Samples




# #- mgx / mbx

@info "samples with techreps:"
@info @chain samplenames(taxa) begin
    DataFrame(; sample=_)
    transform("sample"=> ByRow(s-> replace(s, r"_S\d+"=>"")) => "sample_base")
    groupby("sample_base")
    combine("sample"=> length=> "n_reps")
    subset("n_reps"=> ByRow(>(1)))
end

#-
@info """
  ## Raw Data

  - N rows in metadata (subject=> timepoint pairs): $(nrow(mdata))
  - N sequenced samples with taxaonomic profiles: $(nsamples(taxa))
  - N sequenced samples with functional profiles: $(nsamples(ecs))
  - N unique child stool samples: $(
      nrow(subset(mdata, "omni"=> ByRow(!ismissing), "isKid" => identity)))
  - N unique maternal stool samples: $(
      nrow(subset(mdata, "omni"=> ByRow(!ismissing), "isMom" => identity)))
  - N child subjects with at least 1 stool sample: $(
      nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity), "subject"))
  )
    - N child subjects with stool sample < 6 months: $(
      nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
                            "ageMonths" => ByRow(a-> !ismissing(a) && a < 6)), "subject"))

    )
    - N child subjects with stool sample < 8 months: $(
      nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
                            "ageMonths" => ByRow(a-> !ismissing(a) && a < 8)), "subject"))

    )
    - N child subjects with stool sample < 12 months: $(
      nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
                            "ageMonths" => ByRow(a-> !ismissing(a) && a < 12)), "subject"))

    )
    - N child subjects with stool sample > 18 months: $(
      nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isKid"=> identity, 
                            "ageMonths" => ByRow(a-> !ismissing(a) && a > 18)), "subject"))

    )
      - N maternal subjects with at least 1 stool sample: $(
      nrow(unique(subset(mdata, "omni"=> ByRow(!ismissing), "isMom"=> identity), "subject"))
  )

"""

#- Important Metadata (kids)

keep_fields = ["subject", "timepoint", "ageMonths", "education", "race", "cogScore", "sex",
               "isKid", "isMom", "omni", "omni_collectionDate", "etoh"]
reduced = @chain mdata begin
    select(keep_fields)
    rename("omni"=> "sample_base", "omni_collectionDate"=>"omni_collectionDate")
    subset("sample_base"=> ByRow(!ismissing))
end
    
taxmd = DataFrame(metadata(taxa))

leftjoin!(taxmd, reduced, on=[:sample_base])
for f in features(taxa)
    taxmd[!, string(f)] = vec(collect(abundances(taxa[f,:])))
end

subset!(taxmd, "subject"=> ByRow(!ismissing))
taxmd.isKid = Vector{Bool}(taxmd.isKid)
taxmd.isMom = Vector{Bool}(taxmd.isMom)
#- 

@info """
## For timpoints with sequenced stool sample

- Total (N): $(nrow(taxmd))
- Kids (N): $(count(taxmd.isKid))
    - Has age: $(count(!ismissing, taxmd.ageMonths[taxmd.isKid]))
    - Has education: $(count(!ismissing, taxmd.education[taxmd.isKid]))
    - Has race: $(count(r-> !ismissing(r) && r != "Unknown", taxmd.race[taxmd.isKid]))
    - Has cogScore: $(count(r-> !ismissing(r) && r != "Unknown", taxmd.cogScore[taxmd.isKid]))
    - Has all: $(count(map(r-> !ismissing(r) && r != "Unknown", taxmd.race[taxmd.isKid]) .& 
                       map(!ismissing, taxmd.ageMonths[taxmd.isKid]) .&
                       map(!ismissing, taxmd.education[taxmd.isKid]) .&
                       map(!ismissing, taxmd.cogScore[taxmd.isKid])
                ))
- Moms (N): $(count(isMom))
    - Has education: $(count(!ismissing, taxmd.education[taxmd.isMom]))
    - Has race: $(count(r-> !ismissing(r) && r != "Unknown", taxmd.race[taxmd.isMom]))
    - Has all: $(count(map(r-> !ismissing(r) && r != "Unknown", taxmd.race[taxmd.isMom]) .& 
                       map(!ismissing, taxmd.education[taxmd.isMom])
                ))
"""

@info """
## Kids with all demographics (cogScore, education, race, sex)
"""
#- Exports

@info "Subsetting Table"


taxakids = @chain taxmd begin
    subset("isKid"=> identity)
    sort(["subject", "timepoint"])
end

kids0120 = @chain taxakids begin
    subset("ageMonths" => ByRow(<=(120)), "cogScore"=> ByRow(!ismissing))
    sort(["subject", "timepoint"])
    groupby("subject")
    transform("subject" => (col -> [[true]; fill(false, length(col) - 1)])=> "first")
end

kids06 = @chain taxakids begin
    subset("ageMonths" => ByRow(<=(6)), "cogScore"=> ByRow(!ismissing))
    sort(["subject", "timepoint"])
    groupby("subject")
    transform("subject" => (col -> [[true]; fill(false, length(col) - 1)])=> "first")
end

kids18120 = @chain taxakids begin
    subset("ageMonths" => ByRow(a-> 18 < a <= 120), "cogScore"=> ByRow(!ismissing))
    sort(["subject", "timepoint"])
    groupby("subject")
    transform("subject" => (col -> [[true]; fill(false, length(col) - 1)])=> "first")
end

isdir(tablefiles()) || mkpath(tablefiles())
CSV.write(tablefiles("taxa_all_wide.csv"), select(taxmd, Cols(keep_fields[1:end-3]...), "sample_base", "omni_collectionDate" => "collectionDate", Cols(r"^s__")))
CSV.write(tablefiles("taxa_kids_all_wide.csv"), select(taxakids, Cols(keep_fields[1:end-3]...), "sample_base", "omni_collectionDate" => "collectionDate", Cols(r"^s__")))
CSV.write(tablefiles("taxa_kids_0to120.csv"), select(kids0120, Cols(keep_fields[1:end-3]...), "sample_base", "omni_collectionDate" => "collectionDate", Cols(r"^s__")))
CSV.write(tablefiles("taxa_kids_0to6.csv"), select(kids06, Cols(keep_fields[1:end-3]...), "sample_base", "omni_collectionDate" => "collectionDate", Cols(r"^s__")))
CSV.write(tablefiles("taxa_kids_18to120.csv"), select(kids18120, Cols(keep_fields[1:end-3]...), "sample_base", "omni_collectionDate" => "collectionDate", Cols(r"^s__")))

candace = parse.(Int, readlines("data/candace_ids.txt"))
CSV.write("data/candace_wide.csv", subset(taxmd, "subject"=> ByRow(s-> s in candace)))


## Uploads

filtcols = @chain kids0120 begin
    sort(["subject", "timepoint"])
    groupby("subject")
    transform("ageMonths" => (ages -> begin
                isfirst = fill(false, length(ages))
                i = findfirst(a-> true, ages)
                !isnothing(i) && (isfirst[i] = true)
                isfirst
            end) => "filter_00to120";
            ungroup = false
    )
    transform!("ageMonths" => (ages -> begin
                isfirst = fill(false, length(ages))
                i = findfirst(a-> a <= 6, ages)
                !isnothing(i) && (isfirst[i] = true)
                isfirst
            end) => "filter_00to06";
            ungroup = false
    )
    transform!("ageMonths" => (ages -> begin
                isfirst = fill(false, length(ages))
                i = findfirst(a-> 18 <= a <= 120, ages)
                !isnothing(i) && (isfirst[i] = true)
                isfirst
            end) => "filter_18to120"
    )
end


CSV.write(tablefiles("metadata_filtered.csv"), filtcols)
Resonance.write_arrow(tablefiles("taxa.arrow"), taxa[taxrank.(features(taxa)) .== :species, filtcols.sample])
Resonance.write_arrow(tablefiles("ecs.arrow"), ecs[:, filtcols.sample])
Resonance.write_arrow(tablefiles("kos.arrow"), kos[:, filtcols.sample])
Resonance.write_arrow(tablefiles("unirefs.arrow"), unirefs[:, filtcols.sample])
