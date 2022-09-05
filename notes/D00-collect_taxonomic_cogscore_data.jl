#####
# Notebook D00: Data collection and preparation notebook for analysis scripts
# NOTE: this notebook will become functions on src/prediction.jl source in the future.
#####

## Loading packages

using Resonance
using CSV
using Arrow
using Dates
using DataFrames

## Collecting the data from local files
resonance_datapath = get(ENV, "EXPORTED_DATA_PATH", "data")

### 1. Sample and subject metadata

timepoints_df = CSV.read(
    joinpath(resonance_datapath, "timepoint_metadata.csv"),
    DataFrame;
    delim = ',',
    types = [
        Int64,      #subject
        Int64,      #timepoint
        Float64,    #ageMonths
        String,     #race
        String,     #maternalEd
        Dates.Date, #assessmentDate
        Dates.Date, #scanDate
        Float64,    #cogScore
        String,     #omni
        String,     #etoh
        Bool,       #has_segmentation
        Float64,    #read_depth
        ]
)

#### Assessing timepoints sanity
# if any(.!ismissing.(timepoints_df.cogScore) .⊻ .!ismissing.(timepoints_df.assessmentDate))
#     xor_idx = findall(.!ismissing.(timepoints_df.cogScore) .⊻ .!ismissing.(timepoints_df.assessmentDate))
#     xor_filtered_df = timepoints_df[xor_idx, :]
#     @warn "There are $(nrow(xor_filtered_df)) timepoint samples with an assessmentDate incompatible with cogScore:\n $xor_filtered_df"
# end

### 2. Taxonomic profiles

taxonomic_data = DataFrame(Arrow.Table(joinpath(resonance_datapath, "taxa.arrow")))
taxonomic_profiles_matrix = zeros(Float64, maximum(taxonomic_data.sidx), maximum(taxonomic_data.fidx))

for (i,j,v) in zip(taxonomic_data.sidx, taxonomic_data.fidx, taxonomic_data.abundance)
    taxonomic_profiles_matrix[i,j] = v / 100.0
end

original_feature_names = readlines(joinpath(resonance_datapath, "taxa_features.txt"))
retained_features_idx = occursin.("s__", original_feature_names)
replaced_filtered_featurenames = replace.(original_feature_names[retained_features_idx], r"(.*?)s__" => "")
taxonomic_profiles_df = DataFrame(taxonomic_profiles_matrix[:, retained_features_idx], replaced_filtered_featurenames)[:, sortperm(replaced_filtered_featurenames)]
insertcols!(taxonomic_profiles_df, 1, :sample => unique(sort(taxonomic_data, :sidx).sample))
retained_featurenames = replaced_filtered_featurenames[sortperm(replaced_filtered_featurenames)]

##### The following block of code resolves the issue with "[Collinsella]_massiliensis" and "Collinsella_massiliensis". Some lines used during elaboration were retained and commented, but will be removed in the future.
#any( (taxonomic_profiles_df.Collinsella_massiliensis .!= 0.0) .& (taxonomic_profiles_df."[Collinsella]_massiliensis" .!= 0.0) )
#colinsella_conflict_idx = (taxonomic_profiles_df.Collinsella_massiliensis .!= 0.0) .& (taxonomic_profiles_df."[Collinsella]_massiliensis" .!= 0.0)
#taxonomic_profiles_df[colinsella_conflict_idx, [ "Collinsella_massiliensis", "[Collinsella]_massiliensis"] ]
#taxonomic_profiles_df.Collinsella_massiliensis .= taxonomic_profiles_df.Collinsella_massiliensis .+ taxonomic_profiles_df."[Collinsella]_massiliensis"
select!(taxonomic_profiles_df, Not("[Collinsella]_massiliensis"))
retained_featurenames = retained_featurenames[retained_featurenames .!= "[Collinsella]_massiliensis"]

#### Assessing taxonomic profiles sanity
if any(.!(timepoints_df.omni[.!(ismissing.(timepoints_df.omni))] .∈ Ref(taxonomic_profiles_df.sample)))
    @warn "There are timepoint samples with a omni sample incompatible with the taxonomic profile data!"
end

### Joining the resulting taxonomic profiles and sample metadata.

cogscore_taxa_df = leftjoin(timepoints_df, taxonomic_profiles_df, on = :omni => :sample, validate = (false, true), matchmissing = :equal) |>
    y -> rename!(y, :omni => :sample) |>
    y -> sort(y, [ :subject, :timepoint ]) |>
    y -> select!(y, [ :subject, :timepoint, :ageMonths, :cogScore, :sample, Symbol.(retained_featurenames)... ])

CSV.write("cogscore_taxa_df.csv", cogscore_taxa_df)


#### Data investigation algorithms

## 1. histogram of first cogscore for every subject

# nonmissing_cogscore = @chain cogscore_taxa_df begin
#     dropmissing(:cogScore)
# #    subset(:sample=> ByRow(!ismissing))
#     sort([:subject, :timepoint])
# end

# unique_first_cogScore = @chain nonmissing_cogscore begin
#     unique(:subject)
# end

# n_outliers, outlier_idx = try_outliers(univariate_tietjenmoore, unique_first_cogScore.cogScore, 30; reverse=false)
# hist(unique_first_cogScore.cogScore, bins = 10)
# unique_first_cogScore[outlier_idx, :]