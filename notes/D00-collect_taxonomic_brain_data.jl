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

#####
# Functions
#####

function map_string_timepoint(concatstring)::Tuple{Int64, Int64}
    letter_to_timepoint = Dict(
        'a' => 1, 'b' => 2, 'c' => 3, 'd' => 4,
        'e' => 5, 'f' => 6, 'g' => 7, 'h' => 8,
        'i' => 9, 'j' => 10,'k' => 11,'l' => 12)

    if isletter(concatstring[end])
        return parse(Int64, concatstring[1:end-1]), letter_to_timepoint[concatstring[end]]
    else
        return parse(Int64, concatstring), 1
    end # endif
end # end function

function split_subject_timepoint(stringsvector)
    mappings = map(map_string_timepoint, stringsvector)
    subjects = map(x -> x[1], mappings)
    timepoints = map(x -> x[2], mappings)
    return subjects, timepoints
end # end function

## Collecting the data from local files
resonance_datapath = get(ENV, "EXPORTED_DATA_PATH", "data")

### 1. Sample and subject metadata

timepoints_df = CSV.read(
    "data/timepoints.csv",
    #joinpath(resonance_datapath, "timepoint_metadata.csv"),
    DataFrame)

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

### 3. Brain data

brain_data = dropmissing(CSV.read(joinpath(resonance_datapath, "brain_measures.csv"), DataFrame))

#### 3.1. Brain data sanity checks
# rowsums = [ sum(brain_data[i, 3:end]) for i in 1:nrow(brain_data)]

brain_variables = names(brain_data)[3:end] # 1 is subject, 2 is timepoint, 3:end are brain values

### Joining the resulting taxonomic profiles and sample metadata.

metadata_taxa_df = innerjoin(timepoints_df, taxonomic_profiles_df, on = :omni => :sample, validate = (false, true), matchmissing = :notequal) |>
    y -> rename!(y, :omni => :sample) |>
    y -> sort(y, [ :subject, :timepoint ]) |>
    y -> select!(y, [ :subject, :timepoint, :ageMonths, :sample, Symbol.(retained_featurenames)... ])
    # y -> select!(y, [ :subject, :timepoint, :ageMonths, :cogScore, :sample, Symbol.(retained_featurenames)... ])

brain_taxa_df = innerjoin(metadata_taxa_df, brain_data, on = [:subject => :subject, :timepoint => :timepoint])

CSV.write("brain_taxa_df.csv", brain_taxa_df)

### 4. Raw brain data

raw_brain_data = CSV.read("/home/guilherme/Downloads/segmentationVolumeMeasurements_oct2021.csv", DataFrame)

raw_brain_subjects, raw_brain_timpoints = split_subject_timepoint(raw_brain_data."Measure:volume")

brain_targets = [
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
    "CSF",
    "IntraCranialVol"
]

raw_brain_data."Left-Thalamus" = raw_brain_data."Left-Thalamus" .+ raw_brain_data."Left-Thalamus-Proper"
raw_brain_data."Right-Thalamus" = raw_brain_data."Right-Thalamus" .+ raw_brain_data."Right-Thalamus-Proper"
raw_brain_data."IntraCranialVol" = raw_brain_data."EstimatedTotalIntraCranialVol" .+ raw_brain_data."IntraCranialVol"

brain_data = sort( 
    hcat(
        DataFrame(:subject => raw_brain_subjects, :timepoint => raw_brain_timpoints ),
        select(raw_brain_data, brain_targets)
        )
    , [:subject, :timepoint]
)

brain_taxa_df = innerjoin(timepoints_df, brain_data, on = [:subject => :subject, :timepoint => :timepoint])

#### Data investigation algorithms

## 1. histograms of each Brain variable

brain_data_matrix = Matrix(brain_data)[:, 3:end]

using Statistics

cor(brain_data_matrix; dims=1)

using CairoMakie

corrgram_figure = Figure(resolution = (900, 800))
ax = Axis(corrgram_figure[1, 1],
    title = "Correlation matrix of brain features",
    titlesize = 24.0f0,

    xticks = (1:length(brain_variables), brain_variables),
    xminorticksvisible = false,
    xticklabelrotation = -pi/2,
    xticklabelsize = 12,

    yticks = (1:length(brain_variables), brain_variables),
    yminorticksvisible = false,
    yticklabelsize = 12,
    yaxisposition = :right,
    yreversed = true
)

hm = heatmap!(ax, cor(brain_data_matrix; dims=1);colorrange = (-1.0, 1.0), colormap=["red", "white", "blue"])
Colorbar(corrgram_figure[:, end+1], hm;)

corrgram_figure