# #Machine Learning Models

#####
# Loading packages
#####

using ProgressMeter
using CairoMakie
using CSV
using DataFrames
using MLJ
using Distributions
using Statistics
using StableRNGs
using Resonance

#####
# Functions to use
#####

#####
# Script
#####

# ## Parsing tidy CSV data

rawDf = CSV.read(
    joinpath("data", "tidy_timepoints_with_brain.csv"),
    DataFrame;
    delim = ',',
    types = [Int64, Int64, String, String],
    header = ["subject", "timepoint", "variable", "value"],
    skipto = 2
)

# ### Fixing some metadata names (see source for details)

#fix_metadata_colnames!(rawDf)

# ## Removing duplicates

# ### Removing actual duplicated lines
check_longdata_metaduplicates!(rawDf; remove_duplicates=true)

# ### Removing technical/biological replicates
retainfirst = parse(Bool, get(ENV, "RETAIN_FIRST", "true")) # if true, retain the first technical/biological replicate; if false, retain the last technical/biological replicate
unstackedDf = unstack_techreplicates(rawDf, retainfirst)

## Prediction workflows

### Problem 2: predict future cogScore from current taxonomic profile

#### Building predictionDf

select_columns = [ "target", "ageMonths", "futureAgeMonths", "ageMonthsDelta", "timepointDelta", microbiome_predictors... ]

##### 1. Removing rows with missing values and applying filters

max_input_age = 12
max_prediction_age = 24
max_ageMonthsDelta = 12
selected_timepoint_delta = 1

predictionDf = @chain unstackedDf begin
    subset(:ageMonths => x -> .!(ismissing.(x)))
    build_future_df(:cogScore)
    DataFrames.select(select_columns)
end

mapcols!(col -> tryparsecol(Float64, col), predictionDf)

##### 2. Filtering by ageMonths < max_prediction_age

predictionDf = @chain predictionDf begin
    subset(:ageMonths => x ->  x .<= max_input_age)
    subset(:futureAgeMonths => x ->  x .<= max_prediction_age)
    subset(:ageMonthsDelta => x ->  x .<= max_ageMonthsDelta)
    subset(:timepointDelta => x ->  x .== selected_timepoint_delta)
    DataFrames.select(Not(:timepointDelta))
#    DataFrames.transform(:ageMonths => ByRow(x -> x / max_prediction_age) => [:ageMonths])
end

##### Removing outliers with the Tietjen-Moore test
n_outliers, outlier_idx = try_outliers(univariate_tietjenmoore, predictionDf.target, 10)
println("Found $n_outliers outliers with the Tietjen-Moore test! Filtering DataFrame")
println("New prediction table was reduced from $(nrow(predictionDf)) to $(nrow(predictionDf) - n_outliers) rows.")
predictionDf = predictionDf[Not(outlier_idx), :]

##### 3. Removing microbiome columns with zero standard deviation

std_cutoff = 0.1
standard_deviations = [ std(predictionDf[:,i]) for i in 1:ncol(predictionDf) ]
is_zerostd_col = standard_deviations .== 0
is_negligiblestd_col = standard_deviations .<= std_cutoff
println("$(sum(is_zerostd_col)) columns have zero standard deviation!")
println("$(sum(is_negligiblestd_col)) columns have a standard deviation lower than $std_cutoff")
zerostd_cols = findall(is_zerostd_col)
negligiblestd_cols = findall(is_negligiblestd_col)
zerostd_colnames = names(predictionDf)[zerostd_cols]
negligiblestd_colnames = names(predictionDf)[negligiblestd_cols]

# Uncomment to remove zerostd columns
predictionDf = DataFrames.select(predictionDf, Not(zerostd_colnames))

#### Checking the schema and MLJ scitypes

schema(predictionDf)