# #Machine Learning Models

# ## Loading packages

using Resonance
using CairoMakie
using CSV
using DataFrames
using MLJ
using Statistics
using StableRNGs
ml_rng = StableRNG(102528)

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

retainfirst = true # if true, retain the first technical/biological replicate; if false, retain the last technical/biological replicate

if retainfirst

    unstackedDf = @chain rawDf begin
        groupby( [:subject, :timepoint, :variable] )
        combine(:value=>first)
        unstack(:variable, :value_first)
    end

else 
    unstackedDf = unstack(rawDf, :variable, :value; allowduplicates=true)
end

## Prediction workflows

### Problem 1: predict current cogScore from current taxonomic profile

#### Building predictionDf

select_columns = [ "cogScore", "ageMonths", microbiome_predictors... ]

##### 1. Removing rows with missing values

predictionDf = @chain unstackedDf begin
    select(select_columns)
    subset(:ageMonths => x -> .!(ismissing.(x)))
    subset(:cogScore => x -> .!(ismissing.(x)))
    subset(microbiome_predictors[1] => x -> .!(ismissing.(x)))
end

mapcols!(col -> parse.(Float64, col), predictionDf)

##### 2. Filtering by ageMonths < max_prediction_age

max_prediction_age = 120

predictionDf = @chain predictionDf begin
    subset(:ageMonths => x ->  x .<= max_prediction_age)
    DataFrames.transform(:ageMonths => ByRow(x -> x / max_prediction_age) => [:ageMonths])
end

##### 3. Removing microbiome columns with zero standard deviation

standard_deviations = [ std(predictionDf[:,i]) for i in 1:ncol(predictionDf) ]
is_zerostd_col = standard_deviations .== 0
is_negligiblestd_col = standard_deviations .<= 0.1
zerostd_cols = findall(is_zerostd_col)
negligiblestd_cols = findall(is_negligiblestd_col)
zerostd_colnames = names(predictionDf)[zerostd_cols]
negligiblestd_colnames = names(predictionDf)[negligiblestd_cols]

predictionDf = select(predictionDf, Not(negligiblestd_colnames))

#### Splitting training data between train and test

train, test = partition(eachindex(1:nrow(predictionDf)), 0.7, shuffle=true, rng=ml_rng)

#### Checking the schema and MLJ scitypes

schema(predictionDf)

#### Separating the predictors and target variable, vewing matching machine models

y, X = unpack(predictionDf, ==(:cogScore); rng=ml_rng);

for m in models(matching(X,y)) println(m) end

#### Building the MLJ machine for target problem.

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree
Standardizer = @load Standardizer pkg=MLJModels

input_standardizer = Standardizer()
input_standardization_machine = machine(input_standardizer, X)
fit!(input_standardization_machine, rows = train)
X_scaled = MLJ.transform(input_standardization_machine, X)

target_standardizer = Standardizer()
target_standardization_machine = machine(target_standardizer, y)
fit!(target_standardization_machine, rows = train)
y_scaled = MLJ.transform(target_standardization_machine, y)

rf_model = RandomForestRegressor(n_subfeatures = 0)

# Hyperparameters of DecisionTree.RandomForestRegressor:
# ```
# max_depth = -1, 
# min_samples_leaf = 1, 
# min_samples_split = 2, 
# min_purity_increase = 0.0, 
# n_subfeatures = -1, 
# n_trees = 10, 
# sampling_fraction = 0.7
# ````

#range_max_depth = range(rf_model, :max_depth; lower=-1, upper=5)
range_samples_leaf = range(rf_model, :min_samples_leaf; lower=1, upper=5)
range_min_samples_split = range(rf_model, :min_samples_split; lower=1, upper=5)
range_min_purity_increase = range(rf_model, :min_purity_increase; lower=0.0, upper=0.3)
range_n_trees = range(rf_model, :n_trees; lower=50, upper=450)
range_sampling_fraction = range(rf_model, :sampling_fraction; lower=0.5, upper=0.9)

ranges_vector = [
#    range_max_depth,
    range_samples_leaf,
#    range_min_samples_split,
#    range_min_purity_increase,
    range_n_trees,
#    range_sampling_fraction
]

tunable_rf_model = TunedModel(
    model = rf_model,
    ranges = ranges_vector,
    resampling = CV(nfolds=3),
    tuning = Grid(resolution=5),
    measure = mae,
    acceleration = CPUThreads()
)

tuning_rf_machine = machine(tunable_rf_model, X_scaled, y_scaled)

#### Training the Machine on `train` rows of `X`

fit!(tuning_rf_machine, rows=train)

tuning_report = report(tuning_rf_machine)
@show tuning_report.best_history_entry
@show tuning_report.best_model

best_rf_model = eval(tuning_report.best_model)

best_rf_machine = machine(best_rf_model, X_scaled, y_scaled)

#### Training the Machine on `train` rows of `X`
    
fit!(best_rf_machine, rows=train)    

#### Evaluating the Machine on `test` rows of `X`

train_y_hat = predict(best_rf_machine, rows=train)

test_y_hat = predict(best_rf_machine, rows=test)

train_abs_errors = MLJ.MLJBase.l1(train_y_hat, y_scaled[train]) # Absolute ErrorS
train_mean_abs_error = mean(train_abs_errors) # Mean Absolute Error
train_sq_errors = MLJ.MLJBase.l2(train_y_hat, y_scaled[train]) # Square ErrorS
train_mean_sq_error = mean(train_sq_errors) # Mean Square Error

test_abs_errors = MLJ.MLJBase.l1(test_y_hat, y_scaled[test]) # Absolute ErrorS
test_mean_abs_error = mean(test_abs_errors) # Mean Absolute Error
test_sq_errors = MLJ.MLJBase.l2(test_y_hat, y_scaled[test]) # Square ErrorS
test_mean_sq_error = mean(test_sq_errors) # Mean Square Error

Statistics.cor(y_scaled[train], train_y_hat)
Statistics.cor(y_scaled[test], test_y_hat)

scatter(y_scaled[train], train_y_hat)
scatter(y_scaled[test], test_y_hat)