# #Machine Learning Models

#####
# Loading packages
#####

using ProgressMeter
using CairoMakie
using CSV
using DataFrames
using MLJ
using Statistics
using StableRNGs
using Resonance
ml_rng = StableRNG(102528)

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

### Problem 2: predict future cogScore from current taxonomic profile

#### Building predictionDf

select_columns = [ "target", "ageMonths", "futureAgeMonths", microbiome_predictors... ]

##### 1. Removing rows with missing values

predictionDf = @chain unstackedDf begin
    subset(:ageMonths => x -> .!(ismissing.(x)))
    build_future_df(:cogScore)
    DataFrames.select(select_columns)
end

function parseif(Type::DataType, col)
    if typeof(col) == Vector{Type}
        return(col)
    else
        return(parse.(Type, col))
    end
end

mapcols!(col -> parseif(Float64, col), predictionDf)

##### 2. Filtering by ageMonths < max_prediction_age

max_input_age = 120
max_prediction_age = 240

predictionDf = @chain predictionDf begin
    subset(:ageMonths => x ->  x .<= max_input_age)
    subset(:futureAgeMonths => x ->  x .<= max_prediction_age)
#    DataFrames.transform(:ageMonths => ByRow(x -> x / max_prediction_age) => [:ageMonths])
end

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

# Uncomment to remove negligible columns
#predictionDf = DataFrames.select(predictionDf, Not(negligiblestd_colnames))

#### Checking the schema and MLJ scitypes

schema(predictionDf)

#### Separating the predictors and target variable, vewing matching machine models

y, X = unpack(predictionDf, ==(:target));
for m in models(matching(X,y)) println(m) end

#### Computing the standardizer parameters

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree

#### Building the Tuning grid for training the Random Forest model
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

grid_resolution = 5
range_max_depth = [ -1 ] # range(-1, 3, length = grid_resolution)
range_samples_leaf = range(1, 5, length = grid_resolution)
range_min_purity_increase = [ 0.0 ] # range(0.0, 0.3, length = grid_resolution)
range_n_trees = range(50, 500, length = 2*grid_resolution)
range_sampling_fraction = range(0.5, 0.9, length = grid_resolution)

tuning_grid = vec(collect(Base.product(
    range_max_depth,
    range_samples_leaf,
    range_min_purity_increase,
    range_n_trees,
    range_sampling_fraction,
)))

#### Tuning loop

#### Splitting training data between train and test
train, test = partition(eachindex(1:nrow(predictionDf)), 0.7, shuffle=true, rng=ml_rng)

machines_vector = Vector{Machine{MLJDecisionTreeInterface.RandomForestRegressor, true}}(undef, length(tuning_grid))
train_mae_vector = Vector{Float64}(undef, length(tuning_grid))
train_cor_vector = Vector{Float64}(undef, length(tuning_grid))
test_mae_vector = Vector{Float64}(undef, length(tuning_grid))
test_cor_vector = Vector{Float64}(undef, length(tuning_grid))

@showprogress for i in 1:length(tuning_grid)

    rf_model = RandomForestRegressor(
        max_depth = tuning_grid[i][1], 
        min_samples_leaf = tuning_grid[i][2], 
        min_purity_increase = tuning_grid[i][3], 
        n_subfeatures = 0, 
        n_trees = tuning_grid[i][4], 
        sampling_fraction = tuning_grid[i][5],
        rng=ml_rng
    )

    rf_machine = machine(rf_model, X[train, :], y[train])
    fit!(rf_machine, verbosity=0)

    train_y_hat = predict(rf_machine, X[train, :]) 
    train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
    train_cor = Statistics.cor(train_y_hat, y[train])

    test_y_hat = predict(rf_machine, X[test, :])
    test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
    test_cor = Statistics.cor(test_y_hat, y[test])

    machines_vector[i] = deepcopy(rf_machine)
    train_mae_vector[i] = train_mae
    train_cor_vector[i] = train_cor
    test_mae_vector[i] = test_mae
    test_cor_vector[i] = test_cor

end

println(minimum(train_mae_vector))
println(minimum(test_mae_vector))
println(maximum(train_cor_vector))
println(maximum(test_cor_vector))

#### Selecting model with the best performance on the test set

selected_machine = machines_vector[findmax(test_cor_vector)[2]]

train_y_hat = predict(selected_machine, X[train, :]) 
train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
train_cor = Statistics.cor(train_y_hat, y[train])

test_y_hat = predict(selected_machine, X[test, :])
test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
test_cor = Statistics.cor(test_y_hat, y[test])

#### Plotting and saving figures

train_set_plot = Figure(resolution = (800, 600))
ax = Axis(train_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Future cogScore for the *TRAIN* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[train], train_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(train_cor; digits = 2))"], [Point(60, 120)], textsize = 40)

save("figures/prediction_futurecogscore_trainset_comparison.png", train_set_plot)

test_set_plot = Figure()
ax = Axis(test_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Future cogScore for the *TEST* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[test], test_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(test_cor; digits = 2))"], [Point(60, 120)], textsize = 40)

save("figures/prediction_futurecogscore_testset_comparison.png", test_set_plot)

using DecisionTree

@assert cor(
    apply_forest(selected_machine.fitresult, Matrix(X[test, :])),
    test_y_hat
) == 1.0

impurity_importances = impurity_importance(selected_machine.fitresult)
impurity_ordered_idx = sortperm(impurity_importances; rev = true)

impurityDF = DataFrame(
    :Variable => names(X)[impurity_ordered_idx],
    :Importance => impurity_importances[impurity_ordered_idx]
)