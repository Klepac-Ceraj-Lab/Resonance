#####
# Notebook D03 - Prediction of futureCogScores from current taxonomic data
#####

using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using ProgressMeter
using MLJ
using CairoMakie
using DecisionTree
using CairoMakie
using GLM
ml_rng = StableRNG(0)

#####
# Functions
#####

# ## 1. Collecting the data to make the future dataframe
include("D00-collect_taxonomic_cogscore_data.jl")

## 2. Building prediction inputs and outputs
prediction_df = build_metadata_prediction_df(cogscore_taxa_df, Symbol.(retained_featurenames), [ :cogScore ])

### 2.1. Data boundary conditions
max_original_age = 12
max_prediction_age = 24
max_ageMonthsDelta = 12
selected_timepoint_delta = 1

### 2.2. Subsetting data
prediction_df = @chain cogscore_taxa_df begin
    subset(:ageMonths => x -> .!(ismissing.(x)))
    subset(:cogScore => x -> .!(ismissing.(x)))
    subset(:sample => x -> .!(ismissing.(x)))
    dropmissing()
    subset(:ageMonths => x -> x .<= max_prediction_ageMonths)
end

X = prediction_df[:, 6:end]
y = prediction_df[:, 4]

for m in models(matching(X,y)) println(m) end

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree

# Declaring the optimization grid
nodesize_range = collect(1:1:15)
maxnodes_range = collect(1:1:20)
sampsize_range = [0.5, 0.6, 0.7, 0.8]#collect(5:2:35) ./ nrow(X)
mtry_range = collect(5:5:100)
ntrees_range = [100, 300, 500, 700]    

tuning_grid = vec(collect(Base.product(
    nodesize_range,
    maxnodes_range,
    sampsize_range,
    mtry_range,
    ntrees_range
)))

#### Tuning loop

#### Splitting training data between train and test
Random.seed!(ml_rng, 1)
train, test = partition(eachindex(1:nrow(X)), 0.75, shuffle=true, rng=ml_rng)

train_mae_vector = Vector{Float64}(undef, length(tuning_grid))
test_mae_vector = Vector{Float64}(undef, length(tuning_grid))

bestmae = +Inf

@showprogress for i in 1:length(tuning_grid)

    global bestmae
    Random.seed!(ml_rng, 0)

    rf_model = RandomForestRegressor(
        n_trees = tuning_grid[i][5],
        min_samples_leaf = tuning_grid[i][1],
        max_depth = tuning_grid[i][2],
        sampling_fraction = tuning_grid[i][3],
        n_subfeatures = tuning_grid[i][4],
        rng=ml_rng
    )

    rf_machine = machine(rf_model, X[train, :], y[train])
    MLJ.fit!(rf_machine, verbosity=0)

    train_y_hat = MLJ.predict(rf_machine, X[train, :])
    train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))

    test_y_hat = MLJ.predict_mode(rf_machine, X[test, :]) 
    test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))

    train_mae_vector[i] = train_mae
    test_mae_vector[i] = test_mae

    if !isnan(test_mae)
        if(test_mae) < bestmae
            global bestmae = test_mae
            global bestmodel = deepcopy(rf_machine)
        end
    end

end

println(bestmae)

selected_machine = bestmodel

train_y_hat = MLJ.predict(selected_machine, X[train, :])
slope_correction_df = DataFrame([ train_y_hat, y[train] ], [:yhat, :y])
slope_correction = lm(@formula(y ~ yhat), slope_correction_df)
train_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => train_y_hat))
train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
train_cor = Statistics.cor(train_y_hat, y[train])

test_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => MLJ.predict(selected_machine, X[test, :])))
test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
test_cor = Statistics.cor(test_y_hat, y[test])

train_set_plot = Figure(resolution = (800, 600))
ax = Axis(train_set_plot[1, 1],
    title = "Comparison of predicted vs ground truth *Current cogScore* for ages 18-24 months, *TRAIN* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[train], train_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(train_cor; digits = 2))"], [Point(60, 100)], textsize = 40)

train_set_plot

test_set_plot = Figure()
ax = Axis(test_set_plot[1, 1],
title = "Comparison of predicted vs ground truth *Current cogScore* for ages 18-24 months, *TEST* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[test], test_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(test_cor; digits = 2))"], [Point(60, 95)], textsize = 40)

test_set_plot

all_set_plot = Figure()
ax = Axis(all_set_plot[1, 1],
    title = "Comparison of predicted vs ground truth *Current cogScore* for ages 18-24 months,\n*FULL SAMPLE SET*",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)

yhat = GLM.predict(slope_correction, DataFrame( :yhat => MLJ.predict(selected_machine, X)))
sc_train = scatter!(ax, y[train], train_y_hat; color=:orange)
sc_test = scatter!(ax, y[test], test_y_hat; color=:blue)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 95)], textsize = 40)
Legend(all_set_plot[1,2], [sc_train, sc_test], ["Train", "Test"])

all_set_plot


## 2. Building prediction inputs and outputs
prediction_df = build_metadata_prediction_df(cogscore_taxa_df, Symbol.(retained_featurenames), [ :cogScore ])

### 2.1. Data boundary conditions
max_original_age = 12
max_prediction_age = 24
max_ageMonthsDelta = 12
selected_timepoint_delta = 1

### 2.2. Subsetting data
prediction_df = @chain prediction_df begin
    subset(:ageMonths => x -> .!(ismissing.(x)))
    subset(:ageMonths => x ->  x .<= max_original_age)
    subset(:futureAgeMonths => x ->  x .<= max_prediction_age)
    subset(:ageMonthsDelta => x ->  x .<= max_ageMonthsDelta)
    subset(:futureCogScore => x ->  x .>= 65.0)
    # unique( [:subject] )
    # DataFrames.select(Not(:timepointDelta))
end
prediction_df.futureCogScore = Float64.(prediction_df.futureCogScore)

### 2.3. Removing outliers with the Tietjen-Moore test
#n_outliers, outlier_idx = try_outliers(univariate_tietjenmoore, Float64.(prediction_df.futureCogScore), 10)
#println("Found $n_outliers outliers with the Tietjen-Moore test! Filtering DataFrame")
#println("New prediction table was reduced from $(nrow(predictionDf)) to $(nrow(predictionDf) - n_outliers) rows.")
#predictionDf = predictionDf[Not(outlier_idx), :]

### 2.4. Separating predictors and target from the additional metadata
y = prediction_df.futureCogScore
X = dropmissing(prediction_df[:, retained_featurenames])
for m in models(matching(X,y)) println(m) end

### 2.5. Removing microbiome columns with zero standard deviation
standard_deviations = [ std(X[:,i]) for i in 1:ncol(X) ]
is_zerostd_col = standard_deviations .== 0
println("$(sum(is_zerostd_col)) columns have zero standard deviation!")
zerostd_cols = findall(is_zerostd_col)
zerostd_colnames = names(X)[zerostd_cols]
X = DataFrames.select(X, Not(zerostd_colnames))

## 3. Training the predictor model

### 3.1. Loading the desired model

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree

### 3.2. Building the Tuning grid for hyperparameter tuning
# Hyperparameters of DecisionTree.RandomForestRegressor:
# ```
# n_trees = 10,                 Number of trees to train. Equivalent to R's randomForest `ntree`
# n_subfeatures = 0             Number of features to consider at random per split. Equivalent to R's randomForest `mtry`
# partial_sampling = 0.7        Fraction of samples to train each tree on. Equivalent to R's randomForest `sampsize`, but R takes an integer and this is a fraction.
# max_depth = -1,               Maximum depth of the decision trees. Equivalent to R's randomForest `maxnodes`
# min_samples_leaf = 1,         Equivalent to R's randomForest `nodesize`
# min_samples_split = 2,        No equivalence.
# min_purity_increase = 0.0,    No equivalence.
# ````

# range_n_trees = collect(50:25:500)
# range_n_subfeatures = collect(10:20:100)
# range_max_depth = collect(-1:2:20)
# range_partial_sampling = collect(0.5:0.1:0.9)
# range_min_samples_leaf = collect(5:2:20)
# range_min_purity_increase = [ 0.0, 0.1, 0.2, 0.3 ]

# # range_n_trees = collect(10:30:310)#collect(10:10:150)
# # range_n_subfeatures = collect(10:20:150)#collect(10:20:250)
# # range_partial_sampling = collect(0.3:0.1:0.7)
# # range_max_depth = [ -1 ]#collect(4:2:20)
# # range_min_samples_leaf = collect(1:1:5)#collect(1:1:10)
# # range_min_samples_split = [ 2, 3, 4 ]#[ 2, 3, 4, 5, 6, 7, 8 ]s
# # range_min_purity_increase = [ 0.0, 0.05, 0.10, 0.15 ]#[ 0.00, 0.05, 0.10, 0.15, 0.20, 0.25 ]

# tuning_grid = vec(collect(Base.product(
#     range_n_trees,
#     range_n_subfeatures,
#     range_max_depth,
#     range_partial_sampling,
#     range_min_samples_leaf,
#     range_min_purity_increase,
# )))

    # Declaring the optimization grid
    nodesize_range = collect(1:2:30)
    maxnodes_range = collect(1:2:30)
    sampsize_range = collect(5:2:80) ./ ncol(X)
    mtry_range = collect(5:2:60)
    n_trees_range = [ 500 ]

    tuning_grid = vec(collect(Base.product(
        nodesize_range,
        maxnodes_range,
        sampsize_range,
        mtry_range,
        n_trees_range
    )))

#### Tuning loop

#### Splitting training data between train and test
Random.seed!(ml_rng, 1)
train, test = partition(eachindex(1:nrow(X)), 0.5, shuffle=true, rng=ml_rng)

train_mae_vector = Vector{Float64}(undef, length(tuning_grid))
train_cor_vector = Vector{Float64}(undef, length(tuning_grid))
test_mae_vector = Vector{Float64}(undef, length(tuning_grid))
test_cor_vector = Vector{Float64}(undef, length(tuning_grid))

bestcor = 0.0

@showprogress for i in 1:length(tuning_grid)

    global bestcor
    Random.seed!(ml_rng, 0)

    rf_model = RandomForestRegressor(
        n_trees = tuning_grid[i][5],
        min_samples_leaf = tuning_grid[i][1],
        max_depth = tuning_grid[i][2],
        sampling_fraction = tuning_grid[i][3],
        n_subfeatures = tuning_grid[i][4],
        rng=ml_rng
    )

    # rf_model = RandomForestRegressor(
    #     n_trees = tuning_grid[i][1],
    #     n_subfeatures = tuning_grid[i][2],
    #     sampling_fraction = tuning_grid[i][4],
    #     max_depth = tuning_grid[i][3],
    #     min_samples_leaf = tuning_grid[i][5],
    #     min_purity_increase = tuning_grid[i][6],
    #     rng=ml_rng
    #     # n_trees = tuning_grid[i][1],
    #     # n_subfeatures = tuning_grid[i][2],
    #     # sampling_fraction = tuning_grid[i][3],
    #     # max_depth = tuning_grid[i][4],
    #     # min_samples_leaf = tuning_grid[i][5],
    #     # min_samples_split = tuning_grid[i][6],
    #     # min_purity_increase = tuning_grid[i][7],
    #     # rng=ml_rng
    # )

    rf_machine = machine(rf_model, X[train, :], y[train])
    fit!(rf_machine, verbosity=0)

    train_y_hat = predict(rf_machine, X[train, :]) 
    train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
    train_cor = Statistics.cor(train_y_hat, y[train])

    test_y_hat = predict(rf_machine, X[test, :])
    test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
    test_cor = Statistics.cor(test_y_hat, y[test])

#    machines_vector[i] = deepcopy(rf_machine)
    train_mae_vector[i] = train_mae
    train_cor_vector[i] = train_cor
    test_mae_vector[i] = test_mae
    test_cor_vector[i] = test_cor

    if !isnan(test_cor)
        if(test_cor) > bestcor
            global bestcor = test_cor
            global bestmodel = deepcopy(rf_machine)
        end
    end

end

println(minimum(train_mae_vector))
println(minimum(test_mae_vector))
println(maximum(train_cor_vector))
println(maximum(test_cor_vector))

#### Selecting model with the best performance on the test set

selected_machine = bestmodel

train_y_hat = predict(selected_machine, X[train, :]) 
train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
train_cor = Statistics.cor(train_y_hat, y[train])

test_y_hat = predict(selected_machine, X[test, :])
test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
test_cor = Statistics.cor(test_y_hat, y[test])

#### Plotting and saving figures

using CairoMakie

train_set_plot = Figure(resolution = (800, 600))
ax = Axis(train_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Future cogScore for the *TRAIN* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[train], train_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(train_cor; digits = 2))"], [Point(60, 100)], textsize = 40)

save("figures/prediction_futurecogscore_trainset_comparison.png", train_set_plot)

test_set_plot = Figure()
ax = Axis(test_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Future cogScore for the *TEST* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[test], test_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(test_cor; digits = 2))"], [Point(60, 95)], textsize = 40)

save("figures/prediction_futurecogscore_testset_comparison.png", test_set_plot)

using DecisionTree

@assert cor(
    apply_forest(selected_machine.fitresult, Matrix(X[test, :])),
    test_y_hat
) == 1.0

impurity_importances = impurity_importance(selected_machine.fitresult)
impurity_ordered_idx = sortperm(impurity_importances; rev = true)

impurityDf = DataFrame(
    :Variable => names(X)[impurity_ordered_idx],
    :Importance => impurity_importances[impurity_ordered_idx]
)

### Using linear regression to adjust future cogScore ~ futureAgeMonths and train on residuals

predictionDf[sortperm(predictionDf.target), :]
predictionDf[predictionDf.ageMonthsDelta .>= 12.0, :]

futurecogscore_futureagemonths_exploration = Figure()
ax = Axis(futurecogscore_futureagemonths_regression[1, 1],
    title = "futureCogScore as a function of futureAgeMonths",
    xlabel = "futureAgeMonths",
    ylabel = "futureCogScore"
)
scatter!(ax, predictionDf.futureAgeMonths, predictionDf.target)

futurecogscore_futureagemonths_exploration

using GLM

fm = @formula(target ~ futureAgeMonths)
linear_regressor = lm(fm, predictionDf)

intercept = linear_regressor.model.pp.beta0[1]
slope = linear_regressor.model.pp.beta0[2]

linear_prediction = GLM.predict(linearRegressor, hcat(ones(Float64, nrow(predictionDf)), predictionDf.futureAgeMonths))
residuals = predictionDf.target .- linear_prediction

futurecogscore_futureagemonths_regression = Figure()
ax = Axis(futurecogscore_futureagemonths_regression[1, 1],
    title = "true futureCogScore as a function of regression output from futureAgeMonths",
    xlabel = "futureCogScore",
    ylabel = "estimation of futureCogScore"
)
scatter!(ax, predictionDf.target, linear_prediction)

futurecogscore_futureagemonths_regression

####
predictionDf.target .= residuals