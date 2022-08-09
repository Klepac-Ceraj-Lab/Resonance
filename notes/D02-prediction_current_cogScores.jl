#####
# Notebook D02 - Prediction of futureCogScores from current taxonomic data
#####

using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using ProgressMeter
using MLJ
ml_rng = StableRNG(0)

# ## 1. Collecting the data to make the future dataframe
include("D00-collect_taxonomic_cogscore_data.jl")

max_prediction_ageMonths = 24

### 2.2. Subsetting data
prediction_df = @chain cogscore_taxa_df begin
    subset(:ageMonths => x -> .!(ismissing.(x)))
    subset(:cogScore => x -> .!(ismissing.(x)))
    subset(:sample => x -> .!(ismissing.(x)))
    dropmissing()
    subset(:ageMonths => x -> x .<= max_prediction_ageMonths)
    # unique( [:subject] )
    # DataFrames.select(Not(:timepointDelta))
end

### 2.4. Separating predictors and target from the additional metadata
y = prediction_df.cogScore
X = dropmissing(prediction_df[:, names(prediction_df) .∈ Ref(retained_featurenames)])
for m in models(matching(X,y)) println(m) end

## 3. Training the predictor model

### 3.1. Loading the desired model
RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree

#### Splitting training data between train and test
Random.seed!(ml_rng, 1)
train, test = partition(eachindex(1:nrow(X)), 0.6, shuffle=true, rng=ml_rng)

### 2.5. Removing microbiome columns with zero standard deviation
standard_deviations = [ std(X[train,i]) for i in 1:ncol(X) ]
is_zerostd_col = standard_deviations .== 0
println("$(sum(is_zerostd_col)) columns have zero standard deviation!")
zerostd_cols = findall(is_zerostd_col)
zerostd_colnames = names(X)[zerostd_cols]
X = DataFrames.select(X, Not(zerostd_colnames))

    # Declaring the optimization grid
    nodesize_range = collect(1:2:25)
    maxnodes_range = collect(1:2:25)
    sampsize_range = collect(10:20:110) ./ nrow(X[train, :])
    mtry_range = collect(5:2:30)
    n_trees_range = [ 500 ]

    tuning_grid = vec(collect(Base.product(
        nodesize_range,
        maxnodes_range,
        sampsize_range,
        mtry_range,
        n_trees_range
    )))

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
sort(train_cor_vector[.!(isnan.(train_cor_vector))])
println(maximum(test_cor_vector))
sort(test_cor_vector[.!(isnan.(test_cor_vector))])

#### Selecting model with the best performance on the test set

using GLM

selected_machine = bestmodel

train_y_hat = predict(selected_machine, X[train, :])
slope_correction_df = DataFrame([ train_y_hat, y[train] ], [:yhat, :y])
slope_correction = lm(@formula(y ~ yhat), slope_correction_df)
train_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => train_y_hat))
train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
train_cor = Statistics.cor(train_y_hat, y[train])

test_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => predict(selected_machine, X[test, :])))
test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
test_cor = Statistics.cor(test_y_hat, y[test])

#### Plotting and saving figures

using CairoMakie

train_set_plot = Figure(resolution = (800, 600))
ax = Axis(train_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Current cogScore for the *TRAIN* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[train], train_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(train_cor; digits = 2))"], [Point(60, 100)], textsize = 40)

train_set_plot

save("figures/prediction_currentcogscore_trainset_comparison.png", train_set_plot)

test_set_plot = Figure()
ax = Axis(test_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Future cogScore for the *TEST* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[test], test_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(test_cor; digits = 2))"], [Point(60, 95)], textsize = 40)

test_set_plot

save("figures/prediction_currentcogscore_testset_comparison.png", test_set_plot)

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