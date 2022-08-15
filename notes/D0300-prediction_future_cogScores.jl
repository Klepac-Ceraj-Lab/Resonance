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

### 2.1. Data boundary conditions
max_original_ageMonths = 12
max_prediction_ageMonths = 24
max_ageMonthsDelta = 12
#selected_timepoint_delta = 1 # Removed in favor of unique([:subject, :timepoint, :futureTimepoint])

prediction_df = @chain cogscore_taxa_df begin
    build_metadata_prediction_df(Symbol.(retained_featurenames), [ :cogScore ])
    select( Not(:cogScore) ) # Toggle this and the next line for dropping/filtering original cogScore
#    subset(:cogScore => x -> .!(ismissing.(x)))
    subset(:ageMonths => x -> .!(ismissing.(x)))
    subset(:futureAgeMonths => x -> x .<= max_prediction_ageMonths)
    subset(:ageMonths => x -> x .<= max_prediction_ageMonths)
    unique([:subject, :timepoint])
    unique([:subject, :futureTimepoint])
    dropmissing()
end

### 2.2. Subsetting data

X = prediction_df[:, 9:end] # Only taxonomic profile as input
# X = prediction_df[:, 9:end] # Taxonmic profile + original ageMonths
# X = prediction_df[:, 9:end] # Taxonmic profile + future ageMonths
# X = prediction_df[:, 9:end] # Taxonmic profile + original ageMonths + future ageMonths
y = prediction_df[:, 8]

for m in models(matching(X,y)) println(m) end

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree

# Declaring the optimization grid
nodesize_range = collect(1:2:15)
maxnodes_range = collect(1:2:25)
sampsize_range = [0.4, 0.5, 0.6, 0.7, 0.8]#collect(5:2:35) ./ nrow(X)
mtry_range = collect(10:10:100)
ntrees_range = [100, 300, 500, 700]

tuning_grid = vec(collect(Base.product(
    nodesize_range,
    maxnodes_range,
    sampsize_range,
    mtry_range,
    ntrees_range
)))

#### Iterative trials for train/test split and hyperparemeter tuning
n_trials = 10
trial_partition = Vector{Tuple{Vector{Int64}, Vector{Int64}}}
trial_machine = Vector{Machine}(undef, n_trials)
trial_slopecorrection = Vector{T where T <: RegressionModel}(undef, n_trials)
trial_plot = Vector{Figure}(undef, n_trials)
trial_train_mae = Vector{Float64}(undef, n_trials)
trial_train_mape = Vector{Float64}(undef, n_trials)
trial_train_cor = Vector{Float64}(undef, n_trials)
trial_test_mae = Vector{Float64}(undef, n_trials)
trial_test_mape = Vector{Float64}(undef, n_trials)
trial_test_cor = Vector{Float64}(undef, n_trials)

for split_trial in 1:n_trials

    #### Trial data partition between train and test
    Random.seed!(ml_rng, split_trial)
    train, test = partition(eachindex(1:nrow(X)), 0.75, shuffle=true, rng=ml_rng) # original:0.75

    train_mae_vector = Vector{Float64}(undef, length(tuning_grid))
    test_mae_vector = Vector{Float64}(undef, length(tuning_grid))

    bestmae = +Inf
    bestmodel = RandomForestRegressor()

    @showprogress for i in 1:length(tuning_grid)

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
                bestmae = test_mae
                bestmodel = deepcopy(rf_machine)
            end
        end

    end

    ## Calculateing Figures of Merit

    selected_machine = bestmodel

    train_y_hat = MLJ.predict(selected_machine, X[train, :])
    slope_correction_df = DataFrame([ train_y_hat, y[train] ], [:yhat, :y])
    slope_correction = lm(@formula(y ~ yhat), slope_correction_df)
    train_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => train_y_hat))
    train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
    train_mape = mape(train_y_hat, y[train])
    train_cor = Statistics.cor(train_y_hat, y[train])

    test_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => MLJ.predict(selected_machine, X[test, :])))
    test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
    test_mape = mape(test_y_hat, y[test])
    test_cor = Statistics.cor(test_y_hat, y[test])

    ## Bulding Plots

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

    # all_set_plot

    ## Collecting results
    
    trial_partition = (train, test)
    trial_machine[split_trial] = deepcopy(bestmodel)
    trial_slopecorrection[split_trial] = deepcopy(slope_correction)
    trial_plot[split_trial] = all_set_plot
    trial_train_mae[split_trial] = train_mae
    trial_train_mape[split_trial] = train_mape
    trial_train_cor[split_trial] = train_cor
    trial_test_mae[split_trial] = test_mae
    trial_test_mape[split_trial] = test_mape
    trial_test_cor[split_trial] = test_cor

end # end for split_trial

merit_report_df = DataFrame(
    :Iteration => collect(1:n_trials),
    :Train_MAE => trial_train_mae,
    :Test_MAE => trial_test_mae,
    :Train_MAPE => trial_train_mape,
    :Test_MAPE => trial_test_mape,
    :Train_COR => trial_train_cor,
    :Test_COR => trial_test_cor
)

# trial_partition
# trial_machine
# trial_slopecorrection
# trial_plot
# trial_train_mae
# trial_train_mape
# trial_train_cor
# trial_test_mae
# trial_test_mape
# trial_test_cor

#####
# Reporting Importances
#####

average_importances = DataFrame(
    :Variable => names(X),
    :Importance => map(mean, eachrow(reduce(hcat, [ impurity_importance(trial_machine[i].fitresult) for i in 1:n_trials ])))
); sort(average_importances, :Importance, rev = true)[1:30, :]

using JLD2

modeling_export = Dict(
        :n_trials => n_trials,
        :selected_trial => 1,
        :models => trial_machine,
        :dataset_partitions => trial_partition,
        :slope_corrections => trial_slopecorrection,
        :train_maes => trial_train_mae,
        :train_mapes => trial_train_mape,
        :train_correlations => trial_train_cor,
        :test_maes => trial_test_mae,
        :test_mapes => trial_test_mape,
        :test_correlations => trial_test_cor,
        :importance_df => average_importances
)

JLD2.@save "models/results_regression_futureCogScores_allAgeMonths_onlytaxa.jld" modeling_export

a = JLD2.load("models/results_regression_futureCogScores_allAgeMonths_onlytaxa.jld", "modeling_export")