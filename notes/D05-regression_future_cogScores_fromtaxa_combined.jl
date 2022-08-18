#####
# Notebook D05 - Regression of continuous future CogScores from current taxonomic data
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
using GLM
ml_rng = StableRNG(0)

#####
# Functions
#####

function train_randomforest_future_regressor(df; n_trials = 2, max_stool_ageMonths = 12.0, max_future_ageMonths=24.0, split_proportion=0.75, train_rng=Random.GLOBAL_RNG)

    # ## 1. Subsetting Data
    prediction_df = @chain df begin
        build_metadata_prediction_df(Symbol.(retained_featurenames), [ :cogScore ])
        select( Not(:cogScore) ) # Toggle this and the next line for dropping/filtering original cogScore
    #    subset(:cogScore => x -> .!(ismissing.(x)))
        subset(:ageMonths => x -> .!(ismissing.(x)))
        subset(:futureAgeMonths => x -> x .<= max_future_ageMonths)
        subset(:ageMonths => x -> x .<= max_stool_ageMonths)
        unique([:subject, :timepoint])
        unique([:subject, :futureTimepoint])
        dropmissing()
    end

    # ## 2. Separating inputs/outputs
    X = prediction_df[:, 9:end]
    y = prediction_df.futureCogScore

    # ## 4. Declaring hyperparameter tuning grid
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

    # ## 5. Initializing meta-arrays to record tuning results for each try

    trial_partitions = Vector{Tuple{Vector{Int64}, Vector{Int64}}}(undef, n_trials)
    trial_machines = Vector{Machine}(undef, n_trials)
    trial_slopecorrections = Vector{T where T <: RegressionModel}(undef, n_trials)
    trial_train_maes = repeat([Inf], n_trials)
    trial_train_mapes = repeat([Inf], n_trials)
    trial_train_cors = repeat([-1.0], n_trials)
    trial_test_maes = repeat([Inf], n_trials)
    trial_test_mapes = repeat([Inf], n_trials)
    trial_test_cors = repeat([-1.0], n_trials)

    # ## 6. Actual training loop

    @info "Performing $(n_trials) different train/test splits and tuning $(length(tuning_grid)) different hyperparmeter combinations\nfor the $(nrow(prediction_df)) samples with stool collection before $(max_stool_ageMonths) and next evaluation at $(max_future_ageMonths) months"

    for this_trial in 1:n_trials

        ## Splitting training data between train and test
        Random.seed!(train_rng, this_trial)
        train, test = partition(eachindex(1:nrow(X)), split_proportion, shuffle=true, rng=train_rng)
        trial_partitions[this_trial] = (train, test)

        @showprogress for i in 1:length(tuning_grid)
            Random.seed!(train_rng, 0)

            rf_model = RandomForestRegressor(
                min_samples_leaf = tuning_grid[i][1],
                max_depth = tuning_grid[i][2],
                sampling_fraction = tuning_grid[i][3],
                n_subfeatures = tuning_grid[i][4],
                n_trees = tuning_grid[i][5],
                rng=train_rng
            )

            rf_machine = machine(rf_model, X[train, :], y[train])
            MLJ.fit!(rf_machine, verbosity=0)

            train_y_hat = MLJ.predict(rf_machine, X[train, :])
            slope_correction_df = DataFrame([ train_y_hat, y[train] ], [:yhat, :y])
            slope_correction = lm(@formula(y ~ yhat), slope_correction_df)
            train_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => train_y_hat))
            train_mae = mae(train_y_hat, y[train])

            test_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => MLJ.predict(rf_machine, X[test, :])))
            test_mae = mae(test_y_hat, y[test])

            test_mae < trial_test_maes[this_trial] ? begin
                trial_train_maes[this_trial] = train_mae
                trial_train_mapes[this_trial] = mape(train_y_hat, y[train])
                trial_train_cors[this_trial] = Statistics.cor(train_y_hat, y[train])
                trial_test_maes[this_trial] = test_mae
                trial_test_mapes[this_trial] = mape(test_y_hat, y[test])
                trial_test_cors[this_trial] = Statistics.cor(test_y_hat, y[test])
                trial_machines[this_trial] = deepcopy(rf_machine)
            end : continue

        end # end for i in 1:length(tuning_grid)

    end # end for this_trial

    average_importances = DataFrame(
        :Variable => names(X),
        :Importance => map(mean, eachrow(reduce(hcat, [ impurity_importance(trial_machines[i].fitresult) for i in 1:n_trials ])))
        ); sort!(average_importances, :Importance, rev = true)

    # ## 7. Returning optimization results
    results = Dict(
        :input_data => prediction_df,
        :n_trials => n_trials,
        :selected_trial => findmax(trial_test_maes),
        :models => trial_machines,
        :dataset_partitions => trial_partitions,
        :slope_corrections => trial_slopecorrections,
        :train_maes => trial_train_maes,
        :train_mapes => trial_train_mapes,
        :train_correlations => trial_train_cors,
        :test_maes => trial_test_maes,
        :test_mapes => trial_test_mapes,
        :test_correlations => trial_test_cors,
        :importance_df => average_importances
    )

    @info "Done!"

    return results

end # end function

#####
# Computation Script
#####

include("D00-collect_taxonomic_cogscore_data.jl")

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree
regression_combined_results = train_randomforest_future_regressor(cogscore_taxa_df; n_trials = 10, split_proportion=0.75, train_rng=ml_rng)

using JLD2
JLD2.@save "models/results_regression_futureCogScores_combined_onlytaxa.jld" regression_combined_results

# merit_report_df = DataFrame( #TODO
#     :Iteration => collect(1:n_trials),
#     :Train_MAE => trial_train_mae,
#     :Test_MAE => trial_test_mae,
#     :Train_MAPE => trial_train_mape,
#     :Test_MAPE => trial_test_mape,
#     :Train_COR => trial_train_cor,
#     :Test_COR => trial_test_cor
# )