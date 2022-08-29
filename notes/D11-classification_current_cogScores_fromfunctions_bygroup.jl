#####
# Notebook D11 - Classification of binary current CogScores above/below average from current taxonomic data
#####

using Resonance
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
using JLD2
ml_rng = StableRNG(0)

#####
# Functions
#####

meanclass(x::Vector{T} where T <: Real) = x .>= mean(x)

function train_randomforest_current_classifier(df; n_trials = 2, min_age = 0.0, max_age=6.0, split_proportion=0.75, train_rng=Random.GLOBAL_RNG)

    # ## 1. Subsetting Data
    prediction_df = @chain df begin
        subset(:ageMonths => x -> x .>= min_age)
        subset(:ageMonths => x -> x .< max_age)    
    end
    
    # ## 2. Separating inputs/outputs
    X = prediction_df[:, 6:end]
    y = @chain prediction_df.cogScore begin
        meanclass
        coerce(OrderedFactor)
    end

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
    trial_train_accuracies = zeros(Float64, n_trials)
    trial_test_accuracies = zeros(Float64, n_trials)

    # ## 6. Actual training loop

    @info "Performing $(n_trials) different train/test splits and tuning $(length(tuning_grid)) different hyperparmeter combinations\nfor the $(nrow(prediction_df)) samples between $(min_age) and $(max_age) months"

    for this_trial in 1:n_trials

        ## Splitting training data between train and test
        Random.seed!(train_rng, this_trial)
        train, test = partition(eachindex(1:nrow(X)), split_proportion, shuffle=true, rng=train_rng)
        trial_partitions[this_trial] = (train, test)

        @showprogress for i in 1:length(tuning_grid)
            Random.seed!(train_rng, 0)

            rf_model = RandomForestClassifier(
                min_samples_leaf = tuning_grid[i][1],
                max_depth = tuning_grid[i][2],
                sampling_fraction = tuning_grid[i][3],
                n_subfeatures = tuning_grid[i][4],
                n_trees = tuning_grid[i][5],
                rng=train_rng
            )

            rf_machine = machine(rf_model, X[train, :], y[train])
            MLJ.fit!(rf_machine, verbosity=0)

            train_y_hat = MLJ.predict_mode(rf_machine, X[train, :]) 
            train_acc = mean(train_y_hat .== y[train])

            test_y_hat = MLJ.predict_mode(rf_machine, X[test, :]) 
            test_acc = mean(test_y_hat .== y[test])

            test_acc > trial_test_accuracies[this_trial] ? begin
                trial_test_accuracies[this_trial] = test_acc
                trial_train_accuracies[this_trial] = train_acc
                trial_machines[this_trial] = deepcopy(rf_machine)
            end : continue

        end # end for i in 1:length(tuning_grid)

    end # end for this_trial

    average_importances = DataFrame(
        :Variable => names(X),
        :Importance => map(mean, eachrow(reduce(hcat, [ impurity_importance(trial_machines[i].fitresult[1]) for i in 1:n_trials ])))
        ); sort!(average_importances, :Importance, rev = true)

    # ## 7. Returning optimization results
    results = Dict(
        :n_trials => n_trials,
        :inputs_outputs => (X,y),
        :selected_trial => findmax(trial_test_accuracies),
        :models => trial_machines,
        :dataset_partitions => trial_partitions,
        :train_accuracies => trial_train_accuracies,
        :test_accuracies => trial_test_accuracies,
        :importance_df => average_importances
    )

    @info "Done!"

    return results

end # end function

#####
# Computation Script
#####

## Getting data

mdata = @chain Resonance.load(Metadata()) begin
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :omni ] )
    dropmissing( [ :omni ] )
end

ecs = @chain Resonance.load(ECProfiles(); timepoint_metadata = mdata) begin
    filter(!hastaxon, _ )
end

ecs_df = DataFrame([eachcol(collect(ecs.abundances'))...], map( x -> string(x), ecs.features), copycols=true)
insertcols!(ecs_df, 1, :sample => collect(keys(ecs.sidx))[collect(values(ecs.sidx))])

cogscore_function_df = leftjoin(mdata, ecs_df, on = :omni => :sample, matchmissing=:error) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

max_prediction_ageMonths = 24
prediction_df = @chain cogscore_function_df begin
    dropmissing()
    subset(:ageMonths => x -> x .<= max_prediction_ageMonths)
end

## Training models and exporting results

RandomForestClassifier= MLJ.@load RandomForestClassifier pkg=DecisionTree

# classification_currentCogScores_00to06_fromfunctions_results = train_randomforest_current_classifier(prediction_df; n_trials = 5, min_age = 0.0, max_age=6.0, split_proportion=0.75, train_rng=ml_rng)
# JLD2.@save "models/classification_currentCogScores_00to06_fromfunctions_results.jld" classification_currentCogScores_00to06_fromfunctions_results
classification_currentCogScores_06to12_fromfunctions_results = train_randomforest_current_classifier(prediction_df; n_trials = 5, min_age = 6.0, max_age=12.0, split_proportion=0.75, train_rng=ml_rng)
JLD2.@save "models/classification_currentCogScores_06to12_fromfunctions_results.jld" classification_currentCogScores_06to12_fromfunctions_results
classification_currentCogScores_12to18_fromfunctions_results = train_randomforest_current_classifier(prediction_df; n_trials = 5, min_age = 12.0, max_age=18.0, split_proportion=0.75, train_rng=ml_rng)
JLD2.@save "models/classification_currentCogScores_12to18_fromfunctions_results.jld" classification_currentCogScores_12to18_fromfunctions_results
classification_currentCogScores_18to24_fromfunctions_results = train_randomforest_current_classifier(prediction_df; n_trials = 5, min_age = 18.0, max_age=24.0, split_proportion=0.75, train_rng=ml_rng)
JLD2.@save "models/classification_currentCogScores_18to24_fromfunctions_results.jld" classification_currentCogScores_18to24_fromfunctions_results