##########
# Concurrent cognitive assessment Score regression notebook
##########

#####
# 0. Loading libraries
#####

using Distributed
addprocs(8)  # 8 procs.

@everywhere begin
    using Chain
    using CSV
    using Statistics
    using Random
    using DataFrames
    using StableRNGs
    using MLJ
    using CairoMakie
    using DecisionTree
    using JLD2
    using CategoricalArrays
    using Resonance
    using GLM
    using ShapML        
    ml_rng = StableRNG(0)
end

#####
# Loading Data
#####

## 1. Metadata
mdata = Resonance.load(Metadata())
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 
transform!(seqs, "sample"=> ByRow(String)=> "sample")
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)

## 2. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
species = filter(f-> taxrank(f) == :species, taxa)
mdata_taxa_df = sort(Resonance.comm2wide(species), [ :subject, :timepoint ]);

## 3. Functional Profiles
ecs = Resonance.load(ECProfiles(); timepoint_metadata = seqs)
mdata_ecs_df = sort(Resonance.comm2wide(ecs), [ :subject, :timepoint ]);

### 3.1. Hashing EC names with fixed-length alphabetic hashes due to length and special characters; storing text files for remapping.
oldnames = names(mdata_ecs_df)[20:end]
for oldname in names(mdata_ecs_df)[20:end]
    rename!(mdata_ecs_df, oldname => randstring(['A':'Z'; 'a':'z'], 12))
end
newnames = names(mdata_ecs_df)[20:end]
open(scratchfiles("longnames.txt"), "w") do io
    for i in oldnames
        println(io, i)
    end
end
open(scratchfiles("hashnames.txt"), "w") do io
    for i in newnames
        println(io, i)
    end
end

#####
# Applying prevalence filters
#####
taxa_prevalence_threshold_lower = 0.10
taxa_prevalence_threshold_upper = 1.00
ecs_prevalence_threshold_lower = 0.25
ecs_prevalence_threshold_upper = 0.85

## 00 to 06 months

filtered_mdata_taxa_df_00to06 = filter_prevalences(
    # dropmissing(mdata_taxa_df[mdata_taxa_df.filter_00to06, :], Not([:race, :read_depth])),
    dropmissing(
        mdata_taxa_df[mdata_taxa_df.filter_00to06, :],
        Not([
            "sample",
            "read_depth",
            "edfloat",
            "race",
            "filter_00to120",
            "filter_00to06",
            "filter_18to120",
            "filter_future_omni",
            "filter_future_cog",
            "omni",
            "Mullen::mullen_EarlyLearningComposite",
            "Mullen::mullen_NonVerbalComposite",
            "Mullen::mullen_VerbalComposite"]
        )),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"]);
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

# filtered_mdata_ecs_df_00to06 = filter_prevalences(
#     dropmissing(mdata_ecs_df[mdata_ecs_df.filter_00to06, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"])),
#     :cogScore,
#     [:subject, :timepoint, :sex, :education, :ageMonths],
#     Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"]);
#     lbound = ecs_prevalence_threshold_lower,
#     ubound = ecs_prevalence_threshold_upper
# )

filtered_mdata_taxa_df_18to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_18to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"]);
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

# filtered_mdata_ecs_df_18to120 = filter_prevalences(
#     dropmissing(mdata_ecs_df[mdata_ecs_df.filter_18to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"])),
#     :cogScore,
#     [:subject, :timepoint, :sex, :education, :ageMonths],
#     Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"]);
#     lbound = ecs_prevalence_threshold_lower,
#     ubound = ecs_prevalence_threshold_upper
# )

filtered_mdata_taxa_df_00to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_00to120, :], Not(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"])),
    :cogScore,
    [:subject, :timepoint, :sex, :education, :ageMonths],
    Symbol.(["sample", "read_depth", "edfloat", "race", "filter_00to120", "filter_00to06", "filter_18to120", "filter_future_omni", "filter_future_cog", "omni", "Mullen::mullen_EarlyLearningComposite", "Mullen::mullen_NonVerbalComposite", "Mullen::mullen_VerbalComposite"]);
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

# ## Exports

isdir(tablefiles("figure3")) || mkdir(tablefiles("figure3"))

CSV.write(tablefiles("figure3", "00to06_taxa.csv"), filtered_mdata_taxa_df_00to06) # @Hugemiler what are these used for?
CSV.write(tablefiles("figure3", "00to06_ecs.csv"), filtered_mdata_ecs_df_00to06)
CSV.write(tablefiles("figure3", "18to120_taxa.csv"), filtered_mdata_taxa_df_18to120)
CSV.write(tablefiles("figure3", "18to120_ecs.csv"), filtered_mdata_ecs_df_18to120)
CSV.write(tablefiles("figure3", "00to120_taxa.csv"), filtered_mdata_taxa_df_00to120)

#####
# Training models
#####

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

onlydemo_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

taxa_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

ecs_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

#####
# 00 to 06 months
#####

function new_randomforest(
    name::String,
    original_df::DataFrame,
    data_preparation_fun::Function,
    input_cols,
    output_col;
    unique_col = :subject,
    n_folds = 10,
    n_replicas = 10,
    n_rngs = 100,
    tuning_space = (;
        maxnodes_range = [ -1 ],
        nodesize_range = [ 5 ],
        sampsize_range = [ 0.7 ],
        mtry_range = [ -1 ],
        ntrees_range = [ 10 ]
    ),
    train_rng=Random.GLOBAL_RNG
    )

    # ## 1. Subsetting Data and shuffling DataFrame
    prediction_df = data_preparation_fun(original_df)
    input_cols = collect(input_cols[1]:ncol(prediction_df))

    print("prediction_df initialized")

    # ## 3. Building hyperparameter tuning grid
    tuning_grid = vec(collect(Base.product(
        tuning_space.maxnodes_range,
        tuning_space.nodesize_range,
        tuning_space.sampsize_range,
        tuning_space.mtry_range,
        tuning_space.ntrees_range
    )))

    dataset_partitions = Vector{Tuple{Vector{Int64}, Vector{Int64}}}()
    for fold_idx in 1:n_folds
        train, test = partitionvec(collect(1:nrow(prediction_df)),floor(Int64, nrow(prediction_df)/n_folds),fold_idx, true)
        push!(dataset_partitions, (train, test))
    end

    print("partitions calculated")

    # ys = Vector{Vector{Any}}()
    model_colnames = Symbol.("model_" .* map(x -> lpad(string(x), 6, "0"), collect(1:(n_folds*n_replicas*n_rngs))) )
    hyperpar_idxes = Vector{Int64}(undef, n_folds*n_replicas*n_rngs)
    replica_idxes = Vector{Int64}(undef, n_folds*n_replicas*n_rngs)
    fold_idxes = Vector{Int64}(undef, n_folds*n_replicas*n_rngs)
    rng_idxes = Vector{Int64}(undef, n_folds*n_replicas*n_rngs)
    train_rmses = Vector{Float64}(undef, n_folds*n_replicas*n_rngs)
    test_rmses = Vector{Float64}(undef, n_folds*n_replicas*n_rngs)
    train_mapes = Vector{Float64}(undef, n_folds*n_replicas*n_rngs)
    test_mapes = Vector{Float64}(undef, n_folds*n_replicas*n_rngs)
    train_cors = Vector{Float64}(undef, n_folds*n_replicas*n_rngs)
    test_cors = Vector{Float64}(undef, n_folds*n_replicas*n_rngs)
    importances = Matrix{Float64}(undef, length(input_cols), n_folds*n_replicas*n_rngs)
    shapley = Matrix{Float64}(undef, length(input_cols), n_folds*n_replicas*n_rngs)
    saved_predictions = Matrix{Float64}(undef, nrow(prediction_df), n_folds*n_replicas*n_rngs)
    saved_pertinences = Matrix{String}(undef, nrow(prediction_df), n_folds*n_replicas*n_rngs)
    # importances = DataFrame(:variable => [])
    # saved_predictions = DataFrame(unique_col => [])
    # saved_pertinences = DataFrame(unique_col => [])

    println("ProbeData portions initialized")

    model_index = 0

    for hp_idx in eachindex(tuning_grid)

        for rep_idx in 1:n_replicas

            Random.seed!(train_rng, rep_idx-1)
            shuffled_df = prediction_df[Random.randperm(train_rng, nrow(prediction_df)), :]
            X = shuffled_df[:, input_cols]
            y = shuffled_df[:, output_col]

            println("Dataset shuffled")

            for fold_idx in 1:n_folds

                train, test = dataset_partitions[fold_idx]

                for rng_idx in 1:n_rngs
    
                    Random.seed!(train_rng, rng_idx-1)

                    # rf_model = MLJScikitLearnInterface.RandomForestRegressor(
                    #     max_depth = nothing,
                    #     min_samples_leaf = tuning_grid[hp_idx][2],
                    #     max_features = 0.7,
                    #     n_estimators = tuning_grid[hp_idx][5]
                    # )

                    rf_model = MLJDecisionTreeInterface.RandomForestRegressor(
                        max_depth = tuning_grid[hp_idx][1],
                        min_samples_leaf = tuning_grid[hp_idx][2],
                        sampling_fraction = tuning_grid[hp_idx][3],
                        n_subfeatures = tuning_grid[hp_idx][4],
                        n_trees = tuning_grid[hp_idx][5],
                        rng = train_rng
                    )

                    # rf_model = NearestNeighborModels.KNNRegressor()

                    println("RF_model declared")
    
                    rf_machine = machine(rf_model, X[train, :], y[train], cache=false)
                    MLJ.fit!(rf_machine, verbosity=0)

                    println("Machine built and fit")

                    train_y_hat = MLJ.predict(rf_machine, X[train, :]) 
                    slopecorrection = lm(@formula(y ~ yhat), DataFrame(:y => y[train], :yhat => train_y_hat))
                    train_y_hat = GLM.predict(slopecorrection, DataFrame(:yhat => MLJ.predict(rf_machine, X[train, :]) ))

                    train_rmse = rmse(train_y_hat, y[train])
                    train_mape = mean(abs.(train_y_hat .- y[train]) ./ y[train])
                    train_cor = Statistics.cor(train_y_hat, y[train])

                    test_y_hat = GLM.predict(slopecorrection, DataFrame(:yhat => MLJ.predict(rf_machine, X[test, :])))
                    test_rmse = rmse(test_y_hat, y[test])
                    test_mape = mean(abs.(test_y_hat .- y[test]) ./ y[test])
                    test_cor = Statistics.cor(test_y_hat, y[test])

                    model_index += 1
                    model_colname = model_colnames[model_index]

                    importances[:, model_index] = impurity_importance(rf_machine.fitresult)
                    # importances[:, model_index] .= 0.0

                    data_shap = ShapML.shap(
                        explain = X[test, :],
                        model = rf_machine,
                        predict_function = ( (model, data) -> DataFrame(:y_pred => MLJ.predict(model, data)) ),
                        sample_size = 20,
                        seed = 1,
                        parallel = :features
                    )

                    #Cols are "index", "feature_name", "feature_value", "shap_effect", "shap_effect_sd", "intercept" from @show names(data_shap)
                    meanShapleyValues = @chain DataFrames.combine(
                        groupby(data_shap, :feature_name),
                        :shap_effect => (x -> mean(x)) => :shap_effect ;
                        renamecols = false
                    ) begin 
                        rename!( :feature_name => :variable, :shap_effect => model_colname)
                        sort!(:variable)
                    end

                    shapley[:, model_index] = meanShapleyValues[:, model_colname]

                    this_predictions = @chain DataFrame(
                        unique_col => vcat(shuffled_df[train, unique_col], shuffled_df[test, unique_col]),
                        model_colname => map(x -> round(x; digits = 1), vcat(train_y_hat, test_y_hat))
                    ) begin sort(unique_col) end
                    saved_predictions[:, model_index] = this_predictions[:, model_colname]

                    this_pertinences = @chain DataFrame(
                        unique_col => vcat(shuffled_df[train, unique_col], shuffled_df[test, unique_col]),
                        model_colname => vcat(repeat([ "train" ], length(train)),repeat([ "test" ], length(test)))
                    ) begin sort(unique_col) end
                    saved_pertinences[:, model_index] = this_pertinences[:, model_colname]

                    println("Merits calculated")

                    hyperpar_idxes[model_index] = hp_idx
                    replica_idxes[model_index] = rep_idx
                    fold_idxes[model_index] = fold_idx
                    rng_idxes[model_index] = rng_idx
                    train_rmses[model_index] = train_rmse
                    test_rmses[model_index] = test_rmse
                    train_mapes[model_index] = train_mape
                    test_mapes[model_index] = test_mape
                    train_cors[model_index] = train_cor
                    test_cors[model_index] = test_cor

                    println("Stuff saved, end of model $(model_index)")

                end
            end

            # push!(ys, y)

        end
    end

    report_df = DataFrame(
        :Hyperpar_Idx => hyperpar_idxes,
        :Replica_Idx => replica_idxes,
        :Fold_Idx => fold_idxes,
        :Rng_Idx => rng_idxes,
        :Train_RMSE=> train_rmses,
        :Test_RMSE => test_rmses,
        :Train_MAPE=> train_mapes,
        :Test_MAPE => test_mapes,
        :Train_Cor=> train_cors,
        :Test_Cor => test_cors
    )

    final_importances = @chain importances begin
        DataFrame(:auto)
        rename!( _, names( _ ) .=> model_colnames )
        hcat( DataFrame(:variable => names(prediction_df[:, input_cols])), _ )
    end

    final_shapley = @chain shapley begin
        DataFrame(:auto)
        rename!( _, names( _ ) .=> model_colnames )
        hcat( DataFrame(:variable => sort(names(prediction_df[:, input_cols]))), _ )
    end

    final_saved_predictions = @chain saved_predictions begin
        DataFrame(:auto)
        rename!( _, names( _ ) .=> model_colnames )
        hcat( DataFrame(unique_col => sort(prediction_df[:, unique_col])), _ )
    end

    final_saved_pertinences = @chain saved_pertinences begin
        DataFrame(:auto)
        rename!( _, names( _ ) .=> model_colnames )
        hcat( DataFrame(unique_col => sort(prediction_df[:, unique_col])), _ )
    end

    results = ProbeData(
        name,
        original_df,
        n_folds,
        n_replicas,
        n_rngs,
        report_df,
        final_importances,
        final_shapley,
        final_saved_predictions,
        final_saved_pertinences
    )

    return results
end

# ## 1. Only SES
# regression_currentCogScores_00to06mo_onlydemo = new_randomforest(
#     "regression_currentCogScores_00to06mo_onlydemo",
#     filtered_mdata_taxa_df_00to06,
#     identity,
#     [4, 5, 6],
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = onlydemo_tuning_space,
#     train_rng=ml_rng
# )
    
# jldsave(modelfiles("regression_currentCogScores_00to06mo_onlydemo.jld"); regression_currentCogScores_00to06mo_onlydemo)

# ## 2. Only taxonomic profiles
# regression_currentCogScores_00to06mo_onlytaxa = new_randomforest(
#     "regression_currentCogScores_00to06mo_onlytaxa",
#     filtered_mdata_taxa_df_00to06,
#     identity,
#     collect(6:ncol(filtered_mdata_taxa_df_00to06)), # 4 to include age
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = taxa_tuning_space,
#     train_rng=ml_rng
# )
    
# jldsave(modelfiles("regression_currentCogScores_00to06mo_onlytaxa.jld"); regression_currentCogScores_00to06mo_onlytaxa)

# ## 3. SES + taxonomic profiles
# regression_currentCogScores_00to06mo_demoplustaxa = new_randomforest(
#     "regression_currentCogScores_00to06mo_demoplustaxa",
#     filtered_mdata_taxa_df_00to06,
#     identity,
#     collect(4:ncol(filtered_mdata_taxa_df_00to06)),
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = taxa_tuning_space,
#     train_rng=ml_rng
# )

# jldsave(modelfiles("regression_currentCogScores_00to06mo_demoplustaxa.jld"); regression_currentCogScores_00to06mo_demoplustaxa)

# ## 4. Only functional profiles
# regression_currentCogScores_00to06mo_onlyecs = new_randomforest(
#     "regression_currentCogScores_00to06mo_onlyecs",
#     filtered_mdata_ecs_df_00to06,
#     identity,
#     collect(4:ncol(filtered_mdata_ecs_df_00to06)),
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = ecs_tuning_space,
#     train_rng=ml_rng
# )

# jldsave(modelfiles("regression_currentCogScores_00to06mo_onlyecs.jld"); regression_currentCogScores_00to06mo_onlyecs)

# ## 5. SES + functional profiles
# regression_currentCogScores_00to06mo_demoplusecs = new_randomforest(
#     "regression_currentCogScores_00to06mo_demoplusecs",
#     filtered_mdata_ecs_df_00to06,
#     identity,
#     collect(4:ncol(filtered_mdata_ecs_df_00to06)),
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = ecs_tuning_space,
#     train_rng=ml_rng
# )

# jldsave(modelfiles("regression_currentCogScores_00to06mo_demoplusecs.jld"); regression_currentCogScores_00to06mo_demoplusecs)

#####
# 18 to 120 months
#####

# ## 6. Only SES
regression_currentCogScores_18to120mo_onlydemo = new_randomforest(
    "regression_currentCogScores_18to120mo_onlydemo",
    filtered_mdata_taxa_df_18to120,
    identity,
    [4, 5, 6],
    :cogScore;
    n_folds = 3,
    n_replicas = 50,
    n_rngs = 5,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_18to120mo_onlydemo.jld"); regression_currentCogScores_18to120mo_onlydemo)

## 7. Only taxonomic profiles
regression_currentCogScores_18to120mo_onlytaxa = new_randomforest(
    "regression_currentCogScores_18to120mo_onlytaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(6:ncol(filtered_mdata_taxa_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 50,
    n_rngs = 5,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_18to120mo_onlytaxa.jld"); regression_currentCogScores_18to120mo_onlytaxa)

## 8. SES + taxonomic profiles
regression_currentCogScores_18to120mo_demoplustaxa = new_randomforest(
    "regression_currentCogScores_18to120mo_demoplustaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 50,
    n_rngs = 5,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_18to120mo_demoplustaxa.jld"); regression_currentCogScores_18to120mo_demoplustaxa)

# ## 9. Only functional profiles
# regression_currentCogScores_18to120mo_onlyecs = new_randomforest(
#     "regression_currentCogScores_18to120mo_onlyecs",
#     filtered_mdata_ecs_df_18to120,
#     identity,
#     collect(6:ncol(filtered_mdata_ecs_df_18to120)),
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = ecs_tuning_space,
#     train_rng=ml_rng
# )

# jldsave(modelfiles("regression_currentCogScores_18to120mo_onlyecs.jld"); regression_currentCogScores_18to120mo_onlyecs)

# ## 10. SES + functional profiles
# regression_currentCogScores_18to120mo_demoplusecs = new_randomforest(
#     "regression_currentCogScores_18to120mo_demoplusecs",
#     filtered_mdata_ecs_df_18to120,
#     identity,
#     collect(4:ncol(filtered_mdata_ecs_df_18to120)),
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = ecs_tuning_space,
#     train_rng=ml_rng
# )

# jldsave(modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld"); regression_currentCogScores_18to120mo_demoplusecs)

#####
# 18 to 120 months
#####

# ## 6. Only SES
regression_currentCogScores_00to120mo_onlydemo = new_randomforest(
    "regression_currentCogScores_00to120mo_onlydemo",
    filtered_mdata_taxa_df_00to120,
    identity,
    [4, 5, 6],
    :cogScore;
    n_folds = 3,
    n_replicas = 50,
    n_rngs = 5,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_00to120mo_onlydemo.jld"); regression_currentCogScores_00to120mo_onlydemo)

## 7. Only taxonomic profiles
regression_currentCogScores_00to120mo_onlytaxa = new_randomforest(
    "regression_currentCogScores_00to120mo_onlytaxa",
    filtered_mdata_taxa_df_00to120,
    identity,
    collect(6:ncol(filtered_mdata_taxa_df_00to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 50,
    n_rngs = 5,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
jldsave(modelfiles("regression_currentCogScores_00to120mo_onlytaxa.jld"); regression_currentCogScores_00to120mo_onlytaxa)

## 8. SES + taxonomic profiles
regression_currentCogScores_00to120mo_demoplustaxa = new_randomforest(
    "regression_currentCogScores_00to120mo_demoplustaxa",
    filtered_mdata_taxa_df_00to120,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_00to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 50,
    n_rngs = 5,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

jldsave(modelfiles("regression_currentCogScores_00to120mo_demoplustaxa.jld"); regression_currentCogScores_00to120mo_demoplustaxa)

# ## 9. Only functional profiles
# regression_currentCogScores_00to120mo_onlyecs = new_randomforest(
#     "regression_currentCogScores_00to120mo_onlyecs",
#     filtered_mdata_ecs_df_00to120,
#     identity,
#     collect(6:ncol(filtered_mdata_ecs_df_00to120)),
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = ecs_tuning_space,
#     train_rng=ml_rng
# )

# jldsave(modelfiles("regression_currentCogScores_00to120mo_onlyecs.jld"); regression_currentCogScores_00to120mo_onlyecs)

# ## 10. SES + functional profiles
# regression_currentCogScores_00to120mo_demoplusecs = new_randomforest(
#     "regression_currentCogScores_00to120mo_demoplusecs",
#     filtered_mdata_ecs_df_00to120,
#     identity,
#     collect(4:ncol(filtered_mdata_ecs_df_00to120)),
#     :cogScore;
#     n_folds = 3,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = ecs_tuning_space,
#     train_rng=ml_rng
# )

# jldsave(modelfiles("regression_currentCogScores_18to120mo_demoplusecs.jld"); regression_currentCogScores_18to120mo_demoplusecs)