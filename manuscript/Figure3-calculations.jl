##########
# Concurrent cognitive assessment Score regression notebook
##########

#####
# 0. Loading libraries
#####
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
using Resonance
ml_rng = StableRNG(0)

#####
# Loading Data
#####

## 2. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles())
mdata_taxa_df = sort(Resonance.comm2wide(taxa), [ :subject, :timepoint ]);

open("taxa_present_echo_all.txt", "w") do io
    for i in [ f.name for f in features(taxa) ]
        println(io, i)
    end
end

open("taxa_present_echo_6mo.txt", "w") do io
    taxa_6mo = @chain mdata_taxa_df begin
        subset(:ageMonths => x -> x .< 6.0)
        select( _ , names( _ )[15:end])
        select( _ , names( _ )[[ sum(coll) > 0.0 for coll in eachcol(_) ]])
    end
    for i in names(taxa_6mo)
        println(io, i)
    end
end

open("taxa_present_echo_18mo.txt", "w") do io
    taxa_6mo = @chain mdata_taxa_df begin
        subset(:ageMonths => x -> 18.0 .< x .< 120.0)
        select( _ , names( _ )[15:end])
        select( _ , names( _ )[[ sum(coll) > 0.0 for coll in eachcol(_) ]])
    end
    for i in names(taxa_6mo)
        println(io, i)
    end
end

## 3. Functional Profiles
ecs = Resonance.load(ECProfiles())
mdata_ecs_df = sort(Resonance.comm2wide(ecs), [ :subject, :timepoint ]);

oldnames = names(mdata_ecs_df)[15:2414]

for oldname in names(mdata_ecs_df)[15:2414]
    rename!(mdata_ecs_df, oldname => randstring(['A':'Z'; 'a':'z'], 12))
end

newnames = names(mdata_ecs_df)[15:2414]

open("longnames.txt", "w") do io
    for i in oldnames
        println(io, i)
    end
end

open("hashnames.txt", "w") do io
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

function filter_prevalences(
    mdata_profile_df::DataFrame,

    label_col::Symbol,
    metadata_cols::Vector{Symbol},
    cols_to_delete::Vector{Symbol};
    lbound = 0.1,
    ubound = 1.0,
    education2int = true,
    sex2int = true
    )

    ## Step 1. Remove unnecessary stuff
    profile_df = select(mdata_profile_df, Not(vcat(label_col, metadata_cols, cols_to_delete)))
    
    ## Step 2. Calculate prevalences
    prevalences = map( x -> mean(x .> 0.0), eachcol(profile_df) )
    columns_to_keep = (prevalences .>= lbound) .& (prevalences .<= ubound)

    ## Step 3. update profile_df
    filtered_profile_df = profile_df[:, columns_to_keep]

    ## Step 4. append label and metadata to the beginning of the dataFrame
    filtered_profile_df = hcat(
        mdata_profile_df[:, vcat(label_col, metadata_cols)],
        filtered_profile_df
    )

    if education2int
        filtered_profile_df.education = coerce(int.(skipmissing(filtered_profile_df.education), type = Int), OrderedFactor)
    end

    if sex2int
        filtered_profile_df.sex = coerce(int.(skipmissing(filtered_profile_df.sex), type = Int), OrderedFactor)
    end

    return (filtered_profile_df)

end

## 00 to 06 months

filtered_mdata_taxa_df_00to06 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_00to06, :], Not([:race, :read_depth])),
    :cogScore,
    [:sex, :education, :ageMonths],
    [:sample, :sample_base, :subject, :timepoint, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_00to06 = filter_prevalences(
    dropmissing(mdata_ecs_df[mdata_ecs_df.filter_00to06, :], Not([:race, :read_depth])),
    :cogScore,
    [:sex, :education, :ageMonths],
    [:sample, :sample_base, :subject, :timepoint, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)

filtered_mdata_taxa_df_18to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_18to120, :], Not([:race, :read_depth])),
    :cogScore,
    [:sex, :education, :ageMonths],
    [:sample, :sample_base, :subject, :timepoint, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_mdata_ecs_df_18to120 = filter_prevalences(
    dropmissing(mdata_ecs_df[mdata_ecs_df.filter_18to120, :], Not([:race, :read_depth])),
    :cogScore,
    [:sex, :education, :ageMonths],
    [:sample, :sample_base, :subject, :timepoint, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = ecs_prevalence_threshold_lower,
    ubound = ecs_prevalence_threshold_upper
)

# ## Exports
CSV.write("00to06_taxa.csv", filtered_mdata_taxa_df_00to06)
CSV.write("00to06_ecs.csv", filtered_mdata_ecs_df_00to06)
CSV.write("18to120_taxa.csv", filtered_mdata_taxa_df_18to120)
CSV.write("18to120_ecs.csv", filtered_mdata_ecs_df_18to120)

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

# 1000 of each for onlydemo, onlytaxa, demoplustaxa

#####
# Functions
#####

struct ProbeData
    name::String
    n_folds::Int64
    n_replicas::Int64
    n_rngs::Int64
    merits::DataFrame
    importances::DataFrame
end

function probe_prod_randomforest(
    type::Resonance.Regression,
    name::String,
    original_df::DataFrame,
    data_preparation_fun::Function,
    input_cols,
    output_col;
    n_folds = 10,
    n_replicas = 10,
    n_rngs = 100,
    tuning_space = (;
        maxnodes_range = [ -1 ],
        nodesize_range = [ 0 ],
        sampsize_range = [ 0.7 ],
        mtry_range = [ 5 ],
        ntrees_range = [ 10 ]
    ),
    train_rng=Random.GLOBAL_RNG
    )

    # ## 1. Subsetting Data and shuffling DataFrame
    prediction_df = data_preparation_fun(original_df)

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
        train, test = Resonance.partitionvec(collect(1:nrow(prediction_df)),floor(Int64, nrow(prediction_df)/n_folds),fold_idx, true)
        push!(dataset_partitions, (train, test))
    end

    ys = Vector{Vector{Any}}()
    hyperpar_idxes = Vector{Int64}()
    replica_idxes = Vector{Int64}()
    fold_idxes = Vector{Int64}()
    rng_idxes = Vector{Int64}()
    train_rmses = Vector{Float64}()
    test_rmses = Vector{Float64}()
    train_mapes = Vector{Float64}()
    test_mapes = Vector{Float64}()
    train_cors = Vector{Float64}()
    test_cors = Vector{Float64}()
    importances = DataFrame(:variable => [])

    model_index = 0

    for hp_idx in eachindex(tuning_grid)

        for rep_idx in 1:n_replicas

            Random.seed!(train_rng, rep_idx-1)
            shuffled_df = prediction_df[Random.randperm(train_rng, nrow(prediction_df)), :]
            X = shuffled_df[:, input_cols]
            y = shuffled_df[:, output_col]

            for fold_idx in 1:n_folds

                train, test = dataset_partitions[fold_idx]

                for rng_idx in 1:n_rngs
    
                    Random.seed!(train_rng, rng_idx-1)
    
                    rf_model = RandomForestRegressor(
                        max_depth = tuning_grid[hp_idx][1],
                        min_samples_leaf = tuning_grid[hp_idx][2],
                        sampling_fraction = tuning_grid[hp_idx][3],
                        n_subfeatures = tuning_grid[hp_idx][4],
                        n_trees = tuning_grid[hp_idx][5],
                        rng = train_rng
                    )
    
                    rf_machine = machine(rf_model, X[train, :], y[train])
                    MLJ.fit!(rf_machine, verbosity=0)
            
                    train_y_hat = MLJ.predict(rf_machine, X[train, :]) 
                    train_rmse = rmse(train_y_hat, y[train])
                    train_mape = mean(abs.(train_y_hat .- y[train]) ./ y[train])
                    train_cor = Statistics.cor(train_y_hat, y[train])
                    test_y_hat = MLJ.predict(rf_machine, X[test, :]) 
                    test_rmse = rmse(test_y_hat, y[test])
                    test_mape = mean(abs.(test_y_hat .- y[test]) ./ y[test])
                    test_cor = Statistics.cor(test_y_hat, y[test])

                    model_index += 1
                    model_colname = Symbol("model_"*lpad(string(model_index),6,"0"))
    
                    push!(hyperpar_idxes, hp_idx)
                    push!(replica_idxes, rep_idx)
                    push!(fold_idxes, fold_idx)
                    push!(rng_idxes, rng_idx)
                    push!(train_rmses, train_rmse)
                    push!(test_rmses, test_rmse)
                    push!(train_mapes, train_mape)
                    push!(test_mapes, test_mape)
                    push!(train_cors, train_cor)
                    push!(test_cors, test_cor)
                    importances = outerjoin(importances, DataFrame(:variable => names(X), model_colname => impurity_importance(rf_machine.fitresult)), on = :variable)
                end
            end

            push!(ys, y)

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

    results = ProbeData(
        name,
        n_folds,
        n_replicas,
        n_rngs,
        report_df,
        importances
    )

    return results
end

#####
# 00 to 06 months
#####

## 1. Only SES

regression_currentCogScores_00to06mo_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_onlydemo",
    filtered_mdata_taxa_df_00to06,
    identity,
    [2, 3, 4],
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
JLD2.@save "models/aregression_currentCogScores_00to06mo_onlydemo.jld" regression_currentCogScores_00to06mo_onlydemo

## 2. Only taxonomic profiles

regression_currentCogScores_00to06mo_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_onlytaxa",
    filtered_mdata_taxa_df_00to06,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_00to06)), # 4 to include age
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
JLD2.@save "models/aregression_currentCogScores_00to06mo_onlytaxa.jld" regression_currentCogScores_00to06mo_onlytaxa

# ## 3. SES + taxonomic profiles

regression_currentCogScores_00to06mo_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_demoplustaxa",
    filtered_mdata_taxa_df_00to06,
    identity,
    collect(2:ncol(filtered_mdata_taxa_df_00to06)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

JLD2.@save "models/aregression_currentCogScores_00to06mo_demoplustaxa.jld" regression_currentCogScores_00to06mo_demoplustaxa

## 4. Only functional profiles

regression_currentCogScores_00to06mo_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_onlyecs",
    filtered_mdata_ecs_df_00to06,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_00to06)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

JLD2.@save "models/aregression_currentCogScores_00to06mo_onlyecs.jld" regression_currentCogScores_00to06mo_onlyecs

## 5. SES + functional profiles

regression_currentCogScores_00to06mo_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_00to06mo_demoplusecs",
    filtered_mdata_ecs_df_00to06,
    identity,
    collect(2:ncol(filtered_mdata_ecs_df_00to06)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

JLD2.@save "models/aregression_currentCogScores_00to06mo_demoplusecs.jld" regression_currentCogScores_00to06mo_demoplusecs

#####
# 18 to 120 months
#####

# ## 6. Only SES

regression_currentCogScores_18to120mo_onlydemo = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_onlydemo",
    filtered_mdata_taxa_df_18to120,
    identity,
    [2, 3, 4],
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = onlydemo_tuning_space,
    train_rng=ml_rng
)
    
JLD2.@save "models/aregression_currentCogScores_18to120mo_onlydemo.jld" regression_currentCogScores_18to120mo_onlydemo

# ## 7. Only taxonomic profiles

regression_currentCogScores_18to120mo_onlytaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_onlytaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(4:ncol(filtered_mdata_taxa_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)
    
JLD2.@save "models/aregression_currentCogScores_18to120mo_onlytaxa.jld" regression_currentCogScores_18to120mo_onlytaxa

# ## 8. SES + taxonomic profiles

regression_currentCogScores_18to120mo_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_demoplustaxa",
    filtered_mdata_taxa_df_18to120,
    identity,
    collect(2:ncol(filtered_mdata_taxa_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

JLD2.@save "models/aregression_currentCogScores_18to120mo_demoplustaxa.jld" regression_currentCogScores_18to120mo_demoplustaxa

## 9. Only functional profiles

regression_currentCogScores_18to120mo_onlyecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_onlyecs",
    filtered_mdata_ecs_df_18to120,
    identity,
    collect(4:ncol(filtered_mdata_ecs_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

JLD2.@save "models/aregression_currentCogScores_18to120mo_onlyecs.jld" regression_currentCogScores_18to120mo_onlyecs

## 10. SES + functional profiles

regression_currentCogScores_18to120mo_demoplusecs = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_currentCogScores_18to120mo_demoplusecs",
    filtered_mdata_ecs_df_18to120,
    identity,
    collect(2:ncol(filtered_mdata_ecs_df_18to120)),
    :cogScore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = ecs_tuning_space,
    train_rng=ml_rng
)

JLD2.@save "models/aregression_currentCogScores_18to120mo_demoplusecs.jld" regression_currentCogScores_18to120mo_demoplusecs