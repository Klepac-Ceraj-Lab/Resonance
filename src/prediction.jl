########################################
# 0. Hierarchical Structures and Types #
########################################

abstract type Prediction end
struct Classification <: Prediction end
struct Regression <: Prediction end

struct ProbeData
    name::String
    original_df::DataFrame
    n_folds::Int64
    n_replicas::Int64
    n_rngs::Int64
    merits::DataFrame
    importances::DataFrame
    saved_predictions::DataFrame
    saved_pertinences::DataFrame
end

########################################
# 1. Preprocessing/Wrangling Functions #
########################################

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

## N-fold CV
function partitionvec(x,stride,i,longtail::Bool=true)
    # longtail=true to lengthen the last entry with the leftovers
    # longtail=false to place the leftovers in their own entry
    stride > 0 || error("stride must be positive") # doesn't handle negative strides
    starts = firstindex(x):stride:(lastindex(x)-longtail*stride)+1 # where to start each subvector
    subvecs = [view(x,starts[i]:get(starts,i+1,lastindex(x)+1)-1) for i in eachindex(starts)]
    train, test = reduce(vcat, [ collect(subvecs[map(x -> x != i, eachindex(subvecs))]) ]...), collect(subvecs[i])
    return train, test
end

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

function probe_prod_randomforest(
    ::Regression,
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
        nodesize_range = [ 5 ],
        sampsize_range = [ 0.7 ],
        mtry_range = [ -1 ],
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
        train, test = partitionvec(collect(1:nrow(prediction_df)),floor(Int64, nrow(prediction_df)/n_folds),fold_idx, true)
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
    saved_predictions = DataFrame(:subject => [])
    saved_pertinences = DataFrame(:subject => [])

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

                    this_predictions = DataFrame(
                        :subject => vcat(shuffled_df[train, :subject], shuffled_df[test, :subject]),
                        model_colname => map(x -> round(x; digits = 0), vcat(train_y_hat, test_y_hat))
                    )

                    this_pertinences = DataFrame(
                        :subject => vcat(shuffled_df[train, :subject], shuffled_df[test, :subject]),
                        model_colname => vcat(repeat([ "train" ], length(train)),repeat([ "test" ], length(test)))
                    )

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
                    saved_predictions = outerjoin(saved_predictions, this_predictions, on = [ :subject => :subject ])
                    saved_pertinences = outerjoin(saved_pertinences, this_pertinences, on = [ :subject => :subject ])
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
        original_df,
        n_folds,
        n_replicas,
        n_rngs,
        report_df,
        importances,
        sort(saved_predictions, :subject),
        sort(saved_pertinences, :subject)
    )

    return results
end

#########################################
# 3. Post-processing/Analysis Functions #
#########################################
struct CustomRangeNormalizer
    rmin::Float64
    rmax::Float64
    tmin::Float64
    tmax::Float64
end

## 3.0. Scale normalization
function compute_custom_scale(rvector::Vector{Float64}, tvector::Vector{Float64})
    return CustomRangeNormalizer(
        minimum(rvector),
        maximum(rvector),
        minimum(tvector),
        maximum(tvector)
    )
end

function normalize_number(m, scale::CustomRangeNormalizer)
    return ( (m - scale.rmin) / (scale.rmax - scale.rmin) * (scale.tmax - scale.tmin) + scale.tmin )
end

function scale_normalization(inputs::Vector{Float64}, scale::CustomRangeNormalizer)
    return map(x -> normalize_number(x, scale), inputs)
end

## 3.1. Generate Predictions (method predict)
