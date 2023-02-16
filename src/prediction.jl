########################################
# 0. Hierarchical Structures and Types #
########################################

struct CustomRangeNormalizer
    rmin::Float64
    rmax::Float64
    tmin::Float64
    tmax::Float64
end

abstract type Prediction end
struct Classification <: Prediction end
struct Regression <: Prediction end

abstract type Predictor end
abstract type UnivariatePredictor <: Resonance.Predictor end
abstract type MultivariatePredictor <: Resonance.Predictor end

mutable struct ExpandedRandomForestClassifier
    model::Machine
    train_accuracy::Float64
    test_accuracy::Float64
    fitness::Float64
end

mutable struct UnivariateRandomForestClassifier <: Resonance.UnivariatePredictor
    name::String
    original_data::DataFrame
    inputs_outputs::Tuple{DataFrame, CategoricalArrays.CategoricalVector{Bool, UInt32, Bool, CategoricalArrays.CategoricalValue{Bool, UInt32}, Union{}}}
    n_splits::Int64
    dataset_partitions::Vector{Tuple{Vector{Int64}, Vector{Int64}}}
    models:: Vector{Vector{ExpandedRandomForestClassifier}}
    selected_split::Tuple{Float64, Int64}
    train_accuracies::Vector{Float64}
    test_accuracies::Vector{Float64}
end

mutable struct UnivariateRandomForestRegressor <: Resonance.UnivariatePredictor
    name::String
    original_data::DataFrame
    inputs_outputs::Tuple{DataFrame, Vector{Float64}}
    n_splits::Int64
    dataset_partitions::Vector{Tuple{Vector{Int64}, Vector{Int64}}}
    models::Vector{Machine}
    scale_corrections::Vector{CustomRangeNormalizer}
    selected_split::Tuple{Float64, Int64}
    train_rmses::Vector{Float64}
    test_rmses::Vector{Float64}
    train_cors::Vector{Float64}
    test_cors::Vector{Float64}
end

mutable struct ExpandedRandomForestRegressor <: Resonance.UnivariatePredictor
    model::Machine
    scale_correction::CustomRangeNormalizer
    train_rmse::Float64
    test_rmse::Float64
    train_cor::Float64
    test_cor::Float64
    fitness::Float64
end

mutable struct UnivariatePredictorEnsemble
    screen_names::Vector{String}
    col_names::Vector{Symbol}
    predictors::Vector{T where T <: Resonance.UnivariatePredictor}
end

struct ProbeData
    name::String
    n_folds::Int64
    n_replicas::Int64
    n_rngs::Int64
    merits::DataFrame
    importances::DataFrame
end

########################################
# 1. Preprocessing/Wrangling Functions #
########################################

## 1.1. Helper/utility functions
dropnan(vv) = vv[.!(isnan.(vv))]

## 1.1.1. \xor used to aggregate columns that contain parts of the same original data
function myxor(a::Float64, b::Float64)
    a == b && (return a)
    a > b && (return a)
    a < b && (return b)
end

## 1.2. Raw data processing functions

### 1.2.1. Age bracket filtering for concurrent prediction
function filter_age_bracket(df, min_age, max_age)
    prediction_df = @chain df begin
        subset(:ageMonths => x -> min_age .<= x .< max_age)    
    end
    return prediction_df
end

### 1.2.2. Build a DataFrame for prediction of future target variables from mdata
function build_metadata_prediction_df(base_df::DataFrame, inputs::Vector{Symbol}, targets::Vector{Symbol}; keep_demographic = true)

    subjects = unique(base_df.subject)
    has_stool = findall(.!(ismissing.(base_df.sample)))
    has_cogscore = findall(.!(ismissing.(base_df.cogScorePercentile)))

    subjects_stool_idx = [ intersect(findall(base_df.subject .== el), has_stool) for el in subjects ]
    subjects_cogscore_idx = [ intersect(findall(base_df.subject .== el), has_cogscore) for el in subjects ]

    line_combinations = vcat( [ vec(collect(Base.product(subjects_stool_idx[i], subjects_cogscore_idx[i]))) for i in 1:length(subjects) ]... )
    line_combinations = line_combinations[ [ el[2] > el[1] for el in line_combinations ] ]

    targets_df = @chain base_df[ [el[2] for el in line_combinations], :] begin
        select( [:subject, :timepoint, :ageMonths, targets...] )
        rename!( [:timepoint => :futureTimepoint, :ageMonths => :futureAgeMonths, (targets .=> Symbol.("future" .* uppercasefirst.(String.(targets))))...] )
    end

    inputs_df = @chain base_df[ [el[1] for el in line_combinations], :] begin
        select( Not(:subject) )
    end

    prediction_df = hcat(targets_df, inputs_df)
    prediction_df.ageMonthsDelta = prediction_df.futureAgeMonths .- prediction_df.ageMonths 
    prediction_df.timepointDelta = prediction_df.futureTimepoint .- prediction_df.timepoint 
   
    if keep_demographic
        select!(
            prediction_df,
            [
                :subject, :timepoint, :ageMonths, 
                :sex, :education,
                :futureAgeMonths, :futureTimepoint, :ageMonthsDelta, :timepointDelta,
                targets..., Symbol.("future" .* uppercasefirst.(String.(targets)))..., inputs...
            ]
        )
    else
        select!(
            prediction_df,
            [
                :subject, :timepoint, :ageMonths,
                :futureAgeMonths, :futureTimepoint, :ageMonthsDelta, :timepointDelta,
                targets..., Symbol.("future" .* uppercasefirst.(String.(targets)))..., inputs...
            ]
        )    end
    sort!(prediction_df, [:subject, :timepoint, :futureTimepoint])
    
    return prediction_df

end

### 1.2.3. process dataframe built by function 1.2.2.
function prepare_future_prediction_df(df::DataFrame, inputs::Vector{Symbol}, targets::Vector{Symbol}, max_stool_ageMonths = 12.0, max_future_ageMonths = 24.0; filterout = [ :cogScorePercentile ] )

    future_prediction_df = @chain df begin
        build_metadata_prediction_df(inputs, targets)
        select( Not(filterout) ) # Toggle this and the next line for dropping/filtering original cogScore
    #    subset(:cogScore => x -> .!(ismissing.(x)))
        subset(:ageMonths => x -> .!(ismissing.(x)))
        subset(:futureAgeMonths => x -> x .<= max_future_ageMonths)
        subset(:ageMonths => x -> x .<= max_stool_ageMonths)
        unique([:subject, :timepoint])
        unique([:subject, :futureTimepoint])
        dropmissing()
    end

    return future_prediction_df

end

### 1.2.4. Meanclass (true/false if above/below vector mean)
meanclass(x::Vector{T} where T <: Real) = coerce(x .>= mean(x), OrderedFactor)

## 1.3. Univariate outlier detection

## 1.3.1. Computation of Tietjen-Moore r for outlier detection
function compute_tietjenmoore(data, k, n)
    r_all = abs.(data .- mean(data))
    filteredData = data[sortperm(r_all)][1:(n-k)]
    ksub = (filteredData .- mean(filteredData))

    return( sum(ksub .^ 2) / sum(r_all .^ 2) )
end

## 1.3.2. Tietjen-Moore outlier detection test
function test_tietjenmoore(dataSeries,k, n, alpha, sim_trials)
    ek = compute_tietjenmoore(dataSeries,k, n)
    simulation = [ compute_tietjenmoore(randn(length(dataSeries)), k, n) for i in 1:sim_trials ]
    Talpha=quantile(simulation,alpha)

    return(ek, Talpha)
end

## 1.3.3. Verbatim execution of Tietjen-Moore outlier test
function univariate_tietjenmoore(values::Vector{T} where T <: Real, k::Int64; alpha = 0.05, sim_trials = 100000)

    @info "----- Begin Tietjen-Moore Outlier test -----\n
        H0: There are no outliers in the data set.\n
        H1: There are exactly k outliers in the data set"

    n = length(values)
    L_set, L_critical = test_tietjenmoore(values, k, n, alpha, sim_trials)

    if L_set < L_critical
        @info "Set L-statistic for $n samples and $k outliers: $(round(L_set, digits = 4))\n
            Critical L for $n samples and $k outliers: $(round(L_critical, digits = 4))\n
            L_set < L_critical\n
            **SUCCESSFUL REJECTION OF H0** with confidence level $alpha" 
        r_all = abs.(values .- mean(values))
        outlier_indexes = sortperm(r_all)[(n-k+1):end]
        return outlier_indexes
    else
        @info """
            Set L-statistic for $n samples and $k outliers: $(round(L_set, digits = 4))
            Critical L for $n samples and $k outliers: $(round(L_critical, digits = 4))
            L_set > L_critical !
            **CANNOT REJECT H0** with confidence level $alpha
            """
        return Int64[]
    end # endif L_set < L_critical
end # end function

## 1.3.4. Automated Tietjen-Moore outlier test for differente numbers of outliers
function try_outliers(f, data, n; reverse=true)

    if reverse

        for i in n:-1:1
            outlier_idx = f(data, i)
            if length(outlier_idx) != 0
                return(i, outlier_idx)
            end
        end
    else

        for i in 1:n
            outlier_idx = f(data, i)
            if length(outlier_idx) != 0
                return(i, outlier_idx)
            end
         end
    end
    return 0, Int64[]
end

####################################
# 2. Training/Processing Functions #
####################################

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree

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

## 2.1. Train Tandom Forest Classifier

function train_randomforest(
    type::Resonance.Classification,
    ref_name::String,
    original_df::DataFrame,
    data_preparation_fun::Function,
    class_function::Function,
    input_cols,
    output_col;
    n_splits = 5,
    tuning_space = (;
        maxnodes_range = [ -1 ],
        nodesize_range = [ 0 ],
        sampsize_range = [ 0.7 ],
        mtry_range = [ 5 ],
        ntrees_range = [ 10 ]
    ),
    stress_draws=51,
    select_stress_percentile = 0.841,
    train_rng=Random.GLOBAL_RNG
    )

    # ## 1. Subsetting Data and shuffling DataFrame
    prediction_df = data_preparation_fun(original_df)
    Random.seed!(train_rng, 0)
    prediction_df = prediction_df[Random.randperm(train_rng, nrow(prediction_df)), :]
    
    # ## 2. Separating inputs/outputs
    X = prediction_df[:, input_cols]
    y = class_function(prediction_df[:, output_col])

    # ## 3. Building hyperparameter tuning grid
    tuning_grid = vec(collect(Base.product(
        tuning_space.maxnodes_range,
        tuning_space.nodesize_range,
        tuning_space.sampsize_range,
        tuning_space.mtry_range,
        tuning_space.ntrees_range
    )))

    # ## 4. Initializing meta-arrays to record tuning results for each split
    trial_partitions = Vector{Tuple{Vector{Int64}, Vector{Int64}}}(undef, n_splits)
    trial_machines = Vector{Vector{ExpandedRandomForestClassifier}}(undef, n_splits)
    trial_train_accuracies = zeros(Float64, n_splits)
    trial_test_accuracies = zeros(Float64, n_splits)

    # ## 5. Tuning/training loop
    @info "Performing $(n_splits) different train/test splits and tuning $(length(tuning_grid)) different hyperparmeter combinations\nfor the $(nrow(prediction_df)) samples."

    for this_trial in 1:n_splits

        fitness_threshold = 0.0

        ## 5.1. Splitting training data between train and test
        train, test = partitionvec(collect(1:nrow(prediction_df)),floor(Int64, nrow(prediction_df)/n_splits),this_trial, true)
        trial_partitions[this_trial] = (train, test)

        ## 5.2. Tuning hyperparameter for this split
        for i in eachindex(tuning_grid)

            expanded_models = try_stress_hyperparameters(
                Resonance.Classification(),
                X, y, train, test,
                (;
                    max_depth = tuning_grid[i][1],
                    min_samples_leaf = tuning_grid[i][2],
                    sampling_fraction = tuning_grid[i][3],
                    n_subfeatures = tuning_grid[i][4],
                    n_trees = tuning_grid[i][5]
                ),
                stress_draws;
                train_rng=train_rng)

            merits = report_merits(expanded_models)
            target_percentile = sort(merits[1:stress_draws, :Fitness])[floor(Int64, select_stress_percentile*stress_draws)]
            selected_fitness_index = findfirst(merits[1:stress_draws, :Fitness] .== target_percentile)
            selected_fitness = merits[selected_fitness_index, :Fitness]

            ## 5.2.4. Compare benchmark results with previous best model
            selected_fitness > fitness_threshold ? begin
                fitness_threshold = selected_fitness
                trial_train_accuracies[this_trial] = merits[selected_fitness_index, :Train_ACC]
                trial_test_accuracies[this_trial] = merits[selected_fitness_index, :Test_ACC]
                trial_machines[this_trial] = deepcopy(expanded_models)
            end : continue

        end # end for i in 1:length(tuning_grid) - each set of hyperparameters

        println("Split done!")

    end # end for this_trial - each different train/test split

    # ## 6. Returning optimization results
    results = UnivariateRandomForestClassifier(
        ref_name,                       #name::String
        prediction_df,                  #original_data::DataFrame
        (X,y),                          #inputs_outputs::Tuple{DataFrame, Vector{Float64}}
        n_splits,                       #n_splits::Int64
        trial_partitions,               #dataset_partitions::Vector{Tuple{Vector{Int64}, Vector{Int64}}}
        trial_machines,                 #models::Vector{Machine}
        findmax(trial_test_accuracies),  #selected_split::Tuple{Float64, Int64}
        trial_train_accuracies,         #train_accuracies::Vector{Float64}
        trial_test_accuracies,          #test_accuracies::Vector{Float64}
    )

    @info "Done!"
    return results

end # end function

## 2.2. Train Tandom Forest Regressor

function train_randomforest(
    type::Resonance.Regression,
    ref_name::String,
    original_df::DataFrame,
    data_preparation_fun::Function,
    input_cols,
    output_col;
    n_splits = 2,
    tuning_space = (;
        maxnodes_range = [ -1 ] ,
        nodesize_range = [ 0 ],
        sampsize_range = [ 0.7 ],
        mtry_range = [ 5 ],
        ntrees_range = [ 10 ]
    ),
    split_proportion=0.6,
    train_rng=Random.GLOBAL_RNG
    )

    # ## 1. Subsetting Data and shuffling DataFrame
    prediction_df = data_preparation_fun(original_df)
    Random.seed!(train_rng, 0)
    prediction_df = prediction_df[Random.randperm(train_rng, nrow(prediction_df)), :]
    
    # ## 2. Separating inputs/outputs
    X = prediction_df[:, input_cols]
    y = prediction_df[:, output_col]

    # ## 3. Building hyperparameter tuning grid
    tuning_grid = vec(collect(Base.product(
        tuning_space.maxnodes_range,
        tuning_space.nodesize_range,
        tuning_space.sampsize_range,
        tuning_space.mtry_range,
        tuning_space.ntrees_range
    )))

    # ## 4. Initializing meta-arrays to record tuning results for each split
    trial_partitions = Vector{Tuple{Vector{Int64}, Vector{Int64}}}(undef, n_splits)
    trial_machines = Vector{Machine}(undef, n_splits)
    trial_scalecorrections = Vector{CustomRangeNormalizer}(undef, n_splits)
    trial_train_rmses = repeat([Inf], n_splits)
    trial_test_rmses = repeat([Inf], n_splits)
    trial_train_cors = repeat([-1.0], n_splits)
    trial_test_cors = repeat([-1.0], n_splits)

    # ## 5. Tuning/training loop
    @info "Performing $(n_splits) different train/test splits and tuning $(length(tuning_grid)) different hyperparmeter combinations\nfor the $(nrow(prediction_df)) samples."

    for this_trial in 1:n_splits

        fitness_threshold = -1.0

        ## 5.1. Splitting training data between train and test
        # Random.seed!(train_rng, this_trial)
        # train, test = partition(eachindex(1:size(X, 1)), split_proportion, shuffle=true, rng=train_rng)
        train, test = partitionvec(collect(1:nrow(prediction_df)),floor(Int64, nrow(prediction_df)/n_splits),this_trial, true)
        trial_partitions[this_trial] = (train, test)

        ## 5.2. Tuning hyperparameter for this split
        for i in eachindex(tuning_grid)
            Random.seed!(train_rng, 0)

            ## 5.2.1. Construct model with a candidate set of hyperparameters
            rf_model = RandomForestRegressor(
                max_depth = tuning_grid[i][1],
                min_samples_leaf = tuning_grid[i][2],
                sampling_fraction = tuning_grid[i][3],
                n_subfeatures = tuning_grid[i][4],
                n_trees = tuning_grid[i][5],
                rng=train_rng
            )

            ## 5.2.2. Fit model on training data
            rf_machine = machine(rf_model, X[train, :], y[train])
            MLJ.fit!(rf_machine, verbosity=0)
            ## 5.2.2.1. Fit scale correction on training data
            train_y_hat = MLJ.predict(rf_machine, X[train, :])
            scale_correction = compute_custom_scale(train_y_hat, y[train])
            train_y_hat = scale_normalization(train_y_hat, scale_correction)
            train_rmse = rms(train_y_hat, y[train])
            train_cor = Statistics.cor(train_y_hat, y[train])

            ## 5.2.3. Benchmark model on independent testing data
            test_y_hat = MLJ.predict(rf_machine, X[test, :])
            test_y_hat = scale_normalization(test_y_hat, scale_correction)
            test_rmse = rms(test_y_hat, y[test])
            test_cor = Statistics.cor(test_y_hat, y[test])

            ## 5.2.4. Compare benchmark results with previous best model
            # model_fitness = sqrt(Complex(train_rmse*test_rmse))
            model_fitness = sqrt((train_cor^2)*(test_cor^2)) # Penalizing test_cor to avoid lucky test sets.
            # model_fitness = Complex(train_cor*test_cor*test_cor).re ^ 1.0/3.0

            ## 5.2.4. Compare benchmark results with previous best model
            # ((model_fitness.im == 0) & (model_fitness.re > fitness_threshold) & (train_cor != -1.0) & (test_cor != -1.0)) ? begin
            # ((model_fitness.re < fitness_threshold) & (train_cor > -0.0)) ? begin
            ((model_fitness > fitness_threshold) & (train_cor > -0.1) & (test_cor > -0.1)) ? begin
                println(model_fitness)
                println(train_cor)
                println(test_cor)
                fitness_threshold = model_fitness
                trial_train_rmses[this_trial] = train_rmse
                trial_train_cors[this_trial] = train_cor
                trial_test_rmses[this_trial] = test_rmse
                trial_test_cors[this_trial] = test_cor
                trial_machines[this_trial] = deepcopy(rf_machine)
                trial_scalecorrections[this_trial] = scale_correction
            end : continue

        end # end for i in 1:length(tuning_grid) - each set of hyperparameters

        println("Split done!")

    end # end for this_trial - each different train/test split

    # ## 6. Returning optimization results
    results = UnivariateRandomForestRegressor(
        ref_name,                   #name::String
        prediction_df,              #original_data::DataFrame
        (X,y),                      #inputs_outputs::Tuple{DataFrame, Vector{Float64}}
        n_splits,                   #n_splits::Int64
        trial_partitions,           #dataset_partitions::Vector{Tuple{Vector{Int64}, Vector{Int64}}}
        trial_machines,             #models::Vector{Machine}
        trial_scalecorrections,     #scale_correction::Vector{}
        findmin(trial_test_rmses),  #selected_split::Tuple{Float64, Int64}
        trial_train_rmses,          #train_rmses::Vector{Float64}
        trial_test_rmses,           #test_rmses::Vector{Float64}
        trial_train_cors,           #train_cor::Vector{Float64}
        trial_test_cors,            #test_cor::Vector{Float64}
    )

    @info "Done!"
    return results

end # end function

function try_stress_hyperparameters(
    ::Resonance.Classification,
    X,
    y,
    train,
    test,
    pars = (;
        maxnodes_range = -1 ,
        nodesize_range = 0,
        sampsize_range = 0.7,
        mtry_range = 5,
        ntrees_range = 1
    ),
    n_expanded_models = 101;
    train_rng=Random.GLOBAL_RNG)

    expanded_models = Vector{ExpandedRandomForestClassifier}(undef, n_expanded_models)

    for i in 1:n_expanded_models
        Random.seed!(train_rng, i-1)

        rf_model = RandomForestClassifier(
            max_depth = pars[1],
            min_samples_leaf = pars[2],
            sampling_fraction = pars[3],
            n_subfeatures = pars[4],
            n_trees = pars[5],
            rng=train_rng
        )
        rf_machine = machine(rf_model, X[train, :], y[train])
        MLJ.fit!(rf_machine, verbosity=0)

        train_y_hat = MLJ.predict_mode(rf_machine, X[train, :]) 
        train_acc = mean(train_y_hat .== y[train])
        test_y_hat = MLJ.predict_mode(rf_machine, X[test, :]) 
        test_acc = mean(test_y_hat .== y[test])
        fitness = sqrt(train_acc*test_acc)

        expanded_models[i] = ExpandedRandomForestClassifier(
            deepcopy(rf_machine),
            train_acc,
            test_acc,
            fitness
        )

    end

    return expanded_models

end

function expand_pretrained_model(
    model::Resonance.UnivariateRandomForestClassifier,
    selected_split::Int64,
    n_expanded_models::Int64;
    train_rng=Random.GLOBAL_RNG
    )

    X = model.inputs_outputs[1]
    y = model.inputs_outputs[2]
    train, test = model.dataset_partitions[selected_split]

    expanded_models = Vector{ExpandedRandomForestClassifier}(undef, n_expanded_models)

    for i in 1:n_expanded_models
        Random.seed!(train_rng, i-1)

        rf_model = RandomForestClassifier(
            max_depth = model.models[selected_split].model.max_depth,
            min_samples_leaf = model.models[selected_split].model.min_samples_leaf,
            sampling_fraction = model.models[selected_split].model.sampling_fraction,
            n_subfeatures = model.models[selected_split].model.n_subfeatures,
            n_trees = model.models[selected_split].model.n_trees,
            rng=train_rng
        )
        rf_machine = machine(rf_model, X[train, :], y[train])
        MLJ.fit!(rf_machine, verbosity=0)

        train_y_hat = MLJ.predict_mode(rf_machine, X[train, :]) 
        train_acc = mean(train_y_hat .== y[train])
        test_y_hat = MLJ.predict_mode(rf_machine, X[test, :]) 
        test_acc = mean(test_y_hat .== y[test])

        expanded_models[i] = ExpandedRandomForestClassifier(
            deepcopy(rf_machine),
            train_acc,
            test_acc
        )

    end

    @info "Done!"
    return expanded_models

end # end function

function expand_pretrained_model(
    model::Resonance.UnivariateRandomForestRegressor,
    selected_split::Int64,
    n_expanded_models::Int64;
    train_rng=Random.GLOBAL_RNG
    )

    X = model.inputs_outputs[1]
    y = model.inputs_outputs[2]
    train, test = model.dataset_partitions[selected_split]

    expanded_models = Vector{ExpandedRandomForestRegressor}(undef, n_expanded_models)

    for i in 1:n_expanded_models
        Random.seed!(train_rng, i-1)

        rf_model = RandomForestRegressor(
            max_depth = model.models[selected_split].model.max_depth,
            min_samples_leaf = model.models[selected_split].model.min_samples_leaf,
            sampling_fraction = model.models[selected_split].model.sampling_fraction,
            n_subfeatures = model.models[selected_split].model.n_subfeatures,
            n_trees = model.models[selected_split].model.n_trees,
            rng=train_rng
        )
        rf_machine = machine(rf_model, X[train, :], y[train])
        MLJ.fit!(rf_machine, verbosity=0)

        train_y_hat = MLJ.predict(rf_machine, X[train, :])
        scale_correction = compute_custom_scale(train_y_hat, y[train])
        train_y_hat = scale_normalization(train_y_hat, scale_correction)
        train_rmse = rms(train_y_hat, y[train])
        train_cor = Statistics.cor(train_y_hat, y[train])

        test_y_hat = MLJ.predict(rf_machine, X[test, :])
        test_y_hat = scale_normalization(test_y_hat, scale_correction)
        test_rmse = rms(test_y_hat, y[test])
        test_cor = Statistics.cor(test_y_hat, y[test])

        expanded_models[i] = ExpandedRandomForestRegressor(
            deepcopy(rf_machine),
            scale_correction,
            train_rmse,
            test_rmse,
            train_cor,
            test_cor
        )

    end

    @info "Done!"
    return expanded_models

end # end function

#########################################
# 3. Post-processing/Analysis Functions #
#########################################

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

## 3.1.1. Generate Predictions with the UnivariateRandomForestClassifier
function predict_proba(model::UnivariateRandomForestClassifier, newx=nothing; split_index=0)

    (split_index > length(model.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = model.selected_split[2])

    if newx isa Nothing
        return MLJ.predict(model.models[split_index], model.inputs_outputs[1])
    else
        return MLJ.predict(model.models[split_index], newx)
    end

end

function predict(model::UnivariateRandomForestClassifier, newx=nothing; split_index=0)

    (split_index > length(model.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = model.selected_split[2])

    if newx isa Nothing
        return MLJ.predict_mode(model.models[split_index], model.inputs_outputs[1])
    else
        return MLJ.predict_mode(model.models[split_index], newx)
    end

end

## 3.1.2. Generate Predictions with the UnivariateRandomForestRegressor
function predict(model::UnivariateRandomForestRegressor, newx=nothing; split_index=0)

    (split_index > length(model.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = model.selected_split[2])

    scale_correction = model.scale_corrections[split_index]

    if newx isa Nothing
        return scale_normalization(MLJ.predict(model.models[split_index], model.inputs_outputs[1]), scale_correction)
    else
        return scale_normalization(MLJ.predict(model.models[split_index], newx), scale_correction)
    end

    # ## Without Scale Correction
    # if newx isa Nothing
    #     return MLJ.predict(model.models[split_index], model.inputs_outputs[1])
    # else
    #     return MLJ.predict(model.models[split_index], newx)
    # end

end

## 3.2. Report Figures of Merit

### 3.2.1. Report Figures of Merit for Classification

function report_merits(res::UnivariateRandomForestClassifier)

    merits = DataFrame(
    :Split => collect(1:res.n_splits),
    :Train_ACC => res.train_accuracies,
    :Test_ACC => res.test_accuracies,
    )

    push!(merits, Dict(
        :Split => 0,
        :Train_ACC => mean(merits.Train_ACC),
        :Test_ACC => mean(merits.Test_ACC)
    ))

    push!(merits, Dict(
        :Split => -1,
        :Train_ACC => median(merits.Train_ACC[1:end-1]),
        :Test_ACC => median(merits.Test_ACC[1:end-1])
    ))

    return merits
end

### 3.2.1. Report Figures of Merit for Regression
function report_merits(res::UnivariateRandomForestRegressor)
    
    merits = DataFrame(
        :Split => collect(1:res.n_splits),
        :Train_RMSE => res.train_rmses,
        :Test_RMSE => res.test_rmses,
        :Train_COR => res.train_cors,
        :Test_COR => res.test_cors
        )

    push!(merits, Dict(
        :Split => 0,
        :Train_RMSE => mean(merits.Train_RMSE),
        :Test_RMSE => mean(merits.Test_RMSE),
        :Train_COR => mean(merits.Train_COR),
        :Test_COR => mean(merits.Test_COR)
    ))

    push!(merits, Dict(
        :Split => -1,
        :Train_RMSE => median(merits.Train_RMSE[1:end-1]),
        :Test_RMSE => median(merits.Test_RMSE[1:end-1]),
        :Train_COR => median(merits.Train_COR[1:end-1]),
        :Test_COR => median(merits.Test_COR[1:end-1])
    ))
    
    return merits
end

function report_merits(res::Vector{ExpandedRandomForestClassifier})

    merits = DataFrame(
    :Trial => collect(1:length(res)),
    :Train_ACC => [ el.train_accuracy for el in res ],
    :Test_ACC => [ el.test_accuracy for el in res ],
    :Fitness => [ el.fitness for el in res ],
    )

    push!(merits, Dict(
        :Trial => 0,
        :Train_ACC => mean(merits.Train_ACC),
        :Test_ACC => mean(merits.Test_ACC),
        :Fitness => mean(merits.Fitness)
    ))

    push!(merits, Dict(
        :Trial => -1,
        :Train_ACC => median(merits.Train_ACC[1:end-1]),
        :Test_ACC => median(merits.Test_ACC[1:end-1]),
        :Fitness => median(merits.Fitness[1:end-1])
    ))

    return merits
end

function report_merits(res::Vector{ExpandedRandomForestRegressor})

    merits = DataFrame(
    :Trial => collect(1:length(res)),
    :Train_RMSE => [ el.train_rmse for el in res ],
    :Test_RMSE => [ el.test_rmse for el in res ],
    :Train_COR => [ el.train_cor for el in res ],
    :Test_COR => [ el.test_cor for el in res ],
    )

    push!(merits, Dict(
        :Trial => 0,
        :Train_RMSE => mean(merits.Train_RMSE),
        :Test_RMSE => mean(merits.Test_RMSE),
        :Train_COR => mean(merits.Train_COR),
        :Test_COR => mean(merits.Test_COR)
    ))

    push!(merits, Dict(
        :Trial => -1,
        :Train_RMSE => median(merits.Train_RMSE[1:end-1]),
        :Test_RMSE => median(merits.Test_RMSE[1:end-1]),
        :Train_COR => median(merits.Train_COR[1:end-1]),
        :Test_COR => median(merits.Test_COR[1:end-1])
    ))

    return merits
end
## 3.3. Importance analysis

### 3.3.1. MDI importance for all features in a single train/test split of a single Classifier model
function singlemodel_singlesplit_importance(res::UnivariateRandomForestClassifier; split_index=0)

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    importance_df = DataFrame(
        :Variable => names(res.inputs_outputs[1]),
        :Importance => DecisionTree.impurity_importance(res.models[split_index].fitresult[1])
    )

    sort!(importance_df, :Importance, rev = true)

    return importance_df

end

### 3.3.2. MDI importance for all features in a single train/test split of a single Regressor model
function singlemodel_singlesplit_importance(res::UnivariateRandomForestRegressor; split_index = 0 )

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    importance_df = DataFrame(
        :Variable => names(res.inputs_outputs[1]),
        :Importance => DecisionTree.impurity_importance(res.models[split_index].fitresult)
    )

    sort!(importance_df, :Importance, rev = true)

    return importance_df
    
end

### 3.3.3. Concatenated (hcat) MDI importances for all features in all train/test splits of a single UnivariatePredictor model
function singlemodel_allsplits_importances(res::T where T <: Resonance.UnivariatePredictor)

    singlesplit_importances = [ singlemodel_singlesplit_importance(res; split_index = i) for i in 1:length(res.models) ]
    concatenated_importances_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlesplit_importances)

    return concatenated_importances_df

end

### 3.3.4. Averaged MDI importances for all features in all train/test splits of a single UnivariatePredictor model
function singlemodel_summary_importances(res::T where T <: Resonance.UnivariatePredictor, colname = :AvgImportance, fun = vv -> mean(dropnan(collect(skipmissing(vv)))))

    concatenated_importances_df = singlemodel_allsplits_importances(res)
    summarised_importances_df = DataFrame(
       :Variable => concatenated_importances_df.Variable,
       colname => map(fun, eachrow(Matrix(concatenated_importances_df[:, 2:end])))
    )

    sort!(summarised_importances_df, colname, rev = true)

    return summarised_importances_df

end

### 3.3.5. Concatenated (hcat) averaged MDI importances for all features in all train/test splits of multiple UnivariatePredictor models in an ensemble
function multimodel_individual_summaryimportances(ens::UnivariatePredictorEnsemble, col_prefix = "meanImportance_", fun = vv -> mean(dropnan(collect(skipmissing(vv)))))

    singlemodel_summaries = [ 
        singlemodel_summary_importances(
            ens.predictors[i],
            Symbol(col_prefix * string(ens.col_names[i])),
            fun
        ) for i in 1:length(ens.predictors)
    ]

    concatenated_summaries_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlemodel_summaries)

    return concatenated_summaries_df

end

### 3.3.6. Averaged MDI importances for all features in all train/test splits of multiple UnivariatePredictor models in an ensemble
function multimodel_aggregate_summaryimportances(
    ens::UnivariatePredictorEnsemble,
    singlemodel_col_prefix = "meanImportance_",
    aggregate_colname = :AvgMultimodelImportance,
    singlemodel_summary_fun = vv -> mean(dropnan(collect(skipmissing(vv)))),
    multimodel_summary_fun = vv -> mean(dropnan(collect(skipmissing(vv))))
    )

    concatenated_summaries_df = multimodel_individual_summaryimportances(ens, singlemodel_col_prefix, singlemodel_summary_fun)

    summarised_importances_df = DataFrame(
        :Variable => concatenated_summaries_df.Variable,
        aggregate_colname => map(multimodel_summary_fun, eachrow(Matrix(concatenated_summaries_df[:, 2:end])))
    )

    sort!(summarised_importances_df, aggregate_colname, rev = true)

    return summarised_importances_df

end

### 3.3.7. Binary column (if feature in the top N average importances) for all features in all train/test splits of a single Regressor model
function singlemodel_binarytopn_importances(
    res::T where T <: Resonance.UnivariatePredictor,
    importance_colname = :AvgImportance,
    topn_colname = :TopN,
    fun = vv -> mean(dropnan(collect(skipmissing(vv))));
    n=50
    )

    summarised_importances_df = @chain singlemodel_summary_importances(res, importance_colname, fun) begin
        insertcols!( _ , 2, topn_colname => vcat(ones(Int64, n), zeros(Int64, nrow( _ )-n)) )
        select!( Not(importance_colname) )
    end
    
    return(summarised_importances_df)

end

### 3.3.8. Concatenated (hcat) binary columns (if feature in the top N average importances) for all features in all train/test splits of multiple UnivariatePredictor models in an ensemble
function multimodel_individual_binarytopns(
    ens::UnivariatePredictorEnsemble,
    col_prefix = "topN_",
    fun = vv -> mean(dropnan(collect(skipmissing(vv))));
    n=50
    )

    singlemodel_summaries = [ 
        singlemodel_binarytopn_importances(
            ens.predictors[i],
            :AvgImportance,
            Symbol(col_prefix * string(ens.col_names[i])),
            fun;
            n=n
        ) for i in 1:length(ens.predictors)
    ]

    concatenated_summaries_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlemodel_summaries)

    return concatenated_summaries_df

end

### 3.3.9. *SUM* of binary columns (if feature in the top N average importances) for all features in all train/test splits of multiple UnivariatePredictor models in an ensemble
function multimodel_aggregate_binarytopns(
    ens::UnivariatePredictorEnsemble,
    singlemodel_col_prefix = "topN_",
    aggregate_colname = :SumTopN,
    singlemodel_summary_fun = vv -> mean(dropnan(collect(skipmissing(vv)))),
    multimodel_summary_fun = sum;
    n=30)

    concatenated_summaries_df = multimodel_individual_binarytopns(ens, singlemodel_col_prefix, singlemodel_summary_fun; n=n)

    summarised_importances_df = DataFrame(
        :Variable => concatenated_summaries_df.Variable,
        aggregate_colname => map(multimodel_summary_fun, eachrow(Matrix(concatenated_summaries_df[:, 2:end])))
    )

    sort!(summarised_importances_df, aggregate_colname, rev = true)

    return summarised_importances_df

end

## 3.4. Input set split-independent statistics

### 3.4.1. General input descriptor
function descript_inputs(res::T where T <: Resonance.Predictor)

    X = res.inputs_outputs[1]

    bug_names = names(X)
    bug_prevalences = map( x -> sum( x .> 0.0) / length(x), eachcol(X) )
    bug_averages_withzero = map(x -> mean(x), eachcol(X) )
    bug_averages_nonzero = map(x -> mean(x[x .> 0.0]), eachcol(X) )
    bug_std_withzeros = map( x -> Statistics.std(x), eachcol(X) )
    bug_std_nonzeros = map( x -> Statistics.std(x[x .> 0.0]), eachcol(X) )

    descript_df = DataFrame(
        :Variable => bug_names,
        :Prevalence => bug_prevalences,
        :MeanAbundanceConsideringZeros => bug_averages_withzero,
        :MeanAbundanceExcludingZeros => bug_averages_nonzero,
        :SdevConsideringZeros => bug_std_withzeros,
        :SdevExcludingZeros => bug_std_nonzeros
    )

    return descript_df

end

### 3.4.2. Collect average bug prevalences on the entire input feature set

function singlemodel_summary_prevalences(res::T where T <: Resonance.UnivariatePredictor, colname = :Prevalence)

    descript_df = descript_inputs(res)
    prevalence_df = DataFrame(
        :Variable => descript_df.Variable,
        colname => descript_df.Prevalence
    )

    return prevalence_df

end

### 3.4.3. Collect average bug abundances on the entire input feature set

function singlemodel_summary_abundances(res::T where T <: Resonance.UnivariatePredictor, colname = :MeanAbundance; exclude_zeros=true)

    descript_df = descript_inputs(res)

    if exclude_zeros
        abundances_df = DataFrame(
            :Variable => descript_df.Variable,
            colname => descript_df.MeanAbundanceExcludingZeros
        )
    else
        abundances_df = DataFrame(
            :Variable => descript_df.Variable,
            colname => descript_df.MeanAbundanceConsideringZeros
        )
    end    

    abundances_df[isnan.(abundances_df[:, colname]), colname] .= 0.0

    return abundances_df

end

### 3.4.4. Collect standard deviations of bug abundances on the entire input feature set

function singlemodel_summary_sdevs(res::T where T <: Resonance.UnivariatePredictor, colname = :Sdev; exclude_zeros=true)

    descript_df = descript_inputs(res)

    if exclude_zeros
        sdevs_df = DataFrame(
            :Variable => descript_df.Variable,
            colname => descript_df.SdevExcludingZeros
        )
    else
        sdevs_df = DataFrame(
            :Variable => descript_df.Variable,
            colname => descript_df.SdevConsideringZeros
        )
    end    

    sdevs_df[isnan.(sdevs_df[:, colname]), colname] .= 0.0

    return sdevs_df
    
end

### 3.4.5. Collect summary of bug prevalences on the entire input feature set, for each model in an ensemble of models
function multimodel_individual_prevalences(ens::UnivariatePredictorEnsemble, col_prefix = "Prevalence_")

    singlemodel_prevalences = [ 
        singlemodel_summary_prevalences(
            ens.predictors[i],
            Symbol(col_prefix * string(ens.col_names[i]))
        ) for i in 1:length(ens.predictors)
    ]

    concatenated_prevalences_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlemodel_prevalences)

    return concatenated_prevalences_df

end

### 3.4.6. Collect summary of bug abundances on the entire input feature set, for each model in an ensemble of models
function multimodel_individual_abundances(ens::UnivariatePredictorEnsemble, col_prefix = "Abundance_"; exclude_zeros=true)

    singlemodel_abundances = [ 
        singlemodel_summary_abundances(
            ens.predictors[i],
            Symbol(col_prefix * string(ens.col_names[i]))
            ;exclude_zeros=exclude_zeros
        ) for i in 1:length(ens.predictors)
    ]

    concatenated_abundances_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlemodel_abundances)

    return concatenated_abundances_df

end

## 3.5. Correlation Analysis

### 3.5.1. Biserial correlation function
function biserial_correlation( ## Borrowed from github.com/m3g/XCorrelation.jl =)
	discreteArray,
	continuousArray :: Array{Float64})

	## Sanity check
	
	length(discreteArray) == length(continuousArray) || throw(DimensionMismatch("Size of discrete Array not equal to size of continuous Array"))

	## Computation

	Ntrue = sum(discreteArray .== 1)
	Nfalse = sum(discreteArray .== 0)

	SmeanTrue = Statistics.mean(continuousArray[discreteArray .== 1])
	SmeanFalse = Statistics.mean(continuousArray[discreteArray .== 0])

	sdev = Statistics.std(continuousArray, corrected=false)

	bisCorr = ((SmeanTrue - SmeanFalse)/sdev)*sqrt(Ntrue*Nfalse/(length(discreteArray)^2))

	(bisCorr == NaN) && (bisCorr = -0.)

	return( bisCorr )

end

### 3.5.2. Point-biserial correlations between a nodesplit result (binary larger than or smaller than threshold), and the target output value,
#   for a specific forest (set of trees), considering only samples that reached a node containing that specific feature for decision.
function feature_split_correlation_analysis(X, y, trees)

    feature_correlations_matrix = Matrix{Union{Float64, Missing}}(undef, length(trees), ncol(X))

    for tree_idx in eachindex(trees) # For each tree

        sample_feature_rel = zeros(Bool, nrow(X), ncol(X)) # Is feature j used to classify sample i in this tree?
        feature_split_rel = Matrix{Union{Int64, Missing}}(undef, nrow(X), ncol(X))  # Is the value of feature j for sample i higher or lower than the decision threshold ?
        decision_split_feature = zeros(Int64, nrow(X))

        for sample_idx in 1:nrow(X) # For each sample, make a tree pass.

            this_tree = deepcopy(trees[tree_idx])

            while( !( typeof(this_tree) <: DecisionTree.Leaf ) )
            # while( !( ( typeof(this_tree.left) <: DecisionTree.Leaf ) & ( typeof(this_tree.right) <: DecisionTree.Leaf ) ) )

                if this_tree.featid == 0
                    this_tree = this_tree.left                
                elseif X[sample_idx, this_tree.featid] < this_tree.featval
                    sample_feature_rel[sample_idx, this_tree.featid] = true
                    feature_split_rel[sample_idx, this_tree.featid] = 0
                    decision_split_feature[sample_idx] = this_tree.featid
                    this_tree = this_tree.left                
                else
                    sample_feature_rel[sample_idx, this_tree.featid] = true
                    feature_split_rel[sample_idx, this_tree.featid] = 1
                    this_tree = this_tree.right
                end # end state update
        
            end # end while still a node

        end # end for sample_idx

        for feat_idx in 1:ncol(X)

            samples_to_consider = sample_feature_rel[:, feat_idx]

            if sum(samples_to_consider) != 0 # if there are no samples that considered this feature
                feature_values = X[samples_to_consider, feat_idx]
                split_values = convert(Vector{Int64}, feature_split_rel[samples_to_consider, feat_idx])
                target_values = convert(Vector{Float64}, y[samples_to_consider])

                if ( (Statistics.std(feature_values) != 0.0) & (Statistics.std(target_values) != 0.0) )
                    feature_correlations_matrix[tree_idx, feat_idx] = biserial_correlation(split_values, target_values)
                end
                
            else
                continue
            end # end if sum(samples_to_consider) != 0

        end # end for this_feature   

    end # end for this_tree

    # @show feature_correlations_matrix
    return feature_correlations_matrix

end

### 3.5.3. Function that averages every column of the result of the 3.5.2. function (where every tree is a row),
#   to produce the average correlation for a certain forest (created via a certain split), for a Classification model.
function singlemodel_singlesplit_correlations(
    res::UnivariateRandomForestClassifier,
    ; split_index=0)

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    model_X, model_y = res.inputs_outputs

    model_trees = res.models[split_index].fitresult[1].trees
    feature_correlations_matrix = feature_split_correlation_analysis(model_X, model_y, model_trees)
    feature_mean_correlations = map( x -> mean(dropnan(collect(skipmissing(x)))), eachcol(feature_correlations_matrix))

    correlations_df = DataFrame(
        :Variable => names(model_X),
        :AvgCorrelation => feature_mean_correlations,
    )

    sort!(correlations_df, :AvgCorrelation; rev=true)

    return correlations_df

end

### 3.5.4. Function that averages every column of the result of the 3.5.2. function (where every tree is a row),
#   to produce the average correlation for a certain forest (created via a certain split), for a Regression model.
function singlemodel_singlesplit_correlations(
    res::UnivariateRandomForestRegressor,
    ; split_index=0)

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    model_X, model_y = res.inputs_outputs

    model_trees = res.models[split_index].fitresult.trees
    feature_correlations_matrix = feature_split_correlation_analysis(model_X, model_y, model_trees)
    feature_mean_correlations = map( x -> mean(dropnan(collect(skipmissing(x)))), eachcol(feature_correlations_matrix))

    correlations_df = DataFrame(
        :Variable => names(model_X),
        :AvgCorrelation => feature_mean_correlations,
    )

    sort!(correlations_df, :AvgCorrelation; rev=true)

    return correlations_df

end

### 3.5.5. Concatenated result of the split-correlation anaysis, along all the forests created by all the train/test splits.
function singlemodel_allsplits_correlations(res::T where T <: Resonance.UnivariatePredictor)

    singlesplit_correlations = [ singlemodel_singlesplit_correlations(res; split_index = i) for i in 1:length(res.models) ]
    concatenated_correlations_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlesplit_correlations)

    return concatenated_correlations_df

end

### 3.5.6. Average result of the split-correlation anaysis, along all the forests created by all the train/test splits.
function singlemodel_allsplits_correlationsummary(res::T where T <: Resonance.UnivariatePredictor, colname = :AvgCorrelation, fun = nonna_nonmissing_mean)

    concatenated_correlations_df = singlemodel_allsplits_correlations(res)
    summarised_correlations_df = DataFrame(
       :Variable => concatenated_correlations_df.Variable,
       colname => map(fun, eachrow(Matrix(concatenated_correlations_df[:, 2:end])))
    )

    sort!(summarised_correlations_df, colname, rev = true)

    return summarised_correlations_df

end

### 3.5.7. Concatenated result of the split-correlation anaysis, along all the forests created by all the train/test splits of all models on an ensemble.
function multimodel_individual_correlations(ens::UnivariatePredictorEnsemble, col_prefix = "Correlation_")

    singlemodel_correlations = [ 
        singlemodel_allsplits_correlationsummary(
            ens.predictors[i],
            Symbol(col_prefix * string(ens.col_names[i]))
        ) for i in 1:length(ens.predictors)
    ]

    concatenated_correlations_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlemodel_correlations)

    return concatenated_correlations_df

end

#########################
# 4. Plotting functions #
#########################

## 4.1. Merit (benchmark) plots

### 4.1.1. Barplot of the confusion matrix, averaged over all train/test splits, for a single model.

#### 4.1.1.1. Function to build the Confusion Matrix for a single train/test split on a single model.
function build_confusion_matrix(res::UnivariateRandomForestClassifier, selected_split::Int)

    y = res.inputs_outputs[2][res.dataset_partitions[selected_split][2]]
    yhat = MLJ.predict_mode(res.models[selected_split], res.inputs_outputs[1][res.dataset_partitions[selected_split][2],:])
    confmat = ConfusionMatrix()(yhat, y)

    return confmat
end

#### 4.1.1.2. Function to average the Confusion Matrices for all train/test splits on a single model.
function average_confusion_matrix(res::UnivariateRandomForestClassifier)

    matrix_vector = [ build_confusion_matrix(res, i).mat for i in 1:res.n_splits ]
    concatenated_matrix = vcat([ vec(mat)' for mat in matrix_vector ]...) # Column order is TN, FP, FN, TP !
    averages = [ sum(concatenated_matrix[:,i])/sum(concatenated_matrix) for i in 1:4] # 4 columns: TN, FP, FN, TP !

    return averages
end

#### 4.1.1.3. Function to build the vectors used for plotting the avereages of confusion matrices
function confmatrix2barplot(res::UnivariateRandomForestClassifier)

    confmat_inputs = (
        x = [1, 2, 2, 1], # Column order is TN, FP, FN, TP !
        value = average_confusion_matrix(res),
        grp = [1, 2, 2, 1],
        color = [1, 2, 3, 4]
    )

    return confmat_inputs
end

### 4.1.1.4. The actual plotting function
function singlemodel_merit_barplot!(
    figure::Figure,
    res::UnivariateRandomForestClassifier,
    pos::Tuple{Int64, Int64},
    plot_title::String;
    confplot_colors = [:blue3, :red3, :tomato2, :dodgerblue2]
    )

    ax = Axis(
        figure[pos[1],pos[2]];
        xticks = (1:2, ["Correct", "Incorrect"]),
        ylabel = "Proportion",
        title = plot_title
    )
    ylims!(ax, [0.0, 1.0])
    
    tbl = confmatrix2barplot(res)
    barplot!(
        ax, tbl.x, tbl.value,
        stack = tbl.grp,
        color = confplot_colors[tbl.color],
        size = 5
    )
    return figure
end # end function

## 4.1.2. Scatterplot of ground truth vs prediction, for the most performant train/test split, for a single model.
function singlemodel_merit_scatterplot!(
    ax::Axis,
    res::UnivariateRandomForestRegressor;
    split_index = 0,
    traincorr_pos = Point(0.0, 1.0),
    testcorr_pos = Point(0.0, 0.9))

    # Plot barplot

    y = res.inputs_outputs[2]
    yhat = Resonance.predict(res, res.inputs_outputs[1]; split_index=split_index)

    if split_index == 0
        train, test = res.dataset_partitions[res.selected_split[2]]
    else
        train, test = res.dataset_partitions[split_index]
    end

    scatter!(ax, y[train], yhat[train]; color=:orange)
    scatter!(ax, y[test], yhat[test]; color=:purple)
    ablines!(ax, 0, 1; color=:grey)

    annotations!(
        ax,
        ["r = $(round(cor(y[train], yhat[train]); digits = 2))"],
        [traincorr_pos],
        textsize = 20,
        color = :orange
    )

    annotations!(
        ax,
        ["r = $(round(cor(y[test], yhat[test]); digits = 2))"],
        [testcorr_pos],
        textsize = 20,
        color = :purple
    )

    return ax

end # end function


## 4.2. Variable imporance plots

### 4.2.1. MDI importance for each bug, averaged over all the train/test splits, for a single model, as a barplot
function singlemodel_avgimportance_barplot!(
    ax::Axis,
    res::T where T <: Resonance.UnivariatePredictor;
    n = 30)

    # Collect the importances
    plot_df = singlemodel_summary_importances(res)

    # Plot barplot
    barplot!(ax, reverse(collect(1:n)), plot_df.AvgImportance[1:n], color = :blue, direction=:x)

    return ax

end # end function

### 4.2.2. MDI importance for each bug, averaged over all the train/test splits, for an ensemble of multiple models,
#   as a barplot colored by the amount of models where the bug was amongst the top N predictor
function multimodel_avgimportance_barplot!(
    figure::Figure,
    ens::UnivariatePredictorEnsemble,
    pos::Tuple{Int64, Int64},
    plot_title::String;
    n_consider = 50,
    n_plot = 50)

    ## Building Figure

    average_importances = multimodel_aggregate_summaryimportances(ens)
    sum_topns = multimodel_aggregate_binarytopns(ens; n = n_consider)
    plot_df = leftjoin(average_importances, sum_topns, on = :Variable)

    dropmissing!(plot_df);
    sort!(plot_df, :AvgMultimodelImportance; rev = true)

    ## Building Figure
    bar_colors = ColorSchemes.viridis.colors[floor.(Int64, collect(range(256, 1, length = 1+maximum(plot_df.SumTopN))))]

    # Build the axis
    ax = Axis(
        figure[pos[1],pos[2]];
        ylabel = "Most important Taxa",
        yticks = (reverse(collect(1:n_plot)), plot_df.Variable[1:n_plot]),
        xlabel = "Average MDI (Gini) Importance over all ($(length(ens.predictors))) models",
        title = plot_title
    )
    
    # Plot barplot
    barplot!(ax, reverse(collect(1:n_plot)), plot_df.AvgMultimodelImportance[1:n_plot], color = bar_colors[plot_df.SumTopN[1:n_plot] .+ 1], direction=:x)
    
    # Plot Legend
    labels = string.( collect(0:1:maximum(plot_df.SumTopN)) )
    elements = [ PolyElement(polycolor = bar_colors[i]) for i in 1:length(labels) ]
    Legend(figure[pos[1]+1,pos[2]], elements, labels, "Predictors for which this taxon is among the top $(n_consider) important factors", orientation=:horizontal)

    return figure

end # end function

## 4.3. Logistic Regressions

### 4.3.1. Classification ground truth versus feature abundance for split-correlation analysis on a specific bug and a specific dataset, with logistic regression line
function singlemodel_logistic_regression!(
    figure::Figure,
    res::UnivariateRandomForestClassifier,
    pos::Tuple{Int64, Int64},
    feature_name::String,
    plot_title::String)

    ## Computation

    X, y = res.inputs_outputs
    regression_df = DataFrame(
        :Input => X[:, feature_name],
        :Output => recode(unwrap.(y), true=>1, false=>0)
    )
    
    fitted_logistic = glm(@formula(Output ~ Input), regression_df, Binomial(), LogitLink())
    
    regressionline_x = collect(range(minimum(regression_df.Input)-Statistics.std(regression_df.Input), maximum(regression_df.Input)+Statistics.std(regression_df.Input), 1000))
    regressionline_y = GLM.predict(fitted_logistic, DataFrame( :Input => regressionline_x))

    ## Plotting

    logistic_colormap=cgrad([:red, :blue], 2) # cgrad and Symbol, mycmap

    ax = Axis(
        figure[pos[1],pos[2]],
        ylabel = "Cognitive Score class",
        xlabel = "Input feature",
        title = "Logistic regression of class as function of a single predictor\n"*plot_title*", "*feature_name,
    )
    
    scatter!(ax, regression_df.Input, regression_df.Output, color = regression_df.Output, colormap = logistic_colormap)
    lines!(ax, regressionline_x, regressionline_y, color = :gray, linestyle = :dash)
    
    return figure

end

### 4.3.2. Classification ground truth versus feature abundance for split-correlation analysis on a set of bugs and a set of data subsets, with logistic regression lines
function multimodel_logistic_regression!(
    figure::Figure,
    ens::UnivariatePredictorEnsemble,
    features_to_consider::Vector{String};
    )

    for j in eachindex(features_to_consider)

        for i in eachindex(ens.predictors)

            singlemodel_logistic_regression!(
                figure,
                ens.predictors[i],
                (j,i),
                features_to_consider[j],
                ens.screen_names[i]
                )
        end # end for each feature to consider
        
    end # end for each predictor
    
    return figure

end