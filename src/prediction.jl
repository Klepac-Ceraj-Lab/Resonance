########################################
# 0. Hierarchical Structures and Types #
########################################

abstract type Prediction end
struct Classification <: Prediction end
struct Regression <: Prediction end

abstract type ResonancePredictor end
abstract type ResonanceUnivariatePredictor <: ResonancePredictor end
abstract type ResonanceMultivariatePredictor <: ResonancePredictor end

mutable struct UnivariateRandomForestClassifier <: ResonanceUnivariatePredictor
    name::String
    inputs_outputs::Tuple{DataFrame, CategoricalArrays.CategoricalVector{Bool, UInt32, Bool, CategoricalArrays.CategoricalValue{Bool, UInt32}, Union{}}}
    n_splits::Int64
    dataset_partitions::Vector{Tuple{Vector{Int64}, Vector{Int64}}}
    models::Vector{Machine}
    selected_split::Tuple{Float64, Int64}
    train_accuracies::Vector{Float64}
    test_accuracies::Vector{Float64}
end

mutable struct UnivariateRandomForestRegressor <: ResonanceUnivariatePredictor
    name::String
    inputs_outputs::Tuple{DataFrame, Vector{Float64}}
    n_splits::Int64
    dataset_partitions::Vector{Tuple{Vector{Int64}, Vector{Int64}}}
    models::Vector{Machine}
    slope_corrections::Vector{}
    selected_split::Tuple{Float64, Int64}
    train_maes::Vector{Float64}
    test_maes::Vector{Float64}
    train_mapes::Vector{Float64}
    test_mapes::Vector{Float64}
    train_cor::Vector{Float64}
    test_cor::Vector{Float64}
end

mutable struct UnivariatePredictorEnsemble
    screen_names::Vector{String}
    col_names::Vector{Symbol}
    predictors::Vector{T where T <: ResonanceUnivariatePredictor}
end

########################################
# 1. Preprocessing/Wrangling Functions #
########################################

non_na_mean(vv) = mean(vv[.!(isnan.(vv))])

function myxor(a::Float64, b::Float64)
    a == b && (return a)
    a > b && (return a)
    a < b && (return b)
end

function tryparsecol(T, col)
    newcol = tryparse.(T, col)
    return any(isnothing, newcol) ? col : newcol
end

function build_metadata_prediction_df(base_df, inputs::Vector{Symbol}, targets::Vector{Symbol})

    subjects = unique(base_df.subject)
    has_stool = findall(.!(ismissing.(base_df.sample)))
    has_cogscore = findall(.!(ismissing.(base_df.cogScore)))

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
    
    select!(prediction_df, [:subject, :ageMonths, :timepoint, :futureAgeMonths, :futureTimepoint, :ageMonthsDelta, :timepointDelta, targets..., Symbol.("future" .* uppercasefirst.(String.(targets)))..., inputs...])
    sort!(prediction_df, [:subject, :timepoint, :futureTimepoint])
    
    return prediction_df

end

function compute_tietjenmoore(data, k, n)
    r_all = abs.(data .- mean(data))
    filteredData = data[sortperm(r_all)][1:(n-k)]
    ksub = (filteredData .- mean(filteredData))

    return( sum(ksub .^ 2) / sum(r_all .^ 2) )
end

function test_tietjenmoore(dataSeries,k, n, alpha, sim_trials)
    ek = compute_tietjenmoore(dataSeries,k, n)
    simulation = [ compute_tietjenmoore(randn(length(dataSeries)), k, n) for i in 1:sim_trials ]
    Talpha=quantile(simulation,alpha)

    return(ek, Talpha)
end

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

function unstack_techreplicates(raw_data, retain_first=true)

    unstackedDf = retain_first ?
    begin 
        @chain raw_data begin
            groupby([:subject, :timepoint, :variable])
            combine(:value=>first)
            unstack(:variable, :value_first)
        end #end @chain
    end : #end if_true
    unstack(raw_data, :variable, :value; allowduplicates=true)

    return unstackedDf
end #end function

####################################
# 2. Training/Processing Functions #
####################################

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree

function train_randomforest(
    type::Regression,
    ref_name::String,
    original_df,
    input_cols,
    output_col;
    n_splits = 2,
    data_preparation_fun = identity,
    tuning_space = (
        maxnodes_range = [ -1 ] ,
        nodesize_range = [ 0 ],
        sampsize_range = [ 0.7 ],
        mtry_range = [ 5 ],
        ntrees_range = [ 10 ]
    ),
    split_proportion=0.75,
    train_rng=Random.GLOBAL_RNG
    )

    # ## 1. Subsetting Data
    prediction_df = data_preparation_fun(original_df)
    
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
    trial_slopecorrections = Vector{T where T <: RegressionModel}(undef, n_splits)
    trial_train_maes = repeat([Inf], n_splits)
    trial_train_mapes = repeat([Inf], n_splits)
    trial_train_cors = repeat([-1.0], n_splits)
    trial_test_maes = repeat([Inf], n_splits)
    trial_test_mapes = repeat([Inf], n_splits)
    trial_test_cors = repeat([-1.0], n_splits)

    # ## 5. Tuning/training loop
    @info "Performing $(n_splits) different train/test splits and tuning $(length(tuning_grid)) different hyperparmeter combinations\nfor the $(nrow(prediction_df)) samples."

    for this_trial in 1:n_splits

        ## 5.1. Splitting training data between train and test
        Random.seed!(train_rng, this_trial)
        train, test = partition(eachindex(1:nrow(X)), split_proportion, shuffle=true, rng=train_rng)
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
            ## 5.2.2.1. Fit slope correction on training data

            train_y_hat = MLJ.predict(rf_machine, X[train, :])
            slope_correction_df = DataFrame([ train_y_hat, y[train] ], [:yhat, :y])
            slope_correction = lm(@formula(y ~ yhat), slope_correction_df)
            train_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => train_y_hat))
            train_mape = mape(train_y_hat, y[train])

            ## 5.2.3. Benchmark model on independent testing data
            test_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => MLJ.predict(rf_machine, X[test, :])))
            test_mape = mape(test_y_hat, y[test])

            ## 5.2.4. Compare benchmark results with previous best model
            test_mape < trial_test_mapes[this_trial] ? begin
                trial_train_maes[this_trial] = mae(train_y_hat, y[train])
                trial_train_mapes[this_trial] = train_mape
                trial_train_cors[this_trial] = Statistics.cor(train_y_hat, y[train])
                trial_test_maes[this_trial] = mae(test_y_hat, y[test])
                trial_test_mapes[this_trial] = test_mape
                trial_test_cors[this_trial] = Statistics.cor(test_y_hat, y[test])
                trial_machines[this_trial] = deepcopy(rf_machine)
                trial_slopecorrections[this_trial] = slope_correction
            end : continue

        end # end for i in 1:length(tuning_grid) - each set of hyperparameters

    end # end for this_trial - each different train/test split

    # ## 6. Returning optimization results
    results = UnivariateRandomForestRegressor(
        ref_name,                   #name::String
        (X,y),                      #inputs_outputs::Tuple{DataFrame, Vector{Float64}}
        n_splits,                   #n_splits::Int64
        trial_partitions,           #dataset_partitions::Vector{Tuple{Vector{Int64}, Vector{Int64}}}
        trial_machines,             #models::Vector{Machine}
        trial_slopecorrections,     #slope_correction::Vector{}
        findmin(trial_test_mapes),  #selected_split::Tuple{Float64, Int64}
        trial_train_maes,           #train_maes::Vector{Float64}
        trial_test_maes,            #test_maes::Vector{Float64}
        trial_train_mapes,          #train_mapes::Vector{Float64}
        trial_test_mapes,           #test_mapes::Vector{Float64}
        trial_train_cors,           #train_cor::Vector{Float64}
        trial_test_cors,            #test_cor::Vector{Float64}
    )

    @info "Done!"
    return results

end # end function

#########################################
# 3. Post-processing/Analysis Functions #
#########################################

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

    return merits
end

function report_merits(res::UnivariateRandomForestRegressor)
    
    merits = DataFrame(
        :Split => collect(1:res.n_splits),
        :Train_MAE => res.train_maes,
        :Test_MAE => res.test_maes,
        :Train_MAPE => res.train_mapes,
        :Test_MAPE => res.test_mapes,
        :Train_COR => res.train_cor,
        :Test_COR => res.test_cor
        )

    push!(merits, Dict(
        :Split => 0,
        :Train_MAE => mean(merits.Train_MAE),
        :Test_MAE => mean(merits.Test_MAE),
        :Train_MAPE => mean(merits.Train_MAPE),
        :Test_MAPE => mean(merits.Test_MAPE),
        :Train_COR => mean(merits.Train_COR),
        :Test_COR => mean(merits.Test_COR)
    ))
    
    return merits
end

function get_singlemodel_singlesplit_importance(res::UnivariateRandomForestClassifier; split_index=0)

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    importance_df = DataFrame(
        :Variable => names(res.inputs_outputs[1]),
        :Importance => DecisionTree.impurity_importance(res.models[split_index].fitresult[1])
    )

    sort!(importance_df, :Importance, rev = true)

    return importance_df

end

function get_singlemodel_singlesplit_importance(res::UnivariateRandomForestRegressor; split_index = 0 )

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    importance_df = DataFrame(
        :Variable => names(res.inputs_outputs[1]),
        :Importance => DecisionTree.impurity_importance(res.models[split_index].fitresult)
    )

    sort!(importance_df, :Importance, rev = true)

    return importance_df
    
end

function get_singlemodel_allsplits_importances(res::T where T <: ResonanceUnivariatePredictor)

    singlesplit_importances = [ get_singlemodel_singlesplit_importance(res; split_index = i) for i in 1:length(res.models) ]
    concatenated_importances_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlesplit_importances)

    return concatenated_importances_df

end

function get_singlemodel_summary_importances(res::T where T <: ResonanceUnivariatePredictor, colname = :AvgImportance, fun = non_na_mean)

    concatenated_importances_df = get_singlemodel_allsplits_importances(res)
    summarised_importances_df = DataFrame(
       :Variable => concatenated_importances_df.Variable,
       colname => map(fun, eachrow(Matrix(concatenated_importances_df[:, 2:end])))
    )

    sort!(summarised_importances_df, colname, rev = true)

    return summarised_importances_df

end

function get_singlemodel_binarytopn_importances(res::T where T <: ResonanceUnivariatePredictor, importance_colname = :AvgImportance, topn_colname = :TopN, fun = non_na_mean; n=50)

    summarised_importances_df = @chain get_singlemodel_summary_importances(res, importance_colname, fun) begin
        insertcols!( _ , 2, topn_colname => vcat(ones(Int64, n), zeros(Int64, nrow( _ )-n)) )
        select!( Not(importance_colname) )
    end
    
    return(summarised_importances_df)

end

function get_multimodel_individual_summaryimportances(ens::UnivariatePredictorEnsemble, col_prefix = "meanImportance_", fun = non_na_mean)

    singlemodel_summaries = [ 
        get_singlemodel_summary_importances(
            ens.predictors[i],
            Symbol(col_prefix * string(ens.col_names[i])),
            fun
        ) for i in 1:length(ens.predictors)
    ]

    concatenated_summaries_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlemodel_summaries)

    return concatenated_summaries_df

end

function get_multimodel_individual_binarytopns(ens::UnivariatePredictorEnsemble, col_prefix = "topN_", fun = non_na_mean; n=50)

    singlemodel_summaries = [ 
        get_singlemodel_binarytopn_importances(
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

function get_multimodel_aggregate_summaryimportances(
    ens::UnivariatePredictorEnsemble,
    singlemodel_col_prefix = "meanImportance_",
    aggregate_colname = :AvgMultimodelImportance,
    singlemodel_summary_fun = non_na_mean,
    multimodel_summary_fun = non_na_mean)

    concatenated_summaries_df = get_multimodel_individual_summaryimportances(ens, singlemodel_col_prefix, singlemodel_summary_fun)

    summarised_importances_df = DataFrame(
        :Variable => concatenated_summaries_df.Variable,
        aggregate_colname => map(multimodel_summary_fun, eachrow(Matrix(concatenated_summaries_df[:, 2:end])))
    )

    sort!(summarised_importances_df, aggregate_colname, rev = true)

    return summarised_importances_df

end

function get_multimodel_aggregate_binarytopns(
    ens::UnivariatePredictorEnsemble,
    singlemodel_col_prefix = "topN_",
    aggregate_colname = :SumTopN,
    singlemodel_summary_fun = non_na_mean,
    multimodel_summary_fun = sum;
    n=30)

    concatenated_summaries_df = get_multimodel_individual_binarytopns(ens, singlemodel_col_prefix, singlemodel_summary_fun; n=n)

    summarised_importances_df = DataFrame(
        :Variable => concatenated_summaries_df.Variable,
        aggregate_colname => map(multimodel_summary_fun, eachrow(Matrix(concatenated_summaries_df[:, 2:end])))
    )

    sort!(summarised_importances_df, aggregate_colname, rev = true)

    return summarised_importances_df

end

function build_confusion_matrix(classification_results::Dict, trial::Int)

    y = classification_results[:inputs_outputs][2][classification_results[:dataset_partitions][trial][2]]
    yhat = MLJ.predict_mode(classification_results[:models][trial], classification_results[:inputs_outputs][1][classification_results[:dataset_partitions][trial][2],:])
    confmat = ConfusionMatrix()(yhat, y)

    return confmat
end

function average_confusion_matrix(classification_results::Dict)

    matrix_vector = [ build_confusion_matrix(classification_results, i).mat for i in 1:classification_results[:n_trials] ]
    concatenated_matrix = vcat([ vec(mat)' for mat in matrix_vector ]...) # Column order is TN, FP, FN, TP !
    averages = [ sum(concatenated_matrix[:,i])/sum(concatenated_matrix) for i in 1:4]

    return averages
end

function confmatrix2barplot(classification_results::Dict)

    confmat_inputs = (
        x = [1, 2, 2, 1], # Column order is TN, FP, FN, TP !
        value = average_confusion_matrix(classification_results),
        grp = [1, 2, 2, 1],
        color = [1, 2, 3, 4]
    )

    return confmat_inputs
end

function predict_bestsplit(regression_results::UnivariateRandomForestRegressor)

    selected_split = regression_results.selected_split[2]
    X, y = regression_results.inputs_outputs
    train, test = regression_results.dataset_partitions[selected_split]
    slope_correction = regression_results.slope_corrections[selected_split]
    selected_machine = regression_results.models.[selected_split]
    yhat = GLM.predict(slope_correction, DataFrame( :yhat => MLJ.predict(selected_machine, X)))

    return y, yhat, train, test
end

function predict_bestsplit(regression_results::UnivariateRandomForestClassifier)

    selected_split = regression_results.selected_split[2]
    X, y = regression_results.inputs_outputs
    train, test = regression_results.dataset_partitions[selected_split]
    selected_machine = regression_results.models.[selected_split]
    yhat = MLJ.predict(selected_machine, X)

    return y, yhat, train, test
end

#########################
# 4. Plotting functions #
#########################

function singlemodel_avgimportance_barplot!(
    figure::Figure,
    res::T where T <: ResonanceUnivariatePredictor,
    pos::Tuple{Int64, Int64},
    plot_title::String;
    n = 30)

    # Collect the importances
    plot_df = get_singlemodel_summary_importances(res)

    # Build the axis
    ax1_1 = Axis(
        figure[pos[1],pos[2]];
        ylabel = "Most important Taxa",
        yticks = (reverse(collect(1:n)), plot_df.Variable[1:n]),
        #yticklabelrotation=-pi/4,
        xlabel = "Mean Decrease in Impurity (Gini) Importance",
        title = plot_title
    )

    # Plot barplot
    barplot!(ax1_1, reverse(collect(1:n)), plot_df.AvgImportance[1:n], color = :blue, direction=:x)

    return figure

end # end function

function singlemodel_merit_barplot!(
    figure::Figure,
    res::UnivariateRandomForestRegressor,
    pos::Tuple{Int64, Int64},
    plot_title::String)

    # Build the axis
    ax = Axis(
        figure[pos[1],pos[2]];
        xlabel = "Ground Truth",
        ylabel = "Prediction",
        title = plot_title
    )

    # Plot barplot
    y, yhat, train, test = predict_bestsplit(res)
    scatter!(ax, y[train], yhat[train]; color=:orange)
    scatter!(ax, y[test], yhat[test]; color=:purple)
    ablines!(ax, 0, 1; color=:grey)
    annotations!( ax, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(1.1*min(y), 0.9*max(yhat))], textsize = 20)

    return figure

end # end function

function singlemodel_merit_scatterplot!(
    figure::Figure,
    res::UnivariateRandomForestRegressor,
    pos::Tuple{Int64, Int64},
    plot_title::String)

    # Build the axis
    ax = Axis(
        figure[pos[1],pos[2]];
        xlabel = "Ground Truth",
        ylabel = "Prediction",
        title = plot_title
    )

    # Plot barplot
    y, yhat, train, test = predict_bestsplit(res)
    scatter!(ax, y[train], yhat[train]; color=:orange)
    scatter!(ax, y[test], yhat[test]; color=:purple)
    ablines!(ax, 0, 1; color=:grey)
    annotations!( ax, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(1.1*min(y), 0.9*max(yhat))], textsize = 20)

    return figure

end # end function

function multimodel_avgimportance_barplot!(
    figure::Figure,
    ens::UnivariatePredictorEnsemble,
    pos::Tuple{Int64, Int64},
    plot_title::String;
    n_consider = 50,
    n_plot = 50)

    ## Building Figure

    average_importances = get_multimodel_aggregate_summaryimportances(ens)
    sum_topns = get_multimodel_aggregate_binarytopns(ens; n = n_consider)
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