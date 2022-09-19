using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
using ProgressMeter

#####
# Loading Data
#####

RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# cogScore classification from taxonomic profiles
JLD2.@load "models/classification_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "models/classification_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "models/classification_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "models/classification_currentCogScores_18to24_fromtaxa_results.jld"
JLD2.@load "models/classification_futureCogScores_allselected_fromtaxa_results.jld"
# cogScore regression from taxonomic profiles
JLD2.@load "models/regression_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_18to24_fromtaxa_results.jld"
JLD2.@load "models/regression_futureCogScores_allselected_fromtaxa_results.jld"

ensemble_classification_models = UnivariatePredictorEnsemble(
    [ "Concurrent, 0-6mo", "Concurrent, 6-12mo", "Concurrent, 12-18mo", "Concurrent, 18-24mo", "Future, all samples" ],
    [ :concurrent_00to06, :concurrent_06to12, :concurrent_12to18, :concurrent_18to24, :future_allsamples ],
    [
        classification_currentCogScores_00to06_fromtaxa_results,
        classification_currentCogScores_06to12_fromtaxa_results,
        classification_currentCogScores_12to18_fromtaxa_results,
        classification_currentCogScores_18to24_fromtaxa_results,
        classification_futureCogScores_allselected_fromtaxa_results
    ]
)

ensemble_regression_models = UnivariatePredictorEnsemble(
    [ "Concurrent, 0-6mo", "Concurrent, 6-12mo", "Concurrent, 12-18mo", "Concurrent, 18-24mo", "Future, all samples" ],
    [ :concurrent_00to06, :concurrent_06to12, :concurrent_12to18, :concurrent_18to24, :future_allsamples ],
    [
        regression_currentCogScores_00to06_fromtaxa_results,
        regression_currentCogScores_06to12_fromtaxa_results,
        regression_currentCogScores_12to18_fromtaxa_results,
        regression_currentCogScores_18to24_fromtaxa_results,
        regression_futureCogScores_allselected_fromtaxa_results
    ]
)

#####
# Functions
#####

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

function singlemodel_singlesplit_correlations(
    res::UnivariateRandomForestClassifier,
    ; split_index=0)

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    model_X, model_y = res.inputs_outputs

    model_trees = res.models[split_index].fitresult[1].trees
    feature_correlations_matrix = feature_split_correlation_analysis(model_X, model_y, model_trees)
    feature_mean_correlations = map( x -> mean(dropnan(dropmissing(x))), eachcol(feature_correlations_matrix))

    correlations_df = DataFrame(
        :Variable => names(model_X),
        :AvgCorrelation => feature_mean_correlations,
    )

    sort!(correlations_df, :AvgCorrelation; rev=true)

    return correlations_df

end

function singlemodel_singlesplit_correlations(
    res::UnivariateRandomForestRegressor,
    ; split_index=0)

    (split_index > length(res.models)) && error("Out of Bounds model/split index selected")
    (split_index == 0) && (split_index = res.selected_split[2])

    model_X, model_y = res.inputs_outputs

    model_trees = res.models[split_index].fitresult.trees
    feature_correlations_matrix = feature_split_correlation_analysis(model_X, model_y, model_trees)
    feature_mean_correlations = map( x -> mean(dropnan(dropmissing(x))), eachcol(feature_correlations_matrix))

    correlations_df = DataFrame(
        :Variable => names(model_X),
        :AvgCorrelation => feature_mean_correlations,
    )

    sort!(correlations_df, :AvgCorrelation; rev=true)

    return correlations_df

end

function singlemodel_allsplits_correlations(res::T where T <: ResonanceUnivariatePredictor)

    singlesplit_correlations = [ singlemodel_singlesplit_correlations(res; split_index = i) for i in 1:length(res.models) ]
    concatenated_correlations_df = reduce((x, y) -> DataFrames.outerjoin(x, y, on = :Variable, makeunique = true), singlesplit_correlations)

    return concatenated_correlations_df

end

function singlemodel_allsplits_correlationsummary(res::T where T <: ResonanceUnivariatePredictor, colname = :AvgCorrelation, fun = nonna_nonmissing_mean)

    concatenated_correlations_df = singlemodel_allsplits_correlations(res)
    summarised_correlations_df = DataFrame(
       :Variable => concatenated_correlations_df.Variable,
       colname => map(fun, eachrow(Matrix(concatenated_correlations_df[:, 2:end])))
    )

    sort!(summarised_correlations_df, colname, rev = true)

    return summarised_correlations_df

end

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





function build_classification_importance_correlations_df(res)

    importance_correlation_dfs = [ process_model_correlations(res, i, model_X, model_y) for i in 1:length(res[:models]) ]
    importance_correlation_dfs = reduce(vcat, importance_correlation_dfs) 

    combined_importance_correlation_df = @chain importance_correlation_dfs begin
        groupby(:Variable)
        combine([:Importance, :MeanCorrelation] .=> non_na_mean; renamecols=false)
    end

    insertcols!(combined_importance_correlation_df, :MeanCorSign => sign.(combined_importance_correlation_df.MeanCorrelation) )
    sort!(combined_importance_correlation_df, :Importance, rev=true)


    return combined_importance_correlation_df

end

function build_classification_summary_plot_df(result_set; result_type="classification", n = 50)

    importance_correlations_dfs = [ build_classification_importance_correlations_df(res) for res in result_set ]
    importance_correlations_dfs = map(x -> insertcols!(x, 3, :TopN => vcat(ones(n), zeros(nrow(x)-n))), importance_correlations_dfs)

    ## Processing Importances

    joined_importances = reduce(
        (x, y) -> outerjoin(x, y, on = :Variable, makeunique = true),
        map(df -> select(df, [:Variable, :Importance]), importance_correlations_dfs)
    )
    average_importances = DataFrame(
        :Variable => joined_importances.Variable,
        :AvgImportance => map(x -> mean(x[.!(ismissing.(x))]), eachrow(Matrix(joined_importances[:, 2:end])))
    )

    ## Processing Correlations

    joined_correlations = reduce(
        (x, y) -> outerjoin(x, y, on = :Variable, makeunique = true),
        map(df -> select(df, [:Variable, :MeanCorrelation]), importance_correlations_dfs)
    )
    average_correlations = DataFrame(
        :Variable => joined_correlations.Variable,
        :AvgCorrelation => map(x -> non_na_mean(x[.!(ismissing.(x))]), eachrow(Matrix(joined_correlations[:, 2:end])))
    )
    n_poscorrs = DataFrame(
        :Variable => joined_correlations.Variable,
        :NPosCorr => map(x -> sum(x .> 0.1), eachrow(Matrix(joined_correlations[:, 2:end])))
    )
    n_zerocorrs = DataFrame(
        :Variable => joined_correlations.Variable,
        :NZeroCorr => map(x -> sum((x .< 0.1) .& (x .> -0.1)), eachrow(Matrix(joined_correlations[:, 2:end])))
    )
    n_negcorrs = DataFrame(
        :Variable => joined_correlations.Variable,
        :NNegCorr => map(x -> sum(x .< -0.1), eachrow(Matrix(joined_correlations[:, 2:end])))
    )
    sd_correlations = DataFrame(
        :Variable => joined_correlations.Variable,
        :SdevCorrelation => map(x -> std(x[.!(isnan.(x))]), eachrow(Matrix(joined_correlations[:, 2:end])))
    )

    joined_topns = reduce(
        (x, y) -> outerjoin(x, y, on = :Variable, makeunique = true),
        map(df -> select(df, [:Variable, :TopN]), importance_correlations_dfs)
    )
    sum_topns = DataFrame(
        :Variable => joined_topns.Variable,
        :SumTopN => map(x -> floor(Int64, sum(x[.!(ismissing.(x))])), eachrow(Matrix(joined_topns[:, 2:end])))
    )

    individual_summary_dfs = [
        average_importances,
        average_correlations,
        sum_topns,
        n_poscorrs,
        n_zerocorrs,
        n_negcorrs,
        sd_correlations,
    ]

    joined_summaries = reduce( (x, y) -> outerjoin(x, y, on = :Variable, makeunique = true), individual_summary_dfs)
    sort!(joined_summaries, :AvgImportance; rev=true)

    return joined_summaries

end

#####
# Computation
#####

RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore classification from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/classification_currentCogScores_18to24_fromtaxa_results.jld"
# concurrent cogScore regression from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_18to24_fromtaxa_results.jld"
# future cogScore classification from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/classification_futureCogScores_allselected_fromtaxa_results.jld"
# future cogScore regression from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/regression_futureCogScores_allselected_fromtaxa_results.jld"

# Axis

classification_result_set = [
    classification_currentCogScores_00to06_fromtaxa_results,
    classification_currentCogScores_06to12_fromtaxa_results,
    classification_currentCogScores_12to18_fromtaxa_results,
    classification_currentCogScores_18to24_fromtaxa_results,
    classification_futureCogScores_allselected_fromtaxa_results,    
]

regression_result_set = [
    regression_currentCogScores_00to06_fromtaxa_results,
    regression_currentCogScores_06to12_fromtaxa_results,
    regression_currentCogScores_12to18_fromtaxa_results,
    regression_currentCogScores_18to24_fromtaxa_results,
    regression_futureCogScores_allselected_fromtaxa_results,    
]

classification_summary_plot = build_classification_summary_plot_df(classification_result_set)

# combined_importance_correlation_df = build_classification_importance_correlations_df(regression_currentCogScores_00to06_fromtaxa_results)

#####
# Plotting
#####

