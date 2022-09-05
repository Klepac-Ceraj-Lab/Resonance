using Resonance
using CairoMakie
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM
using ProgressMeter

non_na_mean(a) = mean(a[.!(isnan.(a))])
non_missing_mean(a) = mean(a[.!(ismissing.(a))])

RandomForestClassifier = MLJ.@load RandomForestClassifier pkg=DecisionTree
RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree
# concurrent cogScore classification from taxonomic profiles
JLD2.@load "/home/guilherme/Documents/models/regression_currentCogScores_00to06_fromtaxa_results.jld"
result_set = regression_currentCogScores_00to06_fromtaxa_results

#m = classification_currentCogScores_00to06_fromtaxa_results[:models][1].fitresult[1]
m = result_set[:models][1].fitresult
model_trees = m.trees
model_X = result_set[:inputs_outputs][1]
model_y = result_set[:inputs_outputs][2]

function biserial_correlation( ## Borrowed from m3g/XCorrelation.jl (!!!)
	discreteVariable,
	continuousArray :: Array{Float64})

	## Sanity check
	
	discreteArray = discreteVariable .* 1.

	length(discreteArray) == length(continuousArray) || throw(DimensionMismatch("Size of discrete Array not equal to size of continuous Array"))

	## Computation

	Ntrue = sum(discreteArray .== 1)
	Nfalse = sum(discreteArray .== 0)

	SmeanTrue = Statistics.mean(continuousArray[discreteArray .== 1])
	SmeanFalse = Statistics.mean(continuousArray[discreteArray .== 0])

	sdev = Statistics.std(continuousArray, corrected=false)

	bisCorr = ((SmeanTrue - SmeanFalse)/sdev)*sqrt(Ntrue*Nfalse/(length(discreteArray)^2))

	(bisCorr == NaN) && (bisCorr = -1.)

	return( bisCorr )

end

function feature_split_correlation_analysis(X, y, trees)

    feature_correlations_matrix = Matrix{Union{Float64, Missing}}(undef, length(trees), ncol(X))

    @showprogress for tree_idx in 1:length(trees) # For each tree

        sample_feature_rel = zeros(Int64, nrow(X), ncol(X))
        feature_split_rel = Matrix{Union{Int64, Missing}}(undef, nrow(X), ncol(X))
        decision_split_feature = zeros(Int64, nrow(X))

        for sample_idx in 1:nrow(X) # For each sample, make a tree pass.

            this_tree = deepcopy(trees[tree_idx])

            while( !( typeof(this_tree) <: DecisionTree.Leaf ) )
            # while( !( ( typeof(this_tree.left) <: DecisionTree.Leaf ) & ( typeof(this_tree.right) <: DecisionTree.Leaf ) ) )

                if this_tree.featid == 0
                    this_tree = this_tree.left                
                elseif X[sample_idx, this_tree.featid] < this_tree.featval
                    sample_feature_rel[sample_idx, this_tree.featid] = 1
                    feature_split_rel[sample_idx, this_tree.featid] = 0
                    decision_split_feature[sample_idx] = this_tree.featid
                    this_tree = this_tree.left                
                else
                    sample_feature_rel[sample_idx, this_tree.featid] = 1
                    feature_split_rel[sample_idx, this_tree.featid] = 1
                    this_tree = this_tree.right
                end # end state update
        
            end # end while still a node

        end # end for sample_idx

        for feat_idx in 1:ncol(X)

            samples_to_consider = Bool.(sample_feature_rel[:, feat_idx])

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

    return feature_correlations_matrix

end

feature_correlations_matrix = feature_split_correlation_analysis(model_X, model_y, model_trees)
feature_mean_correlations = map(non_missing_mean, eachcol(feature_correlations_matrix))

importance_correlations_df = DataFrame(
    :Variable => names(X),
    :Importance => DecisionTree.impurity_importance(m),
    :MeanCorrelation => feature_mean_correlations,
    :MeanCorSign => sign.(feature_mean_correlations)
)

sort!(importance_correlations_df, :Importance, rev=true)

importance_correlations_df = importance_correlations_df[1:30, :]


#####
# plotting
#####

