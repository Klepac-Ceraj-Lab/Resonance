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

function descript_inputs(res::T where T <: ResonancePredictor)

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
        :MeanSdevConsideringZeros => bug_std_withzeros,
        :MeanSdevExcludingZeros => bug_std_nonzeros
    )

    return descript_df

end

function singlemodel_summary_prevalences(res::T where T <: ResonanceUnivariatePredictor, colname = :Prevalence)

    descript_df = descript_inputs(res)
    prevalence_df = DataFrame(
        :Variable => descript_df.Variable,
        colname => descript_df.Prevalence
    )

    return prevalence_df

end

function singlemodel_summary_abundances(res::T where T <: ResonanceUnivariatePredictor, colname = :MeanAbundance; exclude_zeros=true)

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

function singlemodel_summary_sdevs(res::T where T <: ResonanceUnivariatePredictor, colname = :MeanSdev; exclude_zeros=true)

    descript_df = descript_inputs(res)

    if exclude_zeros
        sdevs_df = DataFrame(
            :Variable => descript_df.Variable,
            colname => descript_df.MeanSdevExcludingZeros
        )
    else
        sdevs_df = DataFrame(
            :Variable => descript_df.Variable,
            colname => descript_df.MeanSdevConsideringZeros
        )
    end    

    sdevs_df[isnan.(sdevs_df[:, colname]), colname] .= 0.0

    return sdevs_df
    
end

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

#####
# Computation
#####

## Classification

classification_individual_abundances = multimodel_individual_abundances(ensemble_classification_models)
classification_individual_prevalences = multimodel_individual_prevalences(ensemble_classification_models)
classification_individual_correlations = multimodel_individual_correlations(ensemble_classification_models)
classification_aggregate_importances = get_multimodel_aggregate_summaryimportances(ensemble_classification_models)
classification_aggregate_importances.Order = 1:nrow(classification_aggregate_importances)

plot_aggregate_df = DataFrames.outerjoin(classification_aggregate_importances, classification_individual_correlations, on = :Variable, makeunique = true)
plot_aggregate_df = DataFrames.outerjoin(plot_aggregate_df, classification_individual_prevalences, on = :Variable, makeunique = true)
plot_aggregate_df = DataFrames.outerjoin(plot_aggregate_df, classification_individual_abundances, on = :Variable, makeunique = true)
sort!(plot_aggregate_df, :Order)


#####
# Plotting
#####

ens = ensemble_classification_models
n_plot = 50
heatmap_figure = Figure(resolution = (1800, 1800) )

## Correlation

correlation_colormap = cgrad(
    [:red, :white, :blue],
    [0.0, 0.5, 1.0]
)

ax_correlation = Axis(
    heatmap_figure[1,1];
    ylabel = "Most important Taxa",
    yticks = (collect(1:n_plot), plot_aggregate_df.Variable[1:n_plot]),
    xlabel = "Model",
    xticks = (collect(1:5), ens.screen_names),
    title = "Average Point-biserial correlation through the $(length(ens.predictors)) models",
    xticklabelrotation = pi/2,
    yreversed = true
)

corr_mat = permutedims(Matrix(plot_aggregate_df[1:n_plot, 4:8]))

correlation_hm = heatmap!(ax_correlation, corr_mat; colormap=correlation_colormap, colorrange=(-0.8, 0.8), highclip = :blue, lowclip = :red)

for ci in CartesianIndices(corr_mat)
    text!(ax_correlation,
    string(round(corr_mat[ci], digits=2)),
    ; position=(ci[1],ci[2]),
    align=(:center, :center)
    )
end

## Prevalence

prevalence_colormap = :viridis

ax_prevalence = Axis(
    heatmap_figure[1,2];
    #ylabel = "Most important Taxa",
    yticks = (collect(1:n_plot), plot_aggregate_df.Variable[1:n_plot]),
    yticksvisible = false,
    yticklabelsvisible = false,
    xlabel = "Model",
    xticks = (collect(1:5), ens.screen_names),
    title = "Average feature Prevalence through the $(length(ens.predictors)) models",
    xticklabelrotation = pi/2,
    yreversed = true
)

prevalence_mat = permutedims(Matrix(plot_aggregate_df[1:n_plot, 9:13]))

prevalence_hm = heatmap!(ax_prevalence, prevalence_mat; colormap=prevalence_colormap)

for ci in CartesianIndices(prevalence_mat)
    text!(ax_prevalence,
    string(round(prevalence_mat[ci], digits=4)),
    ; position=(ci[1],ci[2]),
    align=(:center, :center)
    )
end

## Prevalence

abundance_colormap = :viridis

ax_abundance = Axis(
    heatmap_figure[1,3];
    #ylabel = "Most important Taxa",
    yticks = (collect(1:n_plot), plot_aggregate_df.Variable[1:n_plot]),
    yticksvisible = false,
    yticklabelsvisible = false,
    xlabel = "Model",
    xticks = (collect(1:5), ens.screen_names),
    title = "Average feature Abundance through the $(length(ens.predictors)) models",
    xticklabelrotation = pi/2,
    yreversed = true
)

abundance_mat = permutedims(Matrix(plot_aggregate_df[1:n_plot, 14:18]))

abundance_hm = heatmap!(ax_abundance, abundance_mat; colormap=abundance_colormap)

for ci in CartesianIndices(abundance_mat)
    text!(ax_abundance,
    string(round(abundance_mat[ci], digits=4)),
    ; position=(ci[1],ci[2]),
    align=(:center, :center)
    )
end

Colorbar(heatmap_figure[2, 1], correlation_hm, label = "Point-Biserial correlation", vertical = false)
Colorbar(heatmap_figure[2, 2], prevalence_hm, label = "Prevalence", vertical = false)
Colorbar(heatmap_figure[2, 3], abundance_hm, label = "Abundance", vertical = false)

heatmap_figure

save("figures/correlation_prevalence_abundance_classification_heatmap.png", heatmap_figure)