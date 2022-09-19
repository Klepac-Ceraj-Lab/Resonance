using Resonance
using CairoMakie
using ColorSchemes
using Statistics
using MLJ
using DecisionTree
using JLD2
using StableRNGs
using GLM

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

# 1. (Vanja) The same figure I sent in the #echo channel yesterday about brain prediction as univariate outputs, but with left alongside right to make for easier visualization
## Done on specific notebook on the brain branch

# 2. (Vanja) A Heatmap with rows as important bugs across an ensemble of models (classification/regression) and columns as each individual model, colored by the sign and value of the correlation (regular for regression, point-biserial for classification) between the predictor and the output variable.


# 3. (Vanja) A heatmap similar ot the last one, but with the cells divided in half, showing average relative abundance of a taxon on the upper half, and the prevalence of a taxon on the lower half.

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



# 4. (Kevin) Scatters of top bugs abundance x cogScore, colored by age bracket (0-6, 6-12, etc)
# 5. (Kevin) For a given classification model, scatter of top 2 important bugs abundances for the samples included in that model, colored by classification result (same colors as your bar plots)
# 6. (Kevin) X-axis = 4 age brackets (1-4), every bug feature importance for that model on the y-axis, each bug connected by line across 4 x coordinates (not sure how to use color on this one, but might want markers / lines to have transparency if it's too dense.
# 7. (Kevin) Same as above, but rank importance for y axis
# 8. (Kevin) Pairwise correlation of features, ordered using hierarchical clustering.
