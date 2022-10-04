#####
# Notebook DXA02 - Regression of current AgeMonths above/below average from current taxonomic data
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

mdata = @chain Resonance.load(Metadata()) begin
    select( [ :subject, :timepoint, :ageMonths, :cogScore, :omni ] )
    dropmissing( [ :omni ] )
end

taxa = @chain Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata) begin
    filter( f-> taxrank(f) == :species, _ )
end

taxa_df = DataFrame([eachcol(collect(taxa.abundances'))...], map( x -> string(x), taxa.features), copycols=true)
taxa_df = taxa_df[ !, sortperm(names(taxa_df)) ]
rename!(taxa_df, [ names(taxa_df) .=> replace.(names(taxa_df), "s__" => "")]...)
taxa_df.Collinsella_massiliensis = myxor.(taxa_df.Collinsella_massiliensis, taxa_df."[Collinsella]_massiliensis")
select!(taxa_df, Not("[Collinsella]_massiliensis"))
insertcols!(taxa_df, 1, :sample => collect(keys(taxa.sidx))[collect(values(taxa.sidx))])

cogscore_taxa_df = leftjoin(mdata, taxa_df, on = :omni => :sample, matchmissing=:error) |>  y -> rename!(y, :omni => :sample) |>  y -> sort(y, [ :subject, :timepoint ]);

#####
# Training and saving models - uncomment to repeat
#####
# RandomForestRegressor= MLJ.@load RandomForestRegressor pkg=DecisionTree
# tuning_space = (
#         maxnodes_range = collect(1:3:15) ,
#         nodesize_range = collect(1:3:20),
#         sampsize_range = [0.5, 0.6, 0.7],
#         mtry_range = collect(5:10:100),
#         ntrees_range = [100, 300, 500]
#     )
# ## 0 to 6 months
# regression_currentAgeMonths_00to06_fromtaxa_results = train_randomforest(
#     Resonance.Regression(),
#     "regression_currentAgeMonths_00to06_fromtaxa_results",
#     cogscore_taxa_df,
#     x -> dropmissing(filter_age_bracket(x, 0.0, 6.0)),
#     6:554,
#     :ageMonths;
#     n_splits = 10,
#     tuning_space = tuning_space,
#     train_rng = ml_rng
# )
# JLD2.@save "models/regression_currentAgeMonths_00to06_fromtaxa_results.jld" regression_currentAgeMonths_00to06_fromtaxa_results
# ## 6 to 12 months
# regression_currentAgeMonths_06to12_fromtaxa_results = train_randomforest(
#     Resonance.Regression(),
#     "regression_currentAgeMonths_06to12_fromtaxa_results",
#     cogscore_taxa_df,
#     x -> dropmissing(filter_age_bracket(x, 6.0, 12.0)),
#     6:554,
#     :ageMonths;
#     n_splits = 10,
#     tuning_space = tuning_space,
#     train_rng = ml_rng
# )
# JLD2.@save "models/regression_currentAgeMonths_06to12_fromtaxa_results.jld" regression_currentAgeMonths_06to12_fromtaxa_results
# ## 12 to 18 months
# regression_currentAgeMonths_12to18_fromtaxa_results = train_randomforest(
#     Resonance.Regression(),
#     "regression_currentAgeMonths_12to18_fromtaxa_results",
#     cogscore_taxa_df,
#     x -> dropmissing(filter_age_bracket(x, 12.0, 18.0)),
#     6:554,
#     :ageMonths;
#     n_splits = 10,
#     tuning_space = tuning_space,
#     train_rng = ml_rng
# )
# JLD2.@save "models/regression_currentAgeMonths_12to18_fromtaxa_results.jld" regression_currentAgeMonths_12to18_fromtaxa_results
# ## 18 to 24 months
# regression_currentAgeMonths_18to24_fromtaxa_results = train_randomforest(
#     Resonance.Regression(),
#     "regression_currentAgeMonths_18to24_fromtaxa_results",
#     cogscore_taxa_df,
#     x -> dropmissing(filter_age_bracket(x, 18.0, 24.0)),
#     6:554,
#     :ageMonths;
#     n_splits = 10,
#     tuning_space = tuning_space,
#     train_rng = ml_rng
# )
# JLD2.@save "models/regression_currentAgeMonths_18to24_fromtaxa_results.jld" regression_currentAgeMonths_18to24_fromtaxa_results
# ## 00 to 120 months
# regression_currentAgeMonths_0to120_fromtaxa_results = train_randomforest(
#     Resonance.Regression(),
#     "regression_currentAgeMonths_0to120_fromtaxa_results",
#     cogscore_taxa_df,
#     x -> dropmissing(filter_age_bracket(x, 0.0, 120.0)),
#     6:554,
#     :ageMonths;
#     n_splits = 10,
#     tuning_space = tuning_space,
#     train_rng = ml_rng
# )
# JLD2.@save "models/regression_currentAgeMonths_0to120_fromtaxa_results.jld" regression_currentAgeMonths_0to120_fromtaxa_results

#####
# Analyzing models
#####

# cogScore regression from taxonomic profiles
JLD2.@load "models/regression_currentAgeMonths_00to06_fromtaxa_results.jld"
JLD2.@load "models/regression_currentAgeMonths_06to12_fromtaxa_results.jld"
JLD2.@load "models/regression_currentAgeMonths_12to18_fromtaxa_results.jld"
JLD2.@load "models/regression_currentAgeMonths_18to24_fromtaxa_results.jld"
JLD2.@load "models/regression_currentAgeMonths_0to120_fromtaxa_results.jld"

ensemble_regression_models = UnivariatePredictorEnsemble(
    [ "Concurrent, 0-6mo", "Concurrent, 6-12mo", "Concurrent, 12-18mo", "Concurrent, 18-24mo", "Concurrent, 0-120mo" ],
    [ :concurrent_00to06, :concurrent_06to12, :concurrent_12to18, :concurrent_18to24, :concurrent_0to120 ],
    [
        regression_currentAgeMonths_00to06_fromtaxa_results,
        regression_currentAgeMonths_06to12_fromtaxa_results,
        regression_currentAgeMonths_12to18_fromtaxa_results,
        regression_currentAgeMonths_18to24_fromtaxa_results,
        regression_currentAgeMonths_0to120_fromtaxa_results
    ]
)

## Plot summary figure for Regression models
figure = Figure( ; resolution = (1000, 1000) )
multimodel_avgimportance_barplot!(
    figure,
    ensemble_regression_models,
    (1,1), "Summary of variable importances throughout all Age regressions models"
)
save("figures/importance_plots/Summarybarplot_Ages_00to24_regression_fromtaxa.png", figure)

## Plot regression scatterplots

figure = Figure(; resolution = (1500, 600))

singlemodel_merit_scatterplot!(figure, regression_currentAgeMonths_00to06_fromtaxa_results, (1,1), "0 to 6 months (n = 73)")
singlemodel_merit_scatterplot!(figure, regression_currentAgeMonths_06to12_fromtaxa_results, (1,2), "6 to 12 months (n = 61)")
singlemodel_merit_scatterplot!(figure, regression_currentAgeMonths_12to18_fromtaxa_results, (1,3), "12 to 18 months (n = 39)")
singlemodel_merit_scatterplot!(figure, regression_currentAgeMonths_18to24_fromtaxa_results, (1,4), "18 to 24 months (n = 50)")
singlemodel_merit_scatterplot!(figure, regression_currentAgeMonths_0to120_fromtaxa_results, (1,5), "0 to 120 months")

labels = ["Training samples", "Test samples"]
elements = [PolyElement(polycolor = el) for el in [:orange, :purple]]
Legend(figure[2,1:5], elements, labels, "Sample set", nbanks = 1, orientation=:horizontal, valign = :top)

Label(figure[1, :, Top()], "Regression results for prediction of ageMonths from present taxonomic profiles for all\nqualifying samples divided in 4 6-month brackets", valign = :bottom,
    # font = "TeX Gyre Heros Bold",
    padding = (0, 50, 40, 0))

save("figures/AgeMonths_merit_scatterplots.png", figure)

#####
# Exclusion analysis
#####

JLD2.@load "models/regression_currentCogScores_00to06_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_06to12_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_12to18_fromtaxa_results.jld"
JLD2.@load "models/regression_currentCogScores_18to24_fromtaxa_results.jld"

ensemble_agemonths_models = UnivariatePredictorEnsemble(
    [ "Concurrent, 0-6mo", "Concurrent, 6-12mo", "Concurrent, 12-18mo", "Concurrent, 18-24mo", "Future, all samples" ],
    [ :concurrent_00to06, :concurrent_06to12, :concurrent_12to18, :concurrent_18to24, :future_allsamples ],
    [
        regression_currentAgeMonths_00to06_fromtaxa_results,
        regression_currentAgeMonths_06to12_fromtaxa_results,
        regression_currentAgeMonths_12to18_fromtaxa_results,
        regression_currentAgeMonths_18to24_fromtaxa_results
    ]
)

ensemble_cogscore_models = UnivariatePredictorEnsemble(
    [ "Concurrent, 0-6mo", "Concurrent, 6-12mo", "Concurrent, 12-18mo", "Concurrent, 18-24mo", "Future, all samples" ],
    [ :concurrent_00to06, :concurrent_06to12, :concurrent_12to18, :concurrent_18to24, :future_allsamples ],
    [
        regression_currentCogScores_00to06_fromtaxa_results,
        regression_currentCogScores_06to12_fromtaxa_results,
        regression_currentCogScores_12to18_fromtaxa_results,
        regression_currentCogScores_18to24_fromtaxa_results
    ]
)

agemonths_importances = multimodel_aggregate_summaryimportances(ensemble_agemonths_models)
cogscore_importances = multimodel_aggregate_summaryimportances(ensemble_cogscore_models)

plot_df = leftjoin(agemonths_importances, cogscore_importances, on = :Variable, makeunique = true)

## Plot summary figure for Regression models
figure = Figure( ; resolution = (800, 600) )
ax = Axis(
    figure[1,1];
    xlabel = "AgeMonths Importance",
    ylabel = "CogScore Importance",
    title = "Comparative importance of taxa between the prediction ageMonths and cogScore"
)

# Plot barplot

scatter!(ax, plot_df.AvgMultimodelImportance, plot_df.AvgMultimodelImportance_1; color=:blue)

figure