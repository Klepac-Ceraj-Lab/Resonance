#####
# Notebook DN01 - Regression of binary current CogScores above/below 50th percentile -- NULL DISTRIBUTION
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

## 1. Taxonomic Profiles
taxa = Resonance.load(TaxonomicProfiles())
mdata_taxa_df = sort(Resonance.comm2wide(taxa), [ :subject, :timepoint ]);

## 2. Brain data
brain = Resonance.load(Neuroimaging())
mdata_brain_df = sort(Resonance.comm2wide(brain; feature_func = string), [ :subject, :timepoint ]);
mdata_brain_df = dropmissing(mdata_brain_df[mdata_brain_df.filter_18to120, :])
TBV = mdata_brain_df."White-matter" .+ mdata_brain_df."Gray-matter"

for col in names(mdata_brain_df)[16:end-2]
    mdata_brain_df[!, col] = 100 * mdata_brain_df[:, col] ./ TBV
end

## Filter taxonomic profile rows to match brain samples from 18 to 20 months

taxa_row_indexes = [ findfirst(mdata_taxa_df.sample .== el) for el in mdata_brain_df.sample ]
mdata_taxa_df = mdata_taxa_df[taxa_row_indexes, :]

#####
# Applying prevalence filters
#####
taxa_prevalence_threshold_lower = 0.10
taxa_prevalence_threshold_upper = 1.00

filtered_mdata_taxa_df_18to120 = filter_prevalences(
    dropmissing(mdata_taxa_df[mdata_taxa_df.filter_18to120, :], Not([:race, :read_depth])),
    :cogScore,
    [:sex, :education, :ageMonths],
    [:sample, :sample_base, :subject, :timepoint, :race, :date, :read_depth, :filter_00to120, :filter_00to06, :filter_18to120];
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)[:, 2:end]

CSV.write("18to120_brain_taxa.csv", hcat(filtered_mdata_taxa_df_18to120, mdata_brain_df[:, 16:end-2]))

#####
# Training models
#####

RandomForestRegressor = MLJ.@load RandomForestRegressor pkg=DecisionTree

taxa_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

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

# #####
# # 18 to 120 months
# #####

# Summary of brain segment aggregations:
# Removed "left-lateral-ventricle", "right-lateral-ventricle", "left-inferior-lateral-ventricle", "right-inferior-lateral-ventricle", because vessels are and not necessarily interesting;
# Removed "left-cerebellum-exterior", "right-cerebellum-exterior"; not useful and bad predictions.
# Removed "left-ventral-DC", "right-ventral-DC". See http://neuromorphometrics.com/Seg/html/segmentation/ventral%20diencephalon.html
# Combined "superior-temporal", "middle-temporal", "transverse-temporal" and "inferior-temporal" into "temporal"
mdata_brain_df."left-temporal" = mdata_brain_df."left-superior-temporal" .+ mdata_brain_df."left-middle-temporal" .+ mdata_brain_df."left-transverse-temporal" .+ mdata_brain_df."left-inferior-temporal"
mdata_brain_df."right-temporal" = mdata_brain_df."right-superior-temporal" .+ mdata_brain_df."right-middle-temporal" .+ mdata_brain_df."right-transverse-temporal" .+ mdata_brain_df."right-inferior-temporal"
# Combined "lateral-orbitofrontal" and "medial-orbitofrontal" into "orbitofrontal"
mdata_brain_df."left-orbitofrontal" = mdata_brain_df."left-lateral-orbitofrontal" .+ mdata_brain_df."left-medial-orbitofrontal"
mdata_brain_df."right-orbitofrontal" = mdata_brain_df."right-lateral-orbitofrontal" .+ mdata_brain_df."right-medial-orbitofrontal"
# Combined "inferior-parietal" and "superior-parietal" into "parietal"
mdata_brain_df."left-parietal" = mdata_brain_df."left-inferior-parietal" .+ mdata_brain_df."left-superior-parietal"
mdata_brain_df."right-parietal" = mdata_brain_df."right-inferior-parietal" .+ mdata_brain_df."right-superior-parietal"
# Combined "caudal-middle-frontal" and "rostral-middle-frontal" into "middle-frontal"
mdata_brain_df."left-middle-frontal" = mdata_brain_df."left-caudal-middle-frontal" .+ mdata_brain_df."left-rostral-middle-frontal"
mdata_brain_df."right-middle-frontal" = mdata_brain_df."right-caudal-middle-frontal" .+ mdata_brain_df."right-rostral-middle-frontal"
# Combined "caudal-anterior-cingulate" and "ostral-anterior-cingulate" into "anterior-cingulate"
mdata_brain_df."left-anterior-cingulate" = mdata_brain_df."left-caudal-anterior-cingulate" .+ mdata_brain_df."left-rostral-anterior-cingulate"
mdata_brain_df."right-anterior-cingulate" = mdata_brain_df."right-caudal-anterior-cingulate" .+ mdata_brain_df."right-rostral-anterior-cingulate"

brain_segments_list = [
    "left-temporal", "right-temporal",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-parietal", "right-parietal",
    "left-middle-frontal", "right-middle-frontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-cerebellum-white-matter", "right-cerebellum-white-matter",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-caudate", "right-caudate",
    "left-putamen", "right-putamen",
    "left-pallidum", "right-pallidum",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-basal-forebrain", "right-basal-forebrain", 
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lingual", "right-lingual",
    "left-parahippocampal", "right-parahippocampal",
    "left-paracentral", "right-paracentral",
    "left-pars-opercularis", "right-pars-opercularis",
    "left-pars-orbitalis", "right-pars-orbitalis",
    "left-pars-triangularis", "right-pars-triangularis",
    "left-pericalcarine", "right-pericalcarine",
    "left-postcentral", "right-postcentral",
    "left-posterior-cingulate", "right-posterior-cingulate",
    "left-precentral", "right-precentral",
    "left-precuneus", "right-precuneus",
    "left-superior-frontal", "right-superior-frontal",
    "left-supramarginal", "right-supramarginal",
    "left-insula", "right-insula",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
    "Brain-stem", "CSF"
]

# ## 7. Only taxonomic profiles

brain_models = Dict{String, ProbeData}()

for brain_segment in brain_segments_list

    tmp_df = deepcopy(filtered_mdata_taxa_df_18to120)
    prediction_df = insertcols!(tmp_df, 1, :target => mdata_brain_df[:, brain_segment])
    model_name = "regression_"*string(brain_segment)*"_18to120mo_onlytaxa"

    push!(
        brain_models,
        brain_segment => probe_prod_randomforest(
            Resonance.Regression(),
            model_name,
            prediction_df,
            identity,
            collect(2:ncol(prediction_df)),
            :target;
            n_folds = 3,
            n_replicas = 100,
            n_rngs = 10,
            tuning_space = taxa_tuning_space,
            train_rng=ml_rng
        )
    )

    println("finished!")

end
    
JLD2.@save "models/brain_models.jld" brain_models

@show [ mean(brain_models[el].merits.Test_Cor) for el in brain_segments_list ]

mean_brain_merits = sort(
    reduce(
        vcat,
        [
            DataFrame(
                :variable => i,
                :Train_RMSE => mean(j.merits.Train_RMSE),
                :Test_RMSE => mean(j.merits.Test_RMSE),
                :Train_MAPE => mean(j.merits.Train_MAPE),
                :Test_MAPE => mean(j.merits.Test_MAPE),
                :Train_Cor => mean(j.merits.Train_Cor),
                :Test_Cor => mean(j.merits.Test_Cor)
            ) for (i,j) in brain_models
        ]
    ), :Test_Cor, rev=true
)

CSV.write("brain_segmented_merits_SuppTable.csv", mean_brain_merits)

combine(mean_brain_merits,
    [ :Test_Cor => mean => :Test_Cor_mean,
    :Test_Cor => std => :Test_Cor_std,
    :Test_Cor => maximum => :Test_Cor_maximum,
    :Test_Cor => (x -> quantile(x, 0.75)) => :Test_Cor_upperquartile,
    :Test_Cor => (x -> quantile(x, 0.50)) => :Test_Cor_median,
    :Test_Cor => (x -> quantile(x, 0.25)) => :Test_Cor_lowerquartile,
    :Test_Cor => minimum => :Test_Cor_minimum ]
    )


combine(mean_brain_merits,
    [
        :Test_MAPE => mean => :Test_MAPE_mean,
        :Test_MAPE => std => :Test_MAPE_std,
        :Test_MAPE => maximum => :Test_MAPE_maximum,
        :Test_MAPE => (x -> quantile(x, 0.75)) => :Test_MAPE_upperquartile,
        :Test_MAPE => (x -> quantile(x, 0.50)) => :Test_MAPE_median,
        :Test_MAPE => (x -> quantile(x, 0.25)) => :Test_MAPE_lowerquartile,
        :Test_MAPE => minimum => :Test_MAPE_minimum ]
)