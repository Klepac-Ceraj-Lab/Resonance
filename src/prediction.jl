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

## 3.1. Importance analysis
function calculate_fitness(train_cor::Vector{Float64}, test_cor::Vector{Float64})

    positive_train_cor = map( x -> maximum([x, 0.0]), train_cor)
    positive_test_cor = map( x -> maximum([x, 0.0]), test_cor)

    return positive_train_cor .* positive_test_cor

end

function weighted_hpimportances(m, hp = 1; change_hashnames=false, hashnamestable::Union{Nothing, DataFrame} = nothing, normalize_importances=true)
    merits = m.merits
    importances = m.importances
    fitnesses = calculate_fitness(merits.Train_Cor, merits.Test_Cor)
    fitness_subset_idx = findall(merits.Hyperpar_Idx .== hp)

    # Option 1: divide by amount of positive-only fitnesses
    # mean_importances = map(x -> sum(x .* fitnesses[fitness_subset_idx])/sum(fitnesses[fitness_subset_idx] .> 0.0), collect(eachrow(Matrix(importances[:, fitness_subset_idx .+ 1]))))
    # Option 2: divide importances by total number of models
    mean_importances = map(x -> sum(x .* fitnesses[fitness_subset_idx])/length(fitness_subset_idx), collect(eachrow(Matrix(importances[:, fitness_subset_idx .+ 1]))))
    mean_importances_df = sort( DataFrame(:variable => string.(importances[:,1]), :weightedImportance => mean_importances), :weightedImportance; rev=true)

    if normalize_importances
        mean_importances_df.weightedImportance = mean_importances_df.weightedImportance ./ sum(mean_importances_df.weightedImportance)
    end

    if change_hashnames
        newnames = leftjoin(mean_importances_df, hashnamestable, on = :variable => :hashname).longname
        mean_importances_df.variable = newnames
    end

    return mean_importances_df
end

function compute_joined_importances(modelwithoutdemo, modelwithdemo; imp_fun = weighted_hpimportances)

    importanceswithoutdemo = imp_fun(modelwithoutdemo)
    importanceswithdemo = imp_fun(modelwithdemo)
    rename!(importanceswithoutdemo, "weightedImportance" => "ImportanceWithoutDemo")
    rename!(importanceswithdemo, "weightedImportance" => "ImportanceWithDemo")

    joined_importances = dropmissing(outerjoin(importanceswithdemo, importanceswithoutdemo, on = [ :variable => :variable ]))
    sort!(joined_importances, :ImportanceWithoutDemo, rev = true)
    insertcols!(joined_importances, 1, :rank => collect(1:nrow(joined_importances)))

    return joined_importances

end

function plot_comparativedemo_importance_barplots!(plot_axis::Axis, joined_importances::DataFrame; n_rows = 40)

    # Plot barplot
    barplot!(plot_axis, reverse(collect(1:n_rows)), joined_importances.ImportanceWithoutDemo[1:n_rows], color = "gray", direction = :x)
    scatter!(plot_axis, joined_importances.ImportanceWithDemo[1:n_rows], reverse(collect(1:n_rows)), marker = :star8, color = "red", direction = :x)
    
    return plot_axis

end

## Comparing LM and RF
function attribute_colors(lm_qvalue, rf_cumImp, plot_colorset = [:dodgerblue, :orange, :green, :purple])
    if (lm_qvalue & rf_cumImp)
        return plot_colorset[4]
    elseif (!lm_qvalue & rf_cumImp)
        return plot_colorset[3]
    elseif (lm_qvalue & !rf_cumImp)
        return plot_colorset[2]
    else
        return plot_colorset[1]
    end
end

function plot_comparative_lmvsrf_scatterplots!(
    plot_axis::Axis,
    rf_model::ProbeData,
    lm_path::String;
    exclude_from_importances= [ "ageMonths" ],
    cumulative_importance_threshold = 0.6,
    plot_colorset = [:dodgerblue, :orange, :green, :purple])

    ## 1. calculate the importances and the cumulative sum, excluding the relevant variables
    rf_model_importances = weighted_hpimportances(rf_model; change_hashnames=false)
    rf_model_importances.relativeWeightedImportance = rf_model_importances.weightedImportance ./ sum(rf_model_importances.weightedImportance[.!(rf_model_importances.variable .∈ Ref(exclude_from_importances))])
    rf_model_importances[rf_model_importances.variable .∈ Ref(exclude_from_importances), :relativeWeightedImportance] .= 0.0
    rf_model_importances.cumulativeWeightedImportance = cumsum(rf_model_importances.relativeWeightedImportance)
    
    ## 2. get the data for linear models, from `Figure2-calculations.jl`
    lm_coefs = CSV.read(lm_path, DataFrame)
    lm_coefs = select(lm_coefs, [:feature, :coef, :pvalue, :qvalue])

    ## 3. Join the data
    plot_comparative_df = dropmissing(outerjoin(rf_model_importances, lm_coefs, on = [ :variable => :feature ]))
    
    ## 4. Attribute color to each point
    is_significative_lm = map(q-> q < 0.2 ? true : false, plot_comparative_df.qvalue)
    is_over_threshold = map(ci-> ci <= cumulative_importance_threshold ? true : false, plot_comparative_df.cumulativeWeightedImportance)
    point_colors = map((a, b) -> attribute_colors(a, b, plot_colorset), is_significative_lm, is_over_threshold)

    # @show DataFrame(:variable => plot_comparative_df.variable, :color => point_colors) # For Debug

    # 5. Plot scatterplot
    scatter!(
        plot_axis,
        log.(plot_comparative_df.pvalue).*(-1), plot_comparative_df.weightedImportance,
        color = point_colors
    )

    return plot_axis
end

## Taxon deepdives
function plot_taxon_deepdive!(figure_layout::GridLayout, figure_col::Int, spec::CommunityProfile, filter_row::Symbol, taxon_to_dive::String;)

    filtered_spec = spec[:, get(spec, filter_row)]
    (lq, uq) = quantile(get(filtered_spec, :cogScore), [0.25, 0.75])
    lqidx = findall(x-> x <= lq, get(filtered_spec, :cogScore))
    mqidx = findall(x-> lq < x <= uq, get(filtered_spec, :cogScore))
    uqidx = findall(x-> x > uq, get(filtered_spec, :cogScore))

    taxab = vec(abundances(filtered_spec[Regex(taxon_to_dive), :]))
    g1 = filter(!=(0), taxab[lqidx])
    g2 = filter(!=(0), taxab[mqidx])
    g3 = filter(!=(0), taxab[uqidx])

    g1p = count(!=(0), taxab[lqidx]) / length(lqidx)
    g2p = count(!=(0), taxab[mqidx]) / length(mqidx)
    g3p = count(!=(0), taxab[uqidx]) / length(uqidx)

    ax1 = Axis(
        figure_layout[1,figure_col];
        ylabel="relative abundance",#\n$taxon_to_dive",
        title = replace(taxon_to_dive, "_"=>" "),
        titlefont = "TeX Gyre Heros Makie Italic")
    ax2 = Axis(
        figure_layout[2,figure_col];
        xticks = (1:3, ["lower", "mid", "upper"]),
        xlabel="quartile",
        ylabel = "Intra-q\nPrevalence")
    ylims!(ax2, [0, 1.0])

    boxplot!(ax1, repeat([ 1.0 ], length(g1)), g1, color = (Makie.wong_colors()[1], 0.5), show_notch = false, show_outliers=false)
    boxplot!(ax1, repeat([ 2.0 ], length(g2)), g2, color = (Makie.wong_colors()[2], 0.5), show_notch = false, show_outliers=false)
    boxplot!(ax1, repeat([ 3.0 ], length(g3)), g3, color = (Makie.wong_colors()[3], 0.5), show_notch = false, show_outliers=false)    
    Random.seed!(0)
    scatter!(ax1, 1 .+ rand(Normal(0, 0.1), length(g1)), g1)
    Random.seed!(0)
    scatter!(ax1, 2 .+ rand(Normal(0, 0.1), length(g2)), g2)
    Random.seed!(0)
    scatter!(ax1, 3 .+ rand(Normal(0, 0.1), length(g3)), g3)
    barplot!(ax2, [1,2,3], [g1p, g2p, g3p]; color=Makie.wong_colors()[1:3])
    
    hidexdecorations!(ax1; grid=false)
    rowsize!(figure_layout, 2, Relative(1/3))

    # Code block to align the y axis labels on ax1 and ax2
    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace

    return figure_layout
end

## Importance supplementary table from ProbeData
function singlemodel_importances_suppltable(rf_model)
    rf_model_importances = weighted_hpimportances(rf_model; normalize_importances=false, change_hashnames=false)
    rf_model_importances.relativeWeightedImportance = rf_model_importances.weightedImportance ./ sum(skipmissing(rf_model_importances.weightedImportance))
    rf_model_importances.cumulativeWeightedImportance = cumsum(rf_model_importances.relativeWeightedImportance)
    return rf_model_importances
end

## Importance Pareto Plot from SuppTable
function plot_importances_pareto!(
    figpos::GridPosition,
    supptbl::DataFrame,
    figtitle = "Pareto Plot";
    barcolor = :lightblue,
    curvecolor = :orange,
    leftlims = [0.0, 10.0],
    rightlims = [0.0, 100.0]
    )
    
    axL = Axis(
        figpos,
        xlabel = "Feature",
        xticks = (collect(1:nrow(supptbl)), replace.(supptbl.variable, "_"=>" ")),
        xticklabelsize=8,
        xticklabelrotation= pi/4,
        xticklabelfont="TeX Gyre Heros Makie Italic",
        title = figtitle,
        ylabel = "Individual relative importance (%)")
    axR = Axis(
        figpos,
        xlabel = "Feature",
        xticks = (collect(1:nrow(supptbl))),
        ylabel = "Cumulative relative importance (%)",
        yticks = 0:10:100,
        yaxisposition = :right)

    ylims!(axL, leftlims)
    ylims!(axR, rightlims)
    hidexdecorations!(axR)
    linkxaxes!(axL, axR)

    barplot!(
        axL,
        collect(1:nrow(supptbl)), (100 .* supptbl.relativeWeightedImportance),
        color = barcolor, strokecolor = :black, strokewidth = 1)
    scatter!(
        axR,
        collect(1:nrow(supptbl)), (100 .* supptbl.cumulativeWeightedImportance),
        color = curvecolor, markersize = 15)
    lines!(
        axR,
        collect(1:nrow(supptbl)), (100 .* supptbl.cumulativeWeightedImportance),
        color = curvecolor, linewidth = 5)

    # tightlimits!(axL, Left(), Right())
    # tightlimits!(axR, Left(), Right()) # Not working as we want it to. Check on a later opportunity.

    return figpos
end