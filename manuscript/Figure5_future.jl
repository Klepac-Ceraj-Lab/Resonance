#####
# This notebook works on Julia 1.8.5  with Resonance, but fails on 1.9.0-rc1 on Leap. Safekeeping to ask later.
#####
using Resonance
using FileIO
using CairoMakie # for plotting
using Distances
using MultivariateStats
using CategoricalArrays
using DataFrames
using CSV
using Chain
using StableRNGs
using MLJ
using CairoMakie
using DecisionTree
using JLD2
using ShapML
using Resonance
ml_rng = StableRNG(0)

## 0. Functions

## This function, when applied on _every_ 2-by-2 combinations of rows on a subdf created by `groupby`, returns whether there is _any_ viable combination, where viable is given by `check_viable_rowcombination`
function check_viable_subdf(df)
    viability = [ check_viable_rowcombination(df, row1, row2) for row1 in 1:nrow(df) for row2 in 1:nrow(df) ]
    return any(viability)
end

## This function checks whether a row combination is valid for future prediction purposes, in other words, if (I) T1 < t2, (II) there is an `omni` in T1 and (III) a `cogScore` in T2.
function check_viable_rowcombination(df, row1, row2)
    if (row1 >= row2)
        return false
    elseif (ismissing(df.omni[row1]))
        return false
    elseif (ismissing(df.cogScore[row2]))
        return false
    else 
        return true
    end
end

## This is a dirty and sleepy way to get all the viable i_1, i_2 combinations of rows for which `check_viable_rowcombination` is true.
function get_future_metadata(df)
    viabilities = @chain DataFrame([ (present = row1, future = row2, viable = check_viable_rowcombination(df, row1, row2)) for row1 in 1:nrow(df) for row2 in 1:nrow(df) ]) begin
        subset(:viable => identity)
    end
    future_df = DataFrame([(
                subject = df.subject[1],
                present_timepoint = df.timepoint[row.present],
                present_ageMonths = df.ageMonths[row.present],
                present_sample = df.omni[row.present],
                future_timepoint = df.timepoint[row.future],
                future_ageMonths = df.ageMonths[row.future],
                future_cogscore = df.cogScore[row.future]
            ) for row in eachrow(viabilities)])
    return(future_df)
end

## This is a function that the ShapML.jl package utilized to predict the shapley values
function predict_function(model, data)
    data_pred = DataFrame(:y_pred => vec(model(Matrix(data)')))
    return data_pred
end

#####

# Loading the raw ECHO metadata
echo_raw_metadata = Resonance.load_raw_metadata()

echo_taxonomic_profiles = @chain Resonance.load_raw_metaphlan() begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
    # select(Not([:file, :sample, :read_depth]))
end

rename_ref_table = CSV.read("oldname_newname_reltable.csv", DataFrame; stringtype = String)

echo_processed_future_data = @chain echo_raw_metadata begin
    select(["subject", "timepoint", "ageMonths", "cogScore", "hhs", "education", "race", "sex", "omni", "omni_collectionDate", "etoh", "etoh_collectionDate"])
    sort([:subject, :timepoint])
    # dropmissing(:omni)
    groupby(:subject)
        filter(:omni => (x -> length(x) .> 1), _ )
        filter(:omni => (x -> any(.!(ismissing.(x)))), _ )
        filter(check_viable_subdf, _ )
    DataFrames.combine(get_future_metadata)
    dropmissing( [ :present_ageMonths, :future_ageMonths ] ) # filters out maternal samples
    innerjoin(rename_ref_table, on = :present_sample => :old_tubeid)
    innerjoin(echo_taxonomic_profiles, on = :new_seqid => :sample_base)
    select!(Not([:present_sample, :file ]))
    subset(:present_ageMonths => x -> x .<= 18.0)
    subset(:future_ageMonths => x -> 12.0 .<= x)
    sort([:subject, :present_timepoint, :future_timepoint]; rev=true)
    unique(:subject)
    dropmissing()
end

echo_processed_future_data = echo_processed_future_data[ findall([ sum(j) for j in eachrow(echo_processed_future_data[:,9:end]) ] .> 0.0), vcat(1:8..., findall([ any( k .> 0.0) for k in eachcol(echo_processed_future_data[:,9:end]) ]) .+ 8) ]
echo_processed_future_data.future_cogscore = parse.(Float64, echo_processed_future_data.future_cogscore)

CSV.write("ECHO_FutureCogscore_inputs.csv", echo_processed_future_data)

#####
# Figure with age of collection and age of assessment
#####

fig = Figure(; resolution = (400,500))
ax = Axis(
    fig[1, 1];
    xlabel = "Timepoint event",
    xticks = ( [ 1, 2 ] , [ "stool\ncollection", "coognitive\nassessment" ] ),
    ylabel = "Age (months)",
    title = "Age distributions for model",
)

for rr in eachrow(echo_processed_future_data)
    scatterlines!(ax, [1,2] , [ rr.present_ageMonths, rr.future_ageMonths ]; color = (:black, 0.3))
end

# boxplot!(ax, [ 1 for _ in eachrow(echo_processed_future_data) ], echo_processed_future_data.present_ageMonths; color = (:blue, 0.5), show_outliers = false )
# boxplot!(ax, [ 2 for _ in eachrow(echo_processed_future_data) ], echo_processed_future_data.future_ageMonths; color = (:red, 0.5), show_outliers = false )

violin!(ax, [ 1 for _ in eachrow(echo_processed_future_data) ], echo_processed_future_data.present_ageMonths; color = (:blue, 0.5), datalimits = extrema )
violin!(ax, [ 2 for _ in eachrow(echo_processed_future_data) ], echo_processed_future_data.future_ageMonths; color = (:red, 0.5), datalimits = extrema  )

save("manuscript/assets/presentstool_futurecogscore_agepathdist.png")

#####
# Applying prevalence filters
#####
taxa_prevalence_threshold_lower = 0.10
taxa_prevalence_threshold_upper = 1.00

## 00 to 06 months

filtered_echo_processed_future_data = filter_prevalences(
    echo_processed_future_data,
    :future_cogscore,
    [:subject, :present_timepoint, :present_ageMonths, :future_timepoint, :future_ageMonths, :sample],
    [:new_seqid];
    education2int = false,
    sex2int = false,
    lbound = taxa_prevalence_threshold_lower,
    ubound = taxa_prevalence_threshold_upper
)

filtered_echo_processed_future_data[!,:present_ageMonths], filtered_echo_processed_future_data[!,:sample] = filtered_echo_processed_future_data[!,:sample], filtered_echo_processed_future_data[!,:present_ageMonths]
rename!(filtered_echo_processed_future_data, [ :present_ageMonths => :sample, :sample => :present_ageMonths ])
#####
# Actual model training
#####

taxa_tuning_space = (
    maxnodes_range = [ -1 ],
    nodesize_range = [ 5 ],
    sampsize_range =  [ 0.7 ],
    mtry_range = [ -1 ],
    ntrees_range = [ 500 ]
)

#####
# 00 to 06 months
#####

## 3. SES + taxonomic profiles
regression_futureCogScores_demoplustaxa = probe_prod_randomforest(
    Resonance.Regression(),
    "regression_futureCogScores_demoplustaxa",
    filtered_echo_processed_future_data,
    identity,
    collect(8:ncol(filtered_echo_processed_future_data)),
    :future_cogscore;
    n_folds = 3,
    n_replicas = 100,
    n_rngs = 10,
    tuning_space = taxa_tuning_space,
    train_rng=ml_rng
)

conf_interval(vv, critical_z=1.96) = critical_z * std(vv) / sqrt(length(vv))
report_regression_merits(m) = DataFrames.combine(
    groupby(m.merits, :Hyperpar_Idx),
    :Train_RMSE => (x -> round(mean(x); digits = 4)) => :Train_RMSE_mean,
    :Train_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Train_RMSE_CI,
    :Train_Cor => (x -> round(mean(x); digits = 4)) => :Train_Cor_mean,
    :Train_Cor => (x -> round(conf_interval(x); digits = 4)) => :Train_Cor_CI,
    :Test_RMSE => (x -> round(mean(x); digits = 4)) => :Test_RMSE_mean,
    :Test_RMSE => (x -> round(conf_interval(x); digits = 4)) => :Test_RMSE_CI,
    :Test_Cor => (x -> round(mean(x); digits = 4)) => :Test_Cor_mean,
    :Test_Cor => (x -> round(conf_interval(x); digits = 4)) => :Test_Cor_CI,
    :Hyperpar_Idx=> (x-> m.name) => :model
    )

report_regression_merits(regression_futureCogScores_demoplustaxa)

weighted_hpimportances(regression_futureCogScores_demoplustaxa)