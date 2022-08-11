#####
# Notebook D02 - Prediction of futureCogScores from current taxonomic data
#####

using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using ProgressMeter
using MLJ
using CairoMakie
using DecisionTree
using CairoMakie
using GLM
ml_rng = StableRNG(0)

# ## 1. Collecting the data to make the future dataframe
include("D00-collect_taxonomic_cogscore_data.jl")

max_prediction_ageMonths = 24

### 2.2. Subsetting data
prediction_df = @chain cogscore_taxa_df begin
    subset(:ageMonths => x -> .!(ismissing.(x)))
    subset(:cogScore => x -> .!(ismissing.(x)))
    subset(:sample => x -> .!(ismissing.(x)))
    dropmissing()
    subset(:ageMonths => x -> x .<= max_prediction_ageMonths)
    # unique( [:subject] )
    # DataFrames.select(Not(:timepointDelta))
end

prediction_df_18_24 = @chain prediction_df begin
    subset(:ageMonths => x -> x .>= 18.0)
    subset(:ageMonths => x -> x .<= 24.0)    
end

X = prediction_df_18_24[:, 6:end]
y = prediction_df_18_24[:, 4]

for m in models(matching(X,y)) println(m) end

RandomForestRegressor= @load RandomForestRegressor pkg=DecisionTree

# Declaring the optimization grid
nodesize_range = collect(1:1:15)
maxnodes_range = collect(1:1:20)
sampsize_range = [0.5, 0.6, 0.7, 0.8]#collect(5:2:35) ./ nrow(X)
mtry_range = collect(5:5:100)
ntrees_range = [100, 300, 500, 700]    

tuning_grid = vec(collect(Base.product(
    nodesize_range,
    maxnodes_range,
    sampsize_range,
    mtry_range,
    ntrees_range
)))

#### Tuning loop

#### Splitting training data between train and test
Random.seed!(ml_rng, 1)
train, test = partition(eachindex(1:nrow(X)), 0.75, shuffle=true, rng=ml_rng)

train_mae_vector = Vector{Float64}(undef, length(tuning_grid))
test_mae_vector = Vector{Float64}(undef, length(tuning_grid))

bestmae = +Inf

@showprogress for i in 1:length(tuning_grid)

    global bestmae
    Random.seed!(ml_rng, 0)

    rf_model = RandomForestRegressor(
        n_trees = tuning_grid[i][5],
        min_samples_leaf = tuning_grid[i][1],
        max_depth = tuning_grid[i][2],
        sampling_fraction = tuning_grid[i][3],
        n_subfeatures = tuning_grid[i][4],
        rng=ml_rng
    )

    rf_machine = machine(rf_model, X[train, :], y[train])
    MLJ.fit!(rf_machine, verbosity=0)

    train_y_hat = MLJ.predict(rf_machine, X[train, :])
    train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))

    test_y_hat = MLJ.predict_mode(rf_machine, X[test, :]) 
    test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))

    train_mae_vector[i] = train_mae
    test_mae_vector[i] = test_mae

    if !isnan(test_mae)
        if(test_mae) < bestmae
            global bestmae = test_mae
            global bestmodel = deepcopy(rf_machine)
        end
    end

end

println(bestmae)

selected_machine = bestmodel

train_y_hat = MLJ.predict(selected_machine, X[train, :])
slope_correction_df = DataFrame([ train_y_hat, y[train] ], [:yhat, :y])
slope_correction = lm(@formula(y ~ yhat), slope_correction_df)
train_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => train_y_hat))
train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
train_cor = Statistics.cor(train_y_hat, y[train])

test_y_hat = GLM.predict(slope_correction, DataFrame(:yhat => MLJ.predict(selected_machine, X[test, :])))
test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
test_cor = Statistics.cor(test_y_hat, y[test])

train_set_plot = Figure(resolution = (800, 600))
ax = Axis(train_set_plot[1, 1],
    title = "Comparison of predicted vs ground truth *Current cogScore* for ages 18-24 months, *TRAIN* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[train], train_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(train_cor; digits = 2))"], [Point(60, 100)], textsize = 40)

train_set_plot

test_set_plot = Figure()
ax = Axis(test_set_plot[1, 1],
title = "Comparison of predicted vs ground truth *Current cogScore* for ages 18-24 months, *TEST* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[test], test_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(test_cor; digits = 2))"], [Point(60, 95)], textsize = 40)

test_set_plot

all_set_plot = Figure()
ax = Axis(all_set_plot[1, 1],
    title = "Comparison of predicted vs ground truth *Current cogScore* for ages 18-24 months,\n*FULL SAMPLE SET*",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)

yhat = GLM.predict(slope_correction, DataFrame( :yhat => MLJ.predict(selected_machine, X)))
sc_train = scatter!(ax, y[train], train_y_hat; color=:orange)
sc_test = scatter!(ax, y[test], test_y_hat; color=:blue)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(cor(y, yhat); digits = 2))"], [Point(60, 95)], textsize = 40)
Legend(all_set_plot[1,2], [sc_train, sc_test], ["Train", "Test"])

all_set_plot