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
ml_rng = StableRNG(0)

#####
# Functions
#####

function get_quantile_class(cogscorevector)
    quantiles = quantile(cogscorevector, [ 0.5, 1.0 ])
    quantile_class = map(x -> findfirst(x .<= quantiles), cogscorevector) .- 1
    return Bool.(quantile_class)
end

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

prediction_df_0_6 = @chain prediction_df begin
    subset(:ageMonths => x -> x .>= 0.0)
    subset(:ageMonths => x -> x .<= 6.0)    
end

hist(prediction_df_0_6.cogScore, bins = 6, color = :blue, strokewidth = 2, strokecolor = :black)

prediction_df_6_12 = @chain prediction_df begin
    subset(:ageMonths => x -> x .>= 6.0)
    subset(:ageMonths => x -> x .<= 12.0)    
end

hist(prediction_df_6_12.cogScore, bins = 8, color = :blue, strokewidth = 2, strokecolor = :black)

prediction_df_12_18 = @chain prediction_df begin
    subset(:ageMonths => x -> x .>= 12.0)
    subset(:ageMonths => x -> x .<= 18.0)    
end

hist(prediction_df_12_18.cogScore, bins = 7, color = :blue, strokewidth = 2, strokecolor = :black)

prediction_df_18_24 = @chain prediction_df begin
    subset(:ageMonths => x -> x .>= 18.0)
    subset(:ageMonths => x -> x .<= 24.0)    
end

hist(prediction_df_18_24.cogScore, bins = 8, color = :blue, strokewidth = 2, strokecolor = :black)

X = prediction_df_6_12[:, 6:end]
y = @chain prediction_df_6_12.cogScore begin
    get_quantile_class()
    coerce(OrderedFactor)
end
for m in models(matching(X,y)) println(m) end

RandomForestClassifier= @load RandomForestClassifier pkg=DecisionTree

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

train_acc_vector = Vector{Float64}(undef, length(tuning_grid))
test_acc_vector = Vector{Float64}(undef, length(tuning_grid))

bestacc = 0.0

@showprogress for i in 1:length(tuning_grid)

    global bestacc
    Random.seed!(ml_rng, 0)

    rf_model = RandomForestClassifier(
        n_trees = tuning_grid[i][5],
        min_samples_leaf = tuning_grid[i][1],
        max_depth = tuning_grid[i][2],
        sampling_fraction = tuning_grid[i][3],
        n_subfeatures = tuning_grid[i][4],
        rng=ml_rng
    )

    rf_machine = machine(rf_model, X[train, :], y[train])
    MLJ.fit!(rf_machine, verbosity=0)

    train_y_hat = MLJ.predict_mode(rf_machine, X[train, :])

    train_acc = mean(train_y_hat .== y[train])

    test_y_hat = MLJ.predict_mode(rf_machine, X[test, :]) 
    test_acc = mean(test_y_hat .== y[test])

    train_acc_vector[i] = train_acc
    test_acc_vector[i] = test_acc

    if !isnan(test_acc)
        if(test_acc) > bestacc
            global bestacc = test_acc
            global bestmodel = deepcopy(rf_machine)
        end
    end

end

println(bestacc)