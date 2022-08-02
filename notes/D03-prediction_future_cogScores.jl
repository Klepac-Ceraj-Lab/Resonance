#####
# Functions to use
#####

function tryparse(Type::DataType, col)
    try
        return(parse.(Type, col))
    catch
        return(col)
    end
end

function univariate_tietjenmoore(values::Vector{T} where T <: Real, k::Int64; alpha = 0.05, sim_trials = 250000, return_indexes = true)

    function compute_tietjenmoore(data, k, n)
        r_all = abs.(data .- mean(data))
        filteredData = data[sortperm(r_all)][1:(n-k)]
        ksub = (filteredData .- mean(filteredData))
    
        return( sum(ksub .^ 2) / sum(r_all .^ 2) )
    end
    
    function test_tietjenmoore(dataSeries,k, n, alpha, sim_trials)
        ek = compute_tietjenmoore(dataSeries,k, n)
        simulation = [ compute_tietjenmoore(randn(length(dataSeries)), k, n) for i in 1:sim_trials ]
        Talpha=round(quantile(Normal(Statistics.mean(simulation), Statistics.std(simulation)),alpha), digits = 3)

        return(ek, Talpha)
    end

    println("----- Begin Tietjen-Moore Outlier test -----")
    println("H0: There are no outliers in the data set.")
    println("H1: There are exactly k outliers in the data set\n")

    n = length(values)

    L_set, L_critical = test_tietjenmoore(values, k, n, alpha, sim_trials)
    println("Set L-statistic for $n samples and $k outliers with mode $(mode): $(round(L_set, digits = 4))")     
    println("Critical L for $n samples and $k outliers with mode $(mode): $(round(L_critical, digits = 4))\n")

    if L_set < L_critical
        println("L_set < L_critical !")
        println("**SUCCESSFUL REJECTION OF H0** with confidence level $alpha" )

        r_all = abs.(values .- mean(values))
        outlier_indexes = sortperm(r_all)[(n-k+1):end]

        return_indexes ? (return(outlier_indexes)) : (return(true))

    else
        println("L_set > L_critical !")
        println("**CANNOT REJECT H0** with confidence level $alpha" )

        return_indexes ? (return([])) : (return(false))    
    
    end # endif L_set < L_critical
end # end function

function try_outliers(f, data, n)
    for i in n:-1:1
        outlier_idx = f(data, i)
        if length(outlier_idx) != 0
            return(i, outlier_idx)
        end
    end
    return(0, Int64[])
end

#####
# Script
#####

# ## Parsing tidy CSV data

rawDf = CSV.read(
    joinpath("data", "tidy_timepoints_with_brain.csv"),
    DataFrame;
    delim = ',',
    types = [Int64, Int64, String, String],
    header = ["subject", "timepoint", "variable", "value"],
    skipto = 2
)

# ### Fixing some metadata names (see source for details)

#fix_metadata_colnames!(rawDf)

# ## Removing duplicates

# ### Removing actual duplicated lines

check_longdata_metaduplicates!(rawDf; remove_duplicates=true)

# ### Removing technical/biological replicates

retainfirst = get(ENV, "RETAIN_FIRST", true) # if true, retain the first technical/biological replicate; if false, retain the last technical/biological replicate

if retainfirst

    unstackedDf = @chain rawDf begin
        groupby( [:subject, :timepoint, :variable] )
        combine(:value=>first)
        unstack(:variable, :value_first)
    end

else 
    unstackedDf = unstack(rawDf, :variable, :value; allowduplicates=true)
end

## Prediction workflows

### Problem 2: predict future cogScore from current taxonomic profile

#### Building predictionDf

select_columns = [ "target", "ageMonths", "futureAgeMonths", "ageMonthsDelta", "timepointDelta", microbiome_predictors... ]

##### 1. Removing rows with missing values and applying filters

max_input_age = 12
max_prediction_age = 24
max_ageMonthsDelta = 12
selected_timepoint_delta = 1

predictionDf = @chain unstackedDf begin
    subset(:ageMonths => x -> .!(ismissing.(x)))
    build_future_df(:cogScore)
    DataFrames.select(select_columns)
end

mapcols!(col -> tryparse(Float64, col), predictionDf)

##### 2. Filtering by ageMonths < max_prediction_age

predictionDf = @chain predictionDf begin
    subset(:ageMonths => x ->  x .<= max_input_age)
    subset(:futureAgeMonths => x ->  x .<= max_prediction_age)
    subset(:ageMonthsDelta => x ->  x .<= max_ageMonthsDelta)
    subset(:timepointDelta => x ->  x .== selected_timepoint_delta)
    DataFrames.select(Not(:timepointDelta))
#    DataFrames.transform(:ageMonths => ByRow(x -> x / max_prediction_age) => [:ageMonths])
end

##### Removing outliers with the Tietjen-Moore test
n_outliers, outlier_idx = try_outliers(univariate_tietjenmoore, predictionDf.target, 10)
println("Found $n_outliers outliers with the Tietjen-Moore test! Filtering DataFrame")
println("New prediction table was reduced from $(nrow(predictionDf)) to $(nrow(predictionDf) - n_outliers) rows.")
predictionDf = predictionDf[Not(outlier_idx), :]

##### 3. Removing microbiome columns with zero standard deviation

std_cutoff = 0.1
standard_deviations = [ std(predictionDf[:,i]) for i in 1:ncol(predictionDf) ]
is_zerostd_col = standard_deviations .== 0
is_negligiblestd_col = standard_deviations .<= std_cutoff
println("$(sum(is_zerostd_col)) columns have zero standard deviation!")
println("$(sum(is_negligiblestd_col)) columns have a standard deviation lower than $std_cutoff")
zerostd_cols = findall(is_zerostd_col)
negligiblestd_cols = findall(is_negligiblestd_col)
zerostd_colnames = names(predictionDf)[zerostd_cols]
negligiblestd_colnames = names(predictionDf)[negligiblestd_cols]

# Uncomment to remove zerostd columns
predictionDf = DataFrames.select(predictionDf, Not(zerostd_colnames))

#### Checking the schema and MLJ scitypes

schema(predictionDf)

#### Separating the predictors and target variable, vewing matching machine models

ml_rng = StableRNG(102528)

y, X = unpack(predictionDf, ==(:target));
for m in models(matching(X,y)) println(m) end

#### Computing the standardizer parameters

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree

#### Building the Tuning grid for training the Random Forest model
# Hyperparameters of DecisionTree.RandomForestRegressor:
# ```
# n_trees = 10,                 Number of trees to train. Equivalent to R's randomForest `ntree`
# n_subfeatures = 0             Number of features to consider at random per split. Equivalent to R's randomForest `mtry`
# partial_sampling = 0.7        Fraction of samples to train each tree on. Equivalent to R's randomForest `sampsize`, but R takes an integer and this is a fraction.
# max_depth = -1,               Maximum depth of the decision trees. Equivalent to R's randomForest `maxnodes`
# min_samples_leaf = 1,         Equivalent to R's randomForest `nodesize`
# min_samples_split = 2,        No equivalence.
# min_purity_increase = 0.0,    No equivalence.
# ````

range_n_trees = collect(50:25:500)
range_max_depth = collect(-1:3)
range_partial_sampling = collect(0.5:0.1:0.9)
range_min_samples_leaf = collect(1:1:5)
range_min_purity_increase = [ 0.0, 0.1, 0.2, 0.3 ]

# range_n_trees = collect(10:30:310)#collect(10:10:150)
# range_n_subfeatures = collect(10:20:150)#collect(10:20:250)
# range_partial_sampling = collect(0.3:0.1:0.7)
# range_max_depth = [ -1 ]#collect(4:2:20)
# range_min_samples_leaf = collect(1:1:5)#collect(1:1:10)
# range_min_samples_split = [ 2, 3, 4 ]#[ 2, 3, 4, 5, 6, 7, 8 ]s
# range_min_purity_increase = [ 0.0, 0.05, 0.10, 0.15 ]#[ 0.00, 0.05, 0.10, 0.15, 0.20, 0.25 ]

tuning_grid = vec(collect(Base.product(
    range_n_trees,
    range_max_depth,
    range_partial_sampling,
    range_min_samples_leaf,
    range_min_purity_increase,
)))

#### Tuning loop

#### Splitting training data between train and test
train, test = partition(eachindex(1:nrow(predictionDf)), 0.6, shuffle=true, rng=ml_rng)

train_mae_vector = Vector{Float64}(undef, length(tuning_grid))
train_cor_vector = Vector{Float64}(undef, length(tuning_grid))
test_mae_vector = Vector{Float64}(undef, length(tuning_grid))
test_cor_vector = Vector{Float64}(undef, length(tuning_grid))

bestcor = 0.0

@showprogress for i in 1:length(tuning_grid)

    global bestcor

    rf_model = RandomForestRegressor(
        n_trees = tuning_grid[i][1],
        sampling_fraction = tuning_grid[i][3],
        max_depth = tuning_grid[i][2],
        min_samples_leaf = tuning_grid[i][4],
        min_purity_increase = tuning_grid[i][5],
        rng=ml_rng
        # n_trees = tuning_grid[i][1],
        # n_subfeatures = tuning_grid[i][2],
        # sampling_fraction = tuning_grid[i][3],
        # max_depth = tuning_grid[i][4],
        # min_samples_leaf = tuning_grid[i][5],
        # min_samples_split = tuning_grid[i][6],
        # min_purity_increase = tuning_grid[i][7],
        # rng=ml_rng
    )

    rf_machine = machine(rf_model, X[train, :], y[train])
    fit!(rf_machine, verbosity=0)

    train_y_hat = predict(rf_machine, X[train, :]) 
    train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
    train_cor = Statistics.cor(train_y_hat, y[train])

    test_y_hat = predict(rf_machine, X[test, :])
    test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
    test_cor = Statistics.cor(test_y_hat, y[test])

#    machines_vector[i] = deepcopy(rf_machine)
    train_mae_vector[i] = train_mae
    train_cor_vector[i] = train_cor
    test_mae_vector[i] = test_mae
    test_cor_vector[i] = test_cor

    if !isnan(test_cor)
        if(test_cor) > bestcor
            global bestcor = test_cor
            global bestmodel = deepcopy(rf_machine)
        end
    end

end

println(minimum(train_mae_vector))
println(minimum(test_mae_vector))
println(maximum(train_cor_vector))
println(maximum(test_cor_vector))

#### Selecting model with the best performance on the test set

selected_machine = bestmodel

train_y_hat = predict(selected_machine, X[train, :]) 
train_mae = mean(MLJ.MLJBase.l1(train_y_hat, y[train]))
train_cor = Statistics.cor(train_y_hat, y[train])

test_y_hat = predict(selected_machine, X[test, :])
test_mae = mean(MLJ.MLJBase.l1(test_y_hat, y[test]))
test_cor = Statistics.cor(test_y_hat, y[test])

#### Plotting and saving figures

train_set_plot = Figure(resolution = (800, 600))
ax = Axis(train_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Future cogScore for the *TRAIN* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[train], train_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(train_cor; digits = 2))"], [Point(60, 100)], textsize = 40)

save("figures/prediction_futurecogscore_trainset_comparison.png", train_set_plot)

test_set_plot = Figure()
ax = Axis(test_set_plot[1, 1],
    title = "Comparison of predicted and ground truth Future cogScore for the *TEST* set",
    xlabel = "Ground Truth",
    ylabel = "Prediction"
)
scatter!(ax, y[test], test_y_hat)
ablines!(ax, 1, 1)
annotations!( ax, ["r = $(round(test_cor; digits = 2))"], [Point(60, 95)], textsize = 40)

save("figures/prediction_futurecogscore_testset_comparison.png", test_set_plot)

using DecisionTree

@assert cor(
    apply_forest(selected_machine.fitresult, Matrix(X[test, :])),
    test_y_hat
) == 1.0

impurity_importances = impurity_importance(selected_machine.fitresult)
impurity_ordered_idx = sortperm(impurity_importances; rev = true)

impurityDf = DataFrame(
    :Variable => names(X)[impurity_ordered_idx],
    :Importance => impurity_importances[impurity_ordered_idx]
)

### Using linear regression to adjust future cogScore ~ futureAgeMonths and train on residuals

predictionDf[sortperm(predictionDf.target), :]
predictionDf[predictionDf.ageMonthsDelta .>= 12.0, :]

futurecogscore_futureagemonths_exploration = Figure()
ax = Axis(futurecogscore_futureagemonths_regression[1, 1],
    title = "futureCogScore as a function of futureAgeMonths",
    xlabel = "futureAgeMonths",
    ylabel = "futureCogScore"
)
scatter!(ax, predictionDf.futureAgeMonths, predictionDf.target)

futurecogscore_futureagemonths_exploration

using GLM

fm = @formula(target ~ futureAgeMonths)
linear_regressor = lm(fm, predictionDf)

intercept = linear_regressor.model.pp.beta0[1]
slope = linear_regressor.model.pp.beta0[2]

linear_prediction = GLM.predict(linearRegressor, hcat(ones(Float64, nrow(predictionDf)), predictionDf.futureAgeMonths))
residuals = predictionDf.target .- linear_prediction

futurecogscore_futureagemonths_regression = Figure()
ax = Axis(futurecogscore_futureagemonths_regression[1, 1],
    title = "true futureCogScore as a function of regression output from futureAgeMonths",
    xlabel = "futureCogScore",
    ylabel = "estimation of futureCogScore"
)
scatter!(ax, predictionDf.target, linear_prediction)

futurecogscore_futureagemonths_regression

####

predictionDf.target .= residuals