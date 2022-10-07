# Redistributable source code for training Univariate Random Forest models
# with the Resonance methodology
#
# Authors: Kevin Bonhan and Guilherme Bottino
# Version: 0.1
# run script with flag -h to see help and instructions
#####
# 0. Command-line Argument parsing function
#####

using Pkg; Pkg.activate(pwd()); Pkg.instantiate();

using ArgParse

function cli_args(args)

    s = ArgParseSettings("command-line arguments test")

    add_arg_group!(s, "Mutually exclusive prediction modes", exclusive=true, required=true)
    @add_arg_table! s begin
        "--classification", "-c"
            action = :store_true
            help = "Runs the Random Forest in Classifier mode"
        "--regression", "-r"
            action = :store_true
            help = "Runs the Random Forest in Regressor mode"
    end

    add_arg_group!(s, "Required arguments", required=true)
    @add_arg_table! s begin
        "input_csv"
            help = "Path to the CSV file containing the input and output features (columns) for each sample (rows), with optional additional metadata"
            required = true        # makes the argument mandatory
        "input_cols"
            help = "The column index range of the input columns (format 5:50)"
            required = true        # makes the argument mandatory
        "output_col"
            help = "The column index of the output column (only supports univariate output)"
            required = true        # makes the argument mandatory
        "output_filename"
            help = "Filename for the output JLD file"
            required = true        # makes the argument mandatory
    end

    s.epilog = """
    examples:\n
    \n
    \ua0\ua0julia $(basename(Base.source_path())) -c examples/example_taxonomic_profile.csv 6:554 4 examples/example_classifier_model.jld
    \n
    \ua0\ua0julia $(basename(Base.source_path())) -r examples/example_taxonomic_profile.csv 6:554 4 examples/example_regressor_model.jld
    \n
    The first example will train a classifier model using dhe data in examples/example_taxonomic_profile.csv.
    Columns 6 to 554 are the input features, while column 4 is the desired output.
    The result will be saved to models/test_output.jld.
    The second example does the same, except that it trains a regressor network.
    """

    parsed_args = parse_args(args, s; as_symbols = true)
    return parsed_args
end

parsed_args = cli_args(ARGS)

println("Parsed args:")
for (key,val) in parsed_args
    println("  $key  =>  $(repr(val))")
end

#####
# 1. Loading required packages
#####

using Chain
using CSV
using Statistics
using Random
using DataFrames
using StableRNGs
using MLJ
using DecisionTree
using JLD2
using Resonance
ml_rng = StableRNG(0)

#####
# 2. Setting up model
#####

input_df = @chain parsed_args[:input_csv] begin
    CSV.read(DataFrame)
    dropmissing()
end 

tuning_space = (
    maxnodes_range = collect(1:3:15) ,
    nodesize_range = collect(1:3:20),
    sampsize_range = [0.6, 0.7],
    mtry_range = collect(5:10:100),
    ntrees_range = [300, 500, 700]
)

if (parsed_args[:classification])

    RandomForestClassifier= MLJ.@load RandomForestClassifier pkg=DecisionTree

    rf_model = train_randomforest(
        Resonance.Classification(),
        parsed_args[:output_filename],
        input_df,
        identity,
        meanclass,
        eval(Meta.parse(parsed_args[:input_cols])),
        eval(Meta.parse(parsed_args[:output_col]));
        n_splits = 5,
        tuning_space = tuning_space,
        train_rng = ml_rng
    )
    JLD2.@save pwd()*"/"*parsed_args[:output_filename] rf_model
    
elseif (parsed_args[:regression])

    RandomForestRegressor= MLJ.@load RandomForestRegressor pkg=DecisionTree

    rf_model = train_randomforest(
        Resonance.Regression(),
        parsed_args[:output_filename],
        input_df,
        identity,
        eval(Meta.parse(parsed_args[:input_cols])),
        eval(Meta.parse(parsed_args[:output_col]));
        n_splits = 5,
        tuning_space = tuning_space,
        train_rng = ml_rng
    )
    JLD2.@save pwd()*"/"*parsed_args[:output_filename] rf_model

end