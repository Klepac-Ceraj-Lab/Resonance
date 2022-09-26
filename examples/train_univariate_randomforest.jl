# Redistributable source code for training Univariate Random Forest models with the Resonance methodology
# Authors: Kevin Bonhan and Guilherme Bottino
# Version: 0.1
#####
# 0. Command-line Argument parsing function
#####

using Pkg; Pkg.activate(replace( ENV["DIRENV_DIR"], '-' => "")); Pkg.instantiate();

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
            help = "Absolute path (ex: `/home/user/documents/data.csv` to the CSV file containing the input and output features (columns) for each sample (rows)"
            required = true        # makes the argument mandatory
        "output_index"
            help = "The column index of the output column (only supports univariate output)"
            required = true        # makes the argument mandatory
        "input_index"
            help = "The column index range of the input columns (format 5:50)"
            required = true        # makes the argument mandatory
        "output_filename"
            help = "Filename for the output JLD file"
            required = true        # makes the argument mandatory
    end

    parsed_args = parse_args(args, s)
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

input_df = CSV.read(parsed_args[:input_csv])