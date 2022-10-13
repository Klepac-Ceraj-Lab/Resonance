# Redistributable source code for predicting outputs using
# pretrained Univariate Random Forest models with the Resonance methodology
#
# Authors: Kevin Bonhan and Guilherme Bottino
# Version: 0.1
# run script with flag -h to see help and instructions
# edit constant `tuning_space` to change hyperparameter search grid
#####
# 0. Command-line Argument parsing function
#####

using Pkg; Pkg.activate(pwd()); Pkg.instantiate();

using ArgParse

function cli_args(args)

    s = ArgParseSettings("command-line arguments test")

    add_arg_group!(s, "Required arguments", required=true)
    @add_arg_table! s begin
        "model_jldpath"
            help = "Path to the JLD file containing a pretrained model from `train_univariate_randomforest.jl`"
            required = true        # makes the argument mandatory
        "input_csv"
            help = "Path to the CSV file containing the input and output features (columns) for each sample (rows), with optional additional metadata"
            required = true        # makes the argument mandatory
        "input_cols"
            help = "The column index range of the input columns (format 5:50)"
            required = true        # makes the argument mandatory
        "output_filename"
            help = "Filename for the output predictions, CSV format"
            required = true        # makes the argument mandatory
    end

    ##[TODO]: add option to output the names of the expected columns for a single file

    s.epilog = """
    examples:\n
    \n
    \ua0\ua0julia examples/$(basename(Base.source_path())) examples/example_classifier_model.jld examples/example_taxonomic_profile.csv 6:554 examples/classifier_output.csv
    \n
    \ua0\ua0julia examples/$(basename(Base.source_path())) examples/example_regressor_model.jld examples/example_taxonomic_profile.csv 6:554 examples/regressor_output.csv
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
# 2. Loading model and input data
#####

## Model
RandomForestClassifier= MLJ.@load RandomForestClassifier pkg=DecisionTree verbosity=0
RandomForestRegressor= MLJ.@load RandomForestRegressor pkg=DecisionTree verbosity=0

model = JLD2.load(parsed_args[:model_jldpath])
model = model[first(keys(model))]

## Data
input_df = @chain parsed_args[:input_csv] begin
    CSV.read(DataFrame)
    dropmissing()
end 

#####
# 3. Generating predictions
#####

function arrange_data_for_model(
    model::T where T <: Resonance.UnivariatePredictor,
    input_df::AbstractDataFrame,
    input_cols:: R where R <: Union{UnitRange, Vector{Int64}};
    print_cols=true
    )

    expected_predictors = names(model.inputs_outputs[1])
    n_expected_predictors = length(expected_predictors)

    received_predictors = names(input_df[!, input_cols])
    n_received_predictors = length(received_predictors)

    names_intersection = received_predictors .∈ Ref(expected_predictors)

    @info "Provided data sanity check:
    You provided $(n_received_predictors) data columns.
    The model expects to receive $(n_expected_predictors).
    From the columns provided, $(sum(names_intersection)) correspond to columns expected by the model.
    $(n_received_predictors - sum(names_intersection)) provided columns were not recognized by the model and will not be considered.
    $(n_expected_predictors - sum(names_intersection)) expected columns were not provided and will be constructed with zeros."

    if print_cols
        @info "Argument print_cols was provided. Printing column names for debug.\nReceived columns matching the model expectations:"
        println(received_predictors[names_intersection])
    end

    arranged_df = DataFrame()
    
    for colname in expected_predictors
        if colname ∈ received_predictors
            insertcols!(arranged_df, colname => input_df[:, colname])
        else
            insertcols!(arranged_df, colname => 0.0)
        end #endif
    end # end for colname

    return arranged_df

end

input_cols = eval(Meta.parse(parsed_args[:input_cols]))

arranged_df = arrange_data_for_model(
    model,
    input_df,
    input_cols
)

predictions = Resonance.predict(model, arranged_df)

#####
# 4. Outputting predictions
#####

output_df = select(input_df, Not(input_cols))

if ncol(output_df) == 0
    output_df = DataFrame(
        :sampleidx => 1:length(predictions),
        :prediction => predictions
    )
else
    insertcols!(output_df, :prediction => predictions)
end

CSV.write(pwd()*"/"*parsed_args[:output_filename], output_df)