#####
# Table 5. Experimental design and input composition for Random Forest experiments
#####

using Resonance
using CategoricalArrays
using Statistics
using Chain
using JLD2
using MLJ
using DataFrames.PrettyTables

isdir(tablefiles()) || mkpath(tablefiles())

## Building the Table

table5 = DataFrame(
    :input_set => 1:10,
    :age_bracket => vcat( repeat([ "0 to 6 months" ], 5), repeat( [ "18 to 120 months" ], 5)),
    :microbial_feature => repeat([ "Not provided", "Taxonomic profile", "Taxonomic profile", "Functional Profile (ECs)", "Functional Profile (ECs)" ], 2),
    :demo_provided => repeat([ "yes", "no", "yes", "no", "yes" ], 2),
)

CSV.write(tablefiles("Table5.csv"), table5)

pretty_table(table5;
    header = [ "Input set", "Age bracket", "Microbiome encoding type", "Demographics Provided? (sex, education)" ],
    backend = Val(:latex)
)