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

# Loading the raw ECHO metadata
echo_raw_metadata = Resonance.load_raw_metadata()

echo_taxonomic_profiles = @chain Resonance.load_raw_metaphlan() begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
end

rename_ref_table = CSV.read("oldname_newname_reltable.csv", DataFrame; stringtype = String)

echo_processed_inputs = @chain echo_raw_metadata begin
    select(["subject", "timepoint", "ageMonths", "cogScore", "hhs", "education", "race", "sex", "omni", "omni_collectionDate", "etoh", "etoh_collectionDate"])
    sort([:subject, :timepoint])
    innerjoin(rename_ref_table, on = :omni => :old_tubeid; matchmissing = :notequal)
    innerjoin(echo_taxonomic_profiles, on = :new_seqid => :sample_base; matchmissing = :notequal)
    select!(Not([ :file ]))
end

echo_processed_inputs = echo_processed_inputs[ findall([ sum(j) for j in eachrow(echo_processed_inputs[:,15:end]) ] .> 0.0), vcat(1:8..., findall([ any( k .> 0.0) for k in eachcol(echo_processed_inputs[:,15:end]) ]) .+ 8) ]

CSV.write("ECHO_metabolomics_preselectedinputs.csv", echo_processed_inputs)