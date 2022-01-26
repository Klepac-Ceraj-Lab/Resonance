using Resonance
using AlgebraOfGraphics

timepoints = CSV.read("data/wrangled/timepoints.csv", DataFrame)

omnisamples = CSV.read("data/wrangled/omnisamples.csv", DataFrame)
rename!(omnisamples, "DOC"=>"date")

covid = CSV.read("data/wrangled/covid.csv", DataFrame)
rename!(covid, ["SampleID"=>"sample", "CollectionDate"=>"date"])

spec = metaphlan_profiles("data/wrangled/species.csv")
spec = spec[:, intersect(samplenames(spec), omnisamples.sample)]

insert!(spec, select(omnisamples, ["sample", "date"]))

insert!(spec, leftjoin(
                    select(omnisamples, ["sample", "subject", "timepoint"]),
                    select(timepoints, ["subject", "timepoint", "ageMonths"]);
                    on = ["subject", "timepoint"])
)

metadf = metadata(spec) |> DataFrame
spec = spec[:, .!ismissing.(metadf.ageMonths)]

p = pcoa(spec)

using Dates
using CategoricalArrays

covdates = [d < Date(2020, 03, 01) ? "pre" : "post" for d in Date.(get(spec, :date))]

df = DataFrame(loadings(p)[:, 1:6], [:pco1, :pco2, :pco3, :pco4, :pco5, :pco6])
df.sample = samplenames(spec)
df.date = Date.(get(spec, :date))
df.covdates = categorical([d < Date(2020, 03, 01) ? "pre" : "post" for d in df.date])
df.ageMonths = get(spec, :ageMonths)

plt = data(df) * mapping(:pco1=> "PCo1 ($(round(varexplained(p)[1] * 100; digits=2))) %", 
                   :pco2=> "PCo2 ($(round(varexplained(p)[2] * 100; digits=2))) %")

draw(plt * mapping(color=:covdates))

plt = data(df) * mapping(:pco3=> "PCo3 ($(round(varexplained(p)[3] * 100; digits=2))) %", 
                         :pco4=> "PCo4 ($(round(varexplained(p)[4] * 100; digits=2))) %")

draw(plt * mapping(color=:covdates))

plt = data(df) * mapping(:pco5=> "PCo5 ($(round(varexplained(p)[5] * 100; digits=2))) %", 
                         :pco6=> "PCo6 ($(round(varexplained(p)[6] * 100; digits=2))) %")

draw(plt * mapping(color=:covdates))
