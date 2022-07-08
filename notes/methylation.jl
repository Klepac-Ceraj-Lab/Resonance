using Resonance

methylation = CSV.read(datafiles("targets_passed.csv"), DataFrame)
wrangled = CSV.read(datafiles("wrangled.csv"), DataFrame)

jnd = leftjoin(methylation, wrangled, on=["Sample_ID"=>"subject", "Timepoint"=>"timepoint"], matchmissing=:notequal)

jnd.Sample_Name = jnd.sample

files = filter(f-> any(s-> !ismissing(s) && contains(f, s), jnd.sample), readdir("/augusta/echo/analysis/biobakery3/links/metaphlan/"))
open(datafiles("methylation_files.txt"), "w") do io
    println.(Ref(io), files)
end

met = metaphlan_profiles(joinpath.(Ref("/augusta/echo/analysis/biobakery3/links/metaphlan/"), files))
CSV.write("/home/kevin/Desktop/methylation_taxa.csv", met)
CSV.write("/home/kevin/Desktop/methylation_taxa_link.csv", select(jnd, ["Sample_ID", "Timepoint", "Sample_Name"]))