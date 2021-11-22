using Resonance

metadata = CSV.read("data/wrangled.csv", DataFrame)
@rsubset!(metadata, !ismissing(:subject), 
                    !ismissing(:sample))

@rsubset!(metadata, !startswith(:sample, "C"),
                    !startswith(:sample, "z"))


volumes = metadata[:, r"^(left|right|sample)"i][:, 2:end]
@rsubset!(volumes, !ismissing(Symbol("Left-Lateral-Ventricle")))
CSV.write("data/volumes.csv", volumes)

met = metaphlan_profiles(filter(f-> contains(f, "profile"), readdir("/grace/echo/biobakery3/links/metaphlan")))