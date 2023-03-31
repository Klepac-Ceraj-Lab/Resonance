using Resonance
using CairoMakie

mdata = Resonance.load(Metadata())
ecs = Resonance.load(ECProfiles(); timepoint_metadata = mdata)
relativeabundance!(ecs)

dihyrofolate = [
    "6.3.2.12",
    "1.5.1.3"
]

pcresol = "1.17.9.1"

dhf1 = vec(abundances(ecs[r"6\.3\.2\.12", :]))
dhf2 = vec(abundances(ecs[r"1\.5\.1\.3:", :]))
# ecs[r"1\.17\.9\.1", :] # not found

scatter(dhf1, get(ecs, :cogScore); axis = (; xlabel = "abundance dihydrofolate synthetase (6.3.2.12)", ylabel="cogScore"))
scatter(dhf2, get(ecs, :cogScore); axis = (; xlabel = "abundance dihydrofolate reductase (1.5.1.3)", ylabel="cogScore"))


ecdf = comm2wide(ecs[r"(6\.3\.2\.12:|1\.5\.1\.3:)", :])

runlms(ecdf[ecdf.filter_00to120, :], tablefiles("lms_ecs_dihydrofolate.csv"), last(names(ecdf), 2))