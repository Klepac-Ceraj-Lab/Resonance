using Resonance
using FileIO
using CairoMakie # for plotting
using Distances
using MultivariateStats
using CategoricalArrays
using Statistics
using Chain

#-

mdata = Resonance.load(Metadata())
echo_tx = Resonance.load_raw_metaphlan()
species = let
    tax = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)
    tax = tax[:, get(tax, :ageMonths) .< 6]
    tax = tax[vec(prevalence(tax) .> 0), :]
end

echo_phy = filter(t-> taxrank(t) == :phylum, echo_tx[:, samplenames(species)])

khuladf = CSV.read("khula_profiles.tsv", DataFrame; comment="#")
khuladf.taxon = taxon.(last.(split.(khuladf.clade_name, '|')))
ss = MicrobiomeSample.(replace.(names(khuladf, r"ZSA"), r"_S\d+_profile"=> ""))
kcomm = filter(t-> taxrank(t) == :species, CommunityProfile(Matrix(khuladf[!, r"ZSA"]), khuladf.taxon, ss))

for s in Microbiome.samples(species)
    set!(s, :cohort, "ECHO")
end

for s in Microbiome.samples(kcomm)
    set!(s, :cohort, "Khula")
end

allcomm = commjoin(species, kcomm)

spedm = Resonance.braycurtis(allcomm)
spepco = fit(MDS, spedm; distances=true)

#-

df = DataFrame(species = featurenames(allcomm), prevalence=vec(prevalence(allcomm)),
                                                 meanab = map(s-> mean([a for a in abundances(allcomm[s, :]) if a != 0]), features(allcomm)),
                                                 prevalence_echo= vec(prevalence(allcomm[:, get(allcomm, :cohort) .== "ECHO"])),
                                                 prevalence_khula=vec(prevalence(allcomm[:, get(allcomm, :cohort) .== "Khula"])),
                                                 meanab_echo = map(s-> mean([a for a in abundances(allcomm[s, get(allcomm, :cohort) .== "ECHO"]) if a != 0]), features(allcomm)),
                                                 meanab_khula = map(s-> mean([a for a in abundances(allcomm[s, get(allcomm, :cohort) .== "Khula"]) if a != 0]), features(allcomm)),
                                                 meanab_echo_all = map(s-> mean(abundances(allcomm[s, get(allcomm, :cohort) .== "ECHO"])), features(allcomm)),
                                                 meanab_khula_all = map(s-> mean(abundances(allcomm[s, get(allcomm, :cohort) .== "Khula"])), features(allcomm)),
)

df.meanab_echo = [isnan(x) ? 0 : x for x in df.meanab_echo]
df.meanab_khula = [isnan(x) ? 0 : x for x in df.meanab_khula]
df.meanab_echo_all = [isnan(x) ? 0 : x for x in df.meanab_echo_all]
df.meanab_khula_all = [isnan(x) ? 0 : x for x in df.meanab_khula_all]

pretty_table(first(sort(df, "meanab_khula"; rev = true), 25); alignment=:l)
pretty_table(first(sort(df, "meanab_echo"; rev = true), 25); alignment=:l)
pretty_table(first(sort(df, "prevalence_khula"; rev = true), 25); alignment=:l)
pretty_table(first(sort(df, "prevalence_echo"; rev = true), 25); alignment=:l)
@chain df begin
    subset(AsTable(["prevalence_echo", "prevalence_khula"]) => ByRow(row-> any(==(0), values(row))))
    select(["species", "prevalence_echo", "prevalence_khula", "meanab_echo", "meanab_khula"])
    sort("prevalence_khula"; rev=true)
    @aside pretty_table(first(_, 10); alignment=:l)    
    sort("prevalence_echo"; rev=true)
    @aside pretty_table(first(_, 10); alignment=:l)
end
#-

figure = Figure()
ax = Axis(figure[1,1])

pc = Resonance.plot_pcoa!(ax, spepco; color = [c == "ECHO" ? :orange : :purple for c in get(allcomm, :cohort)])
figure

scatter(df.meanab_echo_all, df.meanab_khula_all; axis=(; xlabel="ECHO mean abundance (where present)", ylabel="Khula mean abundance (where present)"))
scatter(df.meanab_echo_all, df.meanab_khula_all; axis=(; xlabel="ECHO mean abundance", ylabel="Khula mean abundance"))

CSV.write("data/khula_echo.csv", sort(df, :meanab_khula); rev=true)

#-

let sdf = sort(df, :prevalence_khula; rev=true)
    scatter(sdf.prevalence_echo, sdf.prevalence_khula; axis=(; xlabel="ECHO prevalence", ylabel="Khula prevalence"))
    annotations!(first(sdf.species, 20), Point2f.(zip(first(sdf.prevalence_echo, 20), first(sdf.prevalence_khula, 20))))
    current_figure()
end


let sdf = sort(df, :prevalence_khula; rev=true)
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel = "Bug", ylabel = "prevalence", xticks = (1:10, first(sdf, 10).species))
    ax.xticklabelrotation = π / 4
    barplot!(ax, repeat(1:10; outer=2), [first(sdf.prevalence_khula,10); first(sdf.prevalence_echo, 10)];
        dodge = repeat([1,2], inner=10), color=dodge = repeat([1,2], inner=10))
    fig
end