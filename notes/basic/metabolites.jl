using Resonance
using DataFramesMeta
using CairoMakie

samplemeta = airtable_metadata()
allmeta = CSV.read("data/wrangled.csv", DataFrame)

subset!(allmeta, :sample=>ByRow(!ismissing))
# transform!(allmeta, :sid_old => ByRow(id-> ismissing(id) ? id : replace(id, r"_(\d+)F_"=>s"_\1E_")) => :sid_old_etoh)
# select(allmeta, r"sid_old")

# Note: had to delete empty "C1162_3E_1A" "M1162_3E_1A" columns from C8-pos (they were duplicated)
#
# Note: "M0855_1E_1A" missing from HILIC pos and neg

c18neg = load_metabolites("data/21_0924_VKC_C18-neg_results.xlsx")
c8pos = load_metabolites("data/21_0924_VKC_C8-pos_results.xlsx")
hneg = load_metabolites("data/21_0924_VKC_HILIC-neg_results.xlsx")
hpos = load_metabolites("data/21_0924_VKC_HILIC-pos_results.xlsx")

for df in [c18neg, c8pos, hneg, hpos]
    select!(df, Not(r"^PREF"))
end

c8pos = c8pos[:, names(c18neg)]
hneg = hneg[:, names(c18neg)]
hpos = hpos[:, names(c18neg)]

metabs = vcat(c18neg, c8pos, hneg, hpos)
metabs[!, :uid] = [join(p, "_") for p in zip(metabs.Method, metabs.Compound_ID)]
select!(metabs, :uid, Cols(:))



# rename to new sample ids
metab_samples = names(metabs)[9:end]
rename_dict = Dict((old=>new for (old, new) in zip(samplemeta.sid_old, samplemeta.sample)))
rename_dict["C1162_3E_1A"] = "FE01852" # should have been 4E
rename_dict["C1227_3E_1A"] = "FE01922" # should have been 4E
rename_dict["M0932_7E_1A"] = "FE01759" # should have been C0932

rename_dict = Dict(k => rename_dict[k] for k in keys(rename_dict) if k in metab_samples)
rename!(metabs, rename_dict)
metab_samples = names(metabs)[9:end]



subset!(samplemeta, :sample => ByRow(!ismissing) )
subset!(samplemeta, :sample => ByRow(s-> s in metab_samples))
subset!(allmeta, :sample => ByRow(s-> s in metab_samples))
unique!(allmeta, :sample)
sort!(allmeta, :sample)

metabs = metabs[:, [(1:8)..., (sortperm(names(metabs)[9:end]) .+ 8)...]]
@assert names(metabs)[9:end] == allmeta.sample


# normalize
fill_vals = select(metabs, AsTable(propertynames(metabs)[9:end]) => ByRow(minimum∘skipmissing) => :fill_vals).fill_vals .÷ 2
metabs[!, 9:end] .= coalesce.(metabs[!, 9:end], fill_vals)

meds = median.(eachcol(select(metabs, names(metabs)[9:end])))
globmed = median(meds)

for i in 9:size(metabs, 2)
    metabs[!, i] ./= (meds[i-8] / globmed)
end

CSV.write("data/metabolites.csv", metabs)

## 

fig = Figure()
ax1 = Axis(fig[1,1])
hist!(ax1, meds)
vlines!(ax1, [globmed])

fig


##


##

hist(collect(skipmissing(allmeta.ageMonths)))
metabs.Metabolite[findall(m-> !ismissing(m) && contains(m, r"sero|dopa|glut|buty|epine|anand|trypto|aspart|olea|kynu"i), metabs.Metabolite)] |> sort