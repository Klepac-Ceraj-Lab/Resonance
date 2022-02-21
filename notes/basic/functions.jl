using Resonance

# humann_join("/grace/echo/analysis/biobakery3/links/humann/genefamilies", 
#             "/grace/echo/analysis/biobakery3/links/humann/all_genefamilies.tsv"; 
#             file_name="genefamilies")

# run(pipeline(`cat "/grace/echo/analysis/biobakery3/links/humann/all_genefamilies.tsv"`,
#              `grep -v '|'`;
#              stdout="/grace/echo/analysis/biobakery3/links/humann/all_genefamilies_unstrat.tsv"
#             )
# )

CSV.read("/grace/echo/analysis/biobakery3/links/humann/all_genefamilies_unstrat.tsv", DataFrame)

sns = map(samplenames(func)) do s
    replace(s, r"_S\d+_genefamilies"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

func = CommunityProfile(abundances(func)[:, idx], features(func), sns[idx])

CSV.write("data/genefamilies.csv", func)