using Resonance

humann_join("/grace/echo/analysis/biobakery3/links/humann/genefamilies", 
            "/grace/echo/analysis/biobakery3/links/humann/all_genefamilies.tsv"; 
            file_name="genefamilies")

func = humann_profiles("/grace/echo/analysis/biobakery3/links/humann/all_genefamilies.tsv", stratified=false)

sns = map(samplenames(func)) do s
    replace(s, r"_S\d+_genefamilies"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

func = CommunityProfile(abundances(func)[:, idx], features(func), sns[idx])

CSV.write("data/genefamilies.csv", func)