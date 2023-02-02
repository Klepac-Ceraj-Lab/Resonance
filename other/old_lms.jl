using Resonance
using GLM
using Statistics
using CategoricalArrays
using MultipleTesting
using Dates

mdata = Resonance.load(Metadata())
species = Resonance.load(TaxonomicProfiles(); timepoint_metadata = mdata)

specdf = comm2wide(species)
specdf.quartile = categorical(let
    l, u = quantile(specdf.cogScore, [0.25, 0.75])
    map(s-> s < l ? "lower" : s > u ? "upper" : "middle", specdf.cogScore)
end; levels=["lower", "middle", "upper"], ordered = true)

non_spec_cols = [
    "sample", "subject", "timepoint", "ageMonths", "sex", "race", "education", "date",
    "cogScore", "sample_base", "read_depth", "filter_00to120", "filter_00to06", "filter_18to120", "quartile"
    ]

stotals = map(i-> sum(specdf[i, Not(non_spec_cols)]), 1:nrow(specdf))
for spc in names(specdf, Not(non_spec_cols))
    specdf[!, spc] ./= stotals
end

mdata_raw = Resonance.load_raw_metadata()
species_raw = Resonance.load_raw_metaphlan()

let md = select(subset(mdata_raw, "omni" => ByRow(!ismissing)), "subject", "timepoint", "omni"=>"sample_base", "ageMonths", "cogScore", "hhs", "omni_collectionDate"=>"collectionDate")
    td = DataFrame("sample"=> samplenames(species_raw), "sample_base"=> get(species_raw, :sample_base))
    set!(species_raw, leftjoin(td, md; on="sample_base"))
end

specdf_raw = let spec = species_raw[taxrank.(features(species_raw)) .== :species, .!ismissing.(get(species_raw, :subject))]
    df = DataFrame(
        sample         = get(spec, :sample_base),
        subject        = get(spec, :subject),
        timepoint      = get(spec, :timepoint),
        collectionDate = get(spec, :collectionDate),
        ageMonths      = get(spec, :ageMonths),
        cogScore       = get(spec, :cogScore),
        education      = get(spec, :hhs),
        read_depth     = get(spec, :read_depth)
    )

    for f in features(spec)
        df[!, name(f)] = vec(abundances(spec[f,:])) ./ 100
    end
    df.education = let
        ed = categorical(df.education; levels=[-8 , 2:7...], ordered=true)
        ed = recode(ed,
        -8 => missing,
        2 => "Junior high school",
        3 => "Some high school",
        4 => "High school grad",
        5 => "Some college",
        6 => "College grad",
        7 => "Grad/professional school")
        ed
    end
    df
end

specdf_raw = subset(specdf_raw, "ageMonths"=> ByRow(!ismissing), "cogScore"=> ByRow(!ismissing))
specdf_raw.cogScore = parse.(Float64, specdf_raw.cogScore)

function getsubset(df, col, pred; minprev = 0.15)
    newdf = subset(df, col=> ByRow(pred))
    keep = findall(col-> !(eltype(col) <: Float64) || prevalence(col) >=  minprev, eachcol(newdf))
    return newdf[:, keep]
end

#- 

speco18_full = unique(getsubset(specdf_raw, "ageMonths", a-> 18 <= a; minprev=0.15), :subject)
speco18_120 = unique(getsubset(specdf_raw, "ageMonths", a-> 18 <= a < 120; minprev=0.15), :subject)

speco18_full_sorted = unique(getsubset(sort(specdf_raw, :timepoint), "ageMonths", a-> 18 <= a; minprev=0.15), :subject)
speco18_120_sorted = unique(getsubset(sort(specdf_raw), "ageMonths", a-> 18 <= a < 120; minprev=0.15), :subject)

spec_new = specdf[specdf.filter_18to120, :]
speco18_120_unified = let subtp = Set(Set(t for t in zip(spec_new.subject, spec_new.timepoint)))
    speco18_full_sorted[map(row-> (row.subject, row.timepoint) in subtp, eachrow(speco18_full_sorted)), :]
end
#-

raw4 = let indf = speco18_full_sorted
    dropmissing!(indf, "read_depth")
    lmresults = DataFrame()

    for spc in names(indf, Not(["subject", "timepoint", "collectionDate", "ageMonths", "cogScore", "sample", "read_depth", "education"]))    
        @info spc
        over0 = indf[!, spc] .> 0
        ab = indf[!, spc]

        sum(over0) / size(indf, 1) > 0.25 || continue
        # ab = log2.(ab .+ (minimum(filter(>(0), ab))/2))
        ab = asin.(sqrt.(ab))

        df = indf[:, ["cogScore", "education", "ageMonths", "read_depth"]]
        df.bug = ab

        mod = lm(@formula(bug ~ cogScore + education + ageMonths + read_depth), df; dropcollinear=false)
        ct = DataFrame(coeftable(mod))
        ct.species .= spc
        ct.kind .= "cogScore"
        append!(lmresults, ct)

    end
    select!(lmresults, Cols(:species, :Name, :))
    rename!(lmresults, "Pr(>|t|)"=>"pvalue");

    @chain lmresults begin
        subset!(:Name => ByRow(x->
            !any(y-> contains(x, y), 
                ("(Intercept)", "ageMonths", "read_depth", "education")
                )
            )
        )

        groupby(:kind)
        transform!(:pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
        sort!(:qvalue)
        select(Not(["Lower 95%", "Upper 95%"]))
    end
end
#-

newdf = let 
    indf = specdf[specdf.filter_18to120, :]
    # outfile = tablefiles("lms_species_18to120.csv")
    lmresults = DataFrame()

    for spc in names(indf, Not(non_spec_cols))
        @info spc  
        over0 = indf[!, spc] .> 0
        ab = indf[!, spc]
        # ab = log2.(ab .+ (minimum(filter(>(0), ab))/2))
        ab = asin.(sqrt.(ab))
        
        sum(over0) / size(indf, 1) > 0.15 || continue
        spc in setdiff(newdf.species, raw.species) && continue

        df = indf[:, ["cogScore", "education", "ageMonths", "read_depth"]]
        df.bug = ab

        mod = lm(@formula(bug ~ cogScore + ageMonths + read_depth), df; dropcollinear=false)
        ct = DataFrame(coeftable(mod))
        ct.species .= spc
        ct.kind .= "cogScore"
        append!(lmresults, ct)
    end

    select!(lmresults, Cols(:species, :Name, :))
    rename!(lmresults, "Pr(>|t|)"=>"pvalue");

    @chain lmresults begin
        subset!(:Name => ByRow(x->
            !any(y-> contains(x, y), 
                ("(Intercept)", "ageMonths", "read_depth", "education")
                )
            )
        )

        groupby(:kind)
        transform!(:pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
        sort!(:qvalue)
    end

    # CSV.write(outfile, lmresults)
    lmresults[:, Not(["Lower 95%", "Upper 95%"])]
end



#-

mysetdiff(df1, df2) = setdiff(Set(t for t in zip(df1.subject, df1.timepoint)), Set(t for t in zip(df2.subject, df2.timepoint)))

mysetdiff(spec_new, speco18_120_sorted)
mysetdiff(spec_new, speco18_120_unified)
mysetdiff(speco18_120_unified, spec_new)


#- 

using CairoMakie

scatter(newdf."Coef.", -1 .* log2.(newdf.qvalue), color = map(q-> q < 0.2 ? :orange : :dodgerblue, newdf.qvalue); 
        axis=(; xlabel="Coeficient", ylabel="-log2(qvalue)"))

scatter(raw."Coef.", -1 .* log2.(raw.qvalue), color = map(q-> q < 0.2 ? :orange : :dodgerblue, raw.qvalue); 
        axis=(; xlabel="Coeficient", ylabel="-log2(qvalue)"))