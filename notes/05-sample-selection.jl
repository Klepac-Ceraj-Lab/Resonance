using Resonance
using CategoricalArrays

omni = CSV.read(datafiles("wrangled", "omnisamples.csv"), DataFrame)
etoh = CSV.read(datafiles("wrangled", "etohsamples.csv"), DataFrame)
tps  = CSV.read(datafiles("wrangled", "timepoints.csv"), DataFrame)

DataFrames.transform!(groupby(tps, :subject), :mother_HHS_Education => (r->coalesce(r...)) => :hhs)

tps.education = let
    ed = categorical(tps.hhs; levels=[-8 , 2:7...], ordered=true)
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

DataFrames.transform!(groupby(tps, :subject), :race => (r->coalesce(r...)) => :race)

tps.race = let
    race = categorical(tps.race; ordered=true)
    race = recode(race, 
        "American Indian or Alaska Native"=> "Other",
        "Some other race"                 => "Other",
        "Asian Indian"                    => "Asian",
        "Black or African American"       => "Black",
        missing                           => "Unknown"
    )
    droplevels!(race)
    levels!(race, ["White","Black","Asian","Mixed","Other","Unknown"])
    race
end


select!(tps, "subject", "timepoint",  
                "ageMonths",
                "childGender" => "sex",
                "race",
                "education",
                "cogScore",
                "cogScorePercentile",
                "has_segmentation",
                "ECHOTPCoded",
                "assessmentDate",
                "scanDate"
)



subset!(omni, "Mgx_batch" => ByRow(!ismissing))
subset!(etoh, "Metabolomics_batch" => ByRow(m-> m==1))

tpstoolset = collect(zip(omni.subject, omni.timepoint))
tpmetabset = collect(zip(etoh.subject, etoh.timepoint))

subj = @chain tps begin
    groupby(:subject)
    transform!(
        :cogScore => ByRow(!ismissing) => :has_cogScore,
        AsTable(r"subject|timepoint"i) => ByRow(s -> begin
            (s.subject, s.timepoint) in tpstoolset
        end) => :has_stool,
        AsTable(r"subject|timepoint"i) => ByRow(s -> begin
            (s.subject, s.timepoint) in tpmetabset
        end) => :has_metabolomics;
        ungroup=false
    )
    transform!(AsTable(r"timepoint|has_stool") => (s-> begin
        findprevstool(s.timepoint, s.has_stool)
    end) => :has_prevstool,
               AsTable(r"timepoint|has_stool") => (s-> begin
        findprevstool(s.timepoint, s.has_stool; rev=true)
    end) => :has_futstool; ungroup=false)
end

tps.has_scan = .!ismissing.(tps.scanDate)


subjs = combine(groupby(tps, :subject), 
    :has_stool=> any => :has_stool,
    :has_scan=> any => :has_scan,
    :has_cogScore => any => :has_cogScore,
    :has_segmentation => any => :has_segmentation,
    :has_metabolomics => any => :has_metabolomics,
)

# ## Plotting set intersections
# 
# Load CairoMakie for plotting.

using CairoMakie

# ### Intersections between stool and brain scans

##

ys = ["stool", "cogScore", "segmentation", "metabolomics"]
ycols = [:has_stool, :has_cogScore, :has_segmentation, :has_metabolomics]


intersects = [
    [1],
    [2],
    [1,2],
    [1,3],
    [1,4],
    [1,2,3], 
    [1,2,3,4], 
]

fig = Resonance.plot_upset(tps, ycols, ys, intersects)

save("figures/upset_cog_stool.png", fig)
fig

#-

fig = Resonance.plot_upset(subjs, ycols, ys, intersects)
save("figures/upset_subj_cog_stool.png", fig)
fig

# ## Sample Selection
#
# We want to keep all subjects that have at least 1 timepoint
# with both a stool sample and a cognitive assessment.

transform!(groupby(tps, :subject), AsTable([:has_stool, :has_cogScore]) => (tbl -> any(tbl.has_stool) && any(tbl.has_cogScore)) => :keeper)

keeptps = subset(tps, :keeper => identity, 
                      :ECHOTPCoded => ByRow(tp-> !startswith(tp, "Pre") # remove moms ("prenatal")
))


# Now, add columns with sample names for omnigene (metagenomics) and ethanol (metabolomics).
# Note, up above we removed samples that haven't (yet) been sequenced or LCMS'ed

keepomni = unique(select(omni, ["subject", "timepoint", "sample"]), ["subject","timepoint"])
rename!(keepomni, "sample"=> "omni")
keepetoh = unique(select(etoh, ["subject", "timepoint", "sample"]), ["subject","timepoint"])
rename!(keepetoh, "sample"=> "etoh")

leftjoin!(keeptps, keepomni; on = ["subject", "timepoint"])
leftjoin!(keeptps, keepetoh; on = ["subject", "timepoint"])

CSV.write("data/timepoints_final.csv", keeptps)

# ## Upset plots on kept samples / subjects
#
# Now that we've got our subjects, lets do the same more visualizations.

keepsubj = groupby(keeptps, :subject)

#-
# ### Number of samples pers subject with different data types:

fig = Figure()
ax1 = Axis(fig[1,1]; xlabel = "N samples", ylabel="N subjects", title="total timepoints", xticks = 0:13)
ax2 = Axis(fig[1,2]; xlabel = "N samples", ylabel="N subjects", title="stool samples", xticks = 0:13)
ax3 = Axis(fig[2,1]; xlabel = "N samples", ylabel="N subjects", title="cogScores", xticks = 0:13)
ax4 = Axis(fig[2,2]; xlabel = "N samples", ylabel="N subjects", title="segmentation", xticks = 0:13)

let toplot = combine(keepsubj, :subject => length => :toplot).toplot
    barplot!(ax1, unique(toplot), [count(==(c), toplot) for c in unique(toplot)]; bins = length(unique(toplot)))
end
let toplot = combine(keepsubj, :has_stool => count => :toplot).toplot
    barplot!(ax2, unique(toplot), [count(==(c), toplot) for c in unique(toplot)]; bins = length(unique(toplot)))
end
let toplot = combine(keepsubj, :has_cogScore => count => :toplot).toplot
    barplot!(ax3, unique(toplot), [count(==(c), toplot) for c in unique(toplot)]; bins = length(unique(toplot)))
end
let toplot = combine(keepsubj, :has_segmentation => count => :toplot).toplot
    barplot!(ax4, unique(toplot), [count(==(c), toplot) for c in unique(toplot)]; bins = length(unique(toplot)))
end

axs = [ax1, ax2, ax3, ax4]
linkxaxes!(axs...)
foreach(ax-> tightlimits!(ax, Left(), Right(), Bottom()), axs)

save("figures/keep_samplecounts.png")
fig

#-
# ### Upset plots (intersections)


ys = ["stool", "cogScore", "segmentation", "metabolomics"]
ycols = [:has_stool, :has_cogScore, :has_segmentation, :has_metabolomics]


intersects = [
    [1],
    [2],
    [1,2],
    [1,3],
    [1,4],
    [1,2,3], 
    [1,2,3,4], 
]


fig = Resonance.plot_upset(keeptps, ycols, ys, intersects)
Label(fig[0, 2], "Timepoints intersections"; textsize = 20)

save("figures/upset_keep_cog_stool.png", fig)
fig
