# #Metadata wrangling

# ## Get sample-specific metadata from airtable database

using Resonance

samplemeta = airtable_metadata() # having set ENV["AIRTABLE_KEY"]

# Now to add important metadata to samples

samples = MicrobiomeSample[]

let fields = [
    "timepoint",
    "subject",
    "sid_old",
    "CovidCollectionNumber",
    "Mgx_batch",
    "MaternalID"]

    for row in eachrow(samplemeta)
        startswith(row.sample, "FE") && continue # skip ethanol samples
        s = MicrobiomeSample(row.sample)
        set!(s, NamedTuple(Symbol(k) => v for (k, v) in pairs(row[Not(:sample)]) if !ismissing(v)))
        push!(samples, s)
    end
end

# ## Metadata stored in filemaker pro database
#
# Read in the different tables as `DataFrame`s,
# then normalize certain columns.
# First, samples

fmp_sample = DataFrame(XLSX.readtable("data/Sample_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_sample, Dict(:SampleID=> :sample, :studyID=>:subject, :collectionNum=> :timepoint))
subset!(fmp_sample, :sample=> ByRow(s-> !startswith(s, "FE")))

# Then, subject-specific data

fmp_subject = DataFrame(XLSX.readtable("data/Subject_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_subject, Dict(:studyID=>:subject))

# Then, timepoint-specific data

fmp_timepoint = DataFrame(XLSX.readtable("data/Timepoint_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_timepoint, Dict(:studyID=>:subject))

# and COVID-specific samples

fmp_covid = DataFrame(XLSX.readtable("data/COVID_Fecal_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_covid, Dict(:studyID=>:subject))

# ## Getting data joined together

fmp_timed = outerjoin(fmp_timepoint,
                      subset(fmp_sample, :timepoint=> ByRow(!ismissing)), # filter out covid samples
                      on=[:subject, :timepoint])

fmp_alltp = leftjoin(fmp_timed, fmp_subject, on=[:subject])
unique!(fmp_alltp, [:subject, :timepoint])

# Add info about brain data (just if it's there)

brain = let 
    fcleft = brain_ingest("data/freesurfer_curvature_leftHemi_oct2021.csv"; label="curvature_leftHemi")
    fcright = brain_ingest("data/freesurfer_curvature_rightHemi_oct2021.csv"; label="curvature_rightHemi")
    ftleft = brain_ingest("data/freesurfer_thickness_leftHemi_oct2021.csv"; label="thickness_leftHemi")
    ftright = brain_ingest("data/freesurfer_thickness_rightHemi_oct2021.csv"; label="thickness_rightHemi")
    seg = brain_ingest("/home/kevin/Repos/Resonance/data/segmentationVolumeMeasurements_oct2021.csv"; label="segmentation")

    outerjoin(fcleft, fcright, ftleft, ftright, seg, on=[:subject, :timepoint], makeunique=true)
end

## Validation

for row in eachrow(brain)
    all(h-> !ismissing(h) && h, row[r"has_"]) || @info row[r"has_"]
end

fmp_alltp = leftjoin(fmp_alltp, brain, on=[:subject, :timepoint])
sort!(fmp_alltp, [:subject, :timepoint])

# ## Getting values for presence of data types
#
# We want to plot set intersections for having various kinds of data.
# In some cases, additional wrangling is necessary

fmp_alltp.simple_race = map(fmp_alltp.simple_race) do r
    ismissing(r) && return missing
    r == "Unknown" && return missing
    r == "Decline to Answer" && return missing
    r ∈ ("Mixed", "Mixed Race") && return "Mixed"
    return r
end

unique!(fmp_alltp, [:subject, :timepoint])

fmp_alltp.has_race = .!ismissing.(fmp_alltp."simple_race")

##

subj = groupby(fmp_alltp, :subject)
transform!(subj, nrow => :n_samples; ungroup=false)

transform!(subj, AsTable(r"subject|breast"i) => (s -> begin
    if any(!ismissing, s.breastFedPercent)
        return fill(true, length(s[1]))
    else
        return fill(false, length(s[1]))
    end
    fill(0, length(s[1]))
end) => :has_bfperc_subj; ungroup=false)

transform!(subj, AsTable(r"subject|breast"i) => (s -> begin
    if any(!ismissing, s.breastFedPercent)
        return fill(true, length(s[1]))
    else
        return fill(false, length(s[1]))
    end
    fill(0, length(s[1]))
end) => :has_bfperc_subj; ungroup=false)

fmp_alltp = transform(subj, AsTable([:timepoint, :sample]) => (s -> begin
    has_stool = .!ismissing.(s.sample)
    has_prevstool = fill(false, length(s[1]))
    for i in 1:length(s[1])
        i == 1 && continue
        any(!ismissing, s.sample[1:i-1]) && (has_prevstool[i] = true)
    end
    (; has_stool, has_prevstool)
end)=> [:has_stool, :has_prevstool])

fmp_alltp.has_everbreast = .!ismissing.(fmp_alltp."everBreastFed")
fmp_alltp.has_bfperc = .!ismissing.(fmp_alltp."breastFedPercent")

# ## Plotting set intersections


##

using CairoMakie

ys = ["scan", "segmentation", "previous stool", "stool", "breastfeeding"]
ycols = [:has_curvature_leftHemi, :has_segmentation, :has_prevstool, :has_stool, :has_bfperc]

##

fig = Figure()

intersection_ax = Axis(fig[1,1], ylabel="intersection size")
dot_ax = Axis(fig[2,1])
set_ax = Axis(fig[2,2], xlabel="set size")

intersects = [
    [1],        # 1. scan only    
    [4],        # 2. stool only
    [1,4],      # 3. scan & stool
    [1,2,4],    # 4. segmentation & stool
    [1,2,3],    # 5. segmentation & prior stool
    [1,2,3,4],  # 6. segmentation & stool & prior stool
    [1,4,5],    # 7. scan & stool & bf
    [1,2,4,5],  # 8. segmentation & stool & bf
    [1,2,3,5],  # 9. segmentation & prior stool & bf
    [1,2,3,4,5] #10. segmentation & stool & prior stool & bf
]

barplot!(intersection_ax, 1:10, [count_set(df, ycols, i) for i in intersects],
         bar_labels=:y, flip_labels_at=400)

barplot!(set_ax, 1:5, [
    count(x-> !ismissing(x) && x, fmp_alltp[!, col]) for col in ycols
], direction=:x, bar_labels=:y, flip_labels_at=1000)

upset_dots!(dot_ax, [[1,], [4,], [1,4], [1,2,4], [1,2,3], [1,2,3,4], [1,4,5], [1,2,4,5], [1,2,3,5], [1,2,3,4,5]])

for i in 1:2:length(ycols)
    poly!(dot_ax,
     BBox(0, length(intersects) + 1, i-0.5, i+0.5),
     color = :gray95
    )
 end

hidexdecorations!(count_ax)
hideydecorations!(setsize_ax)

save("figures/upset.pdf", fig)
fig


##

ycols = [:a,:b,:c,:d]
df = DataFrame(rand(Bool, 1000, 4), ycols)

fig = Figure()
intersection_ax = Axis(fig[1,1], ylabel="intersection size")
dot_ax = Axis(fig[2,1])
set_ax = Axis(fig[2,2], xlabel="set size")

intersects = [
    [1],       # 1. a only    
    [2],       # 2. b only
    [3],       # 3. c only
    [4],       # 4. d only
    [1,2],     # 5. a and b only
    [1,3],     # 6. a and c only
    [1,2,3],   # 7. a, b, and c only
    [1,2,3,4], # 8. a, b, c, and d
]

barplot!(intersection_ax, 1:length(intersects),
                          [count_set(df, ycols, i) for i in intersects],
         bar_labels=:y)

barplot!(set_ax, 1:4, [count(df[!, col]) for col in ycols], 
                 direction=:x, bar_labels=:y)

upset_dots!(dot_ax, intersects)
hidexdecorations!(intersection_ax)
hideydecorations!(set_ax)
fig