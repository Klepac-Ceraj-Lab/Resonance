using Resonance

tps = CSV.read("data/wrangled/timepoints.csv", DataFrame)
fss = CSV.read("data/wrangled/omnisamples.csv", DataFrame)
subtpset = collect(zip(fss.subject, fss.timepoint))


subj = @chain tps begin
    groupby(:subject)
    transform!(
        nrow => :n_timepoints,
        :cogScore => ByRow(!ismissing) => :has_cogScore,
        AsTable(r"subject|breast"i) => (s -> begin
            if any(!ismissing, s.breastFedPercent)
                return fill(true, length(s[1]))
            else
                return fill(false, length(s[1]))
            end
            fill(0, length(s[1]))
        end) => :has_bfperc_subj,
        AsTable(r"subject|timepoint"i) => ByRow(s -> begin
            (s.subject, s.timepoint) in subtpset
        end) => :has_stool;
        ungroup=false
    )
    transform!(AsTable(r"timepoint|has_stool") => (s-> begin
        findprevstool(s.timepoint, s.has_stool)
    end) => :has_prevstool,
               AsTable(r"timepoint|has_stool") => (s-> begin
        findprevstool(s.timepoint, s.has_stool; rev=true)
    end) => :has_futstool; ungroup=false)
end


# ## Plotting set intersections
# 
# Load CairoMakie

using CairoMakie

# ### Intersections between stool and brain scans

##

ys = ["scan",  "stool", "prev stool"]
ycols = [:has_segmentation, :has_stool, :has_prevstool]

fig = Figure()

intersection_ax = Axis(fig[1,1:2], ylabel="intersection size", yautolimitmargin = (0, 0.15))
dot_ax = Axis(fig[2,1:2], yticklabelpad = 10, yticks = (1:length(ys), ys))
set_ax = Axis(fig[2,3], xlabel="set size", xticklabelrotation=π/4,
                 xautolimitmargin = (0, 0.25),  xgridvisible = false)

intersects = [
    [1],        # 1. scan only    
    [2],        # 2. stool only
    [1,2],      # 3. scan & stool
    [1,3],      # 4. scan & prior stool
    [2,3],      # 5. stool & prior stool
    # [2,4],      # 5. stool & future stool
    [1,2,3],    # 5. scan & stool & prior stool
    # [1,2,5],    # 6. scan & stool & bf
    # [1,2,3,5]   # 7. scan & stool & prior stool && bf
]

barplot!(intersection_ax, 1:length(intersects), [count_set(tps, ycols, i) for i in intersects],
            bar_labels=:y, color = :gray20,
            label_size = 14, label_formatter = x -> string(Int(x)))

barplot!(set_ax, 1:length(ys), [count(x-> !ismissing(x) && x, tps[!, col]) for col in ycols], 
            direction=:x, bar_labels=:y, label_size = 14, color = :gray20,
            label_formatter = x -> string(Int(x)))

            
for i in 1:2:length(ycols)
    poly!(dot_ax,
    BBox(0, length(intersects) + 1, i-0.5, i+0.5),
    color = :gray95
    )
end

upset_dots!(dot_ax, intersects)

hidexdecorations!(intersection_ax)
hideydecorations!(set_ax)

rowgap!(fig.layout, 0)
linkyaxes!(dot_ax, set_ax)
linkxaxes!(dot_ax, intersection_ax)
hidespines!(intersection_ax, :t, :r, :b)
hidespines!(set_ax, :t, :r, :l)

save("figures/upset_scan_stool.pdf", fig)
fig

##

# ### Age distributions
#
# First, some histograms

using Dates

fss.collectionDate = [ismissing(d) ? missing : Date(d) for d in fss.collectionDate]
colldates = sort(fss.collectionDate |> skipmissing |> collect)

drange = range(Date(2017,07,01), today(), step=Month(1))

colldatecount = map(drange) do dt
    count(cd-> dt <= cd < dt+Month(1), colldates)
end

assdates = sort(tps.assessmentDate |> skipmissing |> collect)

assdatecount = map(drange) do dt
    count(cd-> dt <= cd < dt+Month(1), assdates)
end

##

fig = Figure()
stool_age = Axis(fig[1,1], title="Age of stool samples", xlabel="age (years)", ylabel="count")
scan_age = Axis(fig[1,2], title="Age of scans", xlabel="age (years)", ylabel="count")
stool_age_young = Axis(fig[2,1], title="Age of stool samples", xlabel="age (months)", ylabel="count")
scan_age_young = Axis(fig[2,2], title="Age of scans", xlabel="age (months)", ylabel="count")
for ax in [stool_age, scan_age, stool_age_young, scan_age_young]
    tightlimits!(ax, Left(), Right(), Bottom())
end

hist!(stool_age, collect(skipmissing(subset(tps, :has_stool=> ByRow(x-> !ismissing(x) && x)).ageMonths)) ./ 12, color=:gray20)
hist!(scan_age, collect(skipmissing(subset(tps, :has_segmentation=> ByRow(x-> !ismissing(x) && x)).ageMonths)) ./ 12, color=:gray20)

hist!(stool_age_young, subset(tps, AsTable([:ageMonths, :has_stool]) => ByRow(
                                row -> !ismissing(row.ageMonths) && row.ageMonths < 24 && !ismissing(row.has_stool) && row.has_stool)
                                ).ageMonths |> skipmissing |> collect,
                                color=:gray20)


hist!(scan_age_young, subset(tps, AsTable([:ageMonths, :has_segmentation]) => ByRow(
                                row -> !ismissing(row.ageMonths) && row.ageMonths < 24 && !ismissing(row.has_segmentation) && row.has_segmentation)
                                ).ageMonths |> skipmissing |> collect,
                                color=:gray20)

Label(fig[1, 0], text = "All kids", tellheight=false)
Label(fig[2, 1], text = "Young kids", tellheight=false)

save("figures/age_hists.pdf")
fig

##
fig = Figure()

stools_bydate = Axis(fig[1, 1], title="Stool samples collected by date",
                        xticks=(1:4:length(drange), string.(drange[1:4:end])),
                        xticklabelrotation=π/4)
barplot!(stools_bydate, 1:length(drange), colldatecount, color=:gray20)
vlines!(stools_bydate, [33], linestyle=:dash, color=:gray20)
tightlimits!(stools_bydate, Left(), Right(), Bottom())

ass_bydate = Axis(fig[2, 1], title="Assessments by date",
                        xticks=(1:4:length(drange), string.(drange[1:4:end])),
                        xticklabelrotation=π/4)
barplot!(ass_bydate, 1:length(drange), assdatecount, color=:gray20)
vlines!(ass_bydate, [33], linestyle=:dash, color=:gray20)
tightlimits!(ass_bydate, Left(), Right(), Bottom())

linkyaxes!(stools_bydate, ass_bydate)

fig
##

@info "Collected in 2019: $(count(c-> Year(c) == Year(2019), colldates))"
@info "Collected in 2020: $(count(c-> Year(c) == Year(2020), colldates))"
@info "Collected in 2021: $(count(c-> Year(c) == Year(2021), colldates))"

##

on_subj = combine(subj, :has_segmentation => (s-> any(skipmissing(s))) => :any_segmentation,
                        :has_stool        => (s-> any(skipmissing(s))) => :any_stool,
                        AsTable([:has_segmentation, :has_prevstool]) => (tab->
                            any(row-> all(t-> !ismissing(t) && t, row), tab)) => :any_scan_stool)

##

ys = ["scan",  "stool", "scan + prev stool"]
ycols = [:any_segmentation, :any_stool, :any_scan_stool]

intersects = [
    [1],
    [2],
    [1,2],
    [1,3],
    [1,2,3]
]


fig = Figure()

intersection_ax = Axis(fig[1,1:2], ylabel="intersection size", yautolimitmargin = (0, 0.15), title="By subject")
dot_ax = Axis(fig[2,1:2], yticklabelpad = 10, yticks = (1:3, ys))
set_ax = Axis(fig[2,3], xlabel="set size", xticklabelrotation=π/4,
                 xautolimitmargin = (0, 0.25),  xgridvisible = false)

barplot!(intersection_ax, 1:5, [count_set(on_subj, ycols, i) for i in intersects],
            bar_labels=:y, color = :gray20,
            label_size = 14, label_formatter = x -> string(Int(x)))

barplot!(set_ax, 1:3, [count(x-> !ismissing(x) && x, on_subj[!, col]) for col in ycols], 
            direction=:x, bar_labels=:y, label_size = 14, color = :gray20,
            label_formatter = x -> string(Int(x)))

            
for i in 1:2:length(ycols)
    poly!(dot_ax,
    BBox(0, length(intersects) + 1, i-0.5, i+0.5),
    color = :gray95
    )
end

upset_dots!(dot_ax, intersects)

hidexdecorations!(intersection_ax)
hideydecorations!(set_ax)

rowgap!(fig.layout, 0)
linkyaxes!(dot_ax, set_ax)
linkxaxes!(dot_ax, intersection_ax)
hidespines!(intersection_ax, :t, :r, :b)
hidespines!(set_ax, :t, :r, :l)

save("figures/upset_scan_stool_bysybject.pdf", fig)
fig

##

# slope charts

fig = Figure(resolution=(800, 1800))
ax = Axis(fig[1,1], xlabel="age (years)", ydectorations=false)

let y = 1
    df = sort(tps, :n_samples)
    gdf = groupby(df, :subject)
    for sub in gdf
        age = sub.ageMonths ./ 12
        col = map(s-> s ? :royalblue4 : :darkgray, sub.has_stool)
        idx = .!ismissing.(age)
        xs = collect(skipmissing(age))
        scatter!(ax, xs, fill(y, sum(idx)), color=col[idx])
        length(xs) > 1 && lines!(ax, [extrema(xs)...], [y, y], color=:black, linewidth=0.5)
        y += 1
    end
end

fig

plot(rand(10), rand(10))


##

tps.has_cogScore = .!ismissing.(tps.cogScore)
ys = ["cogScore",  "stool", "prev stool"]
ycols = [:has_segmentation, :has_stool, :has_prevstool]

fig = Figure()

intersection_ax = Axis(fig[1,1:2], ylabel="intersection size", yautolimitmargin = (0, 0.15))
dot_ax = Axis(fig[2,1:2], yticklabelpad = 10, yticks = (1:length(ys), ys))
set_ax = Axis(fig[2,3], xlabel="set size", xticklabelrotation=π/4,
                 xautolimitmargin = (0, 0.25),  xgridvisible = false)

intersects = [
    [1],        # 1. score only    
    [2],        # 2. stool only
    [1,2],      # 3. score & stool
    [1,3],      # 4. score & prior stool
    [1,2,3],      # 4. score & prior stool
]

barplot!(intersection_ax, 1:length(intersects), [count_set(tps, ycols, i) for i in intersects],
            bar_labels=:y, color = :gray20,
            label_size = 14, label_formatter = x -> string(Int(x)))

barplot!(set_ax, 1:length(ycols), [count(x-> !ismissing(x) && x, tps[!, col]) for col in ycols], 
            direction=:x, bar_labels=:y, label_size = 14, color = :gray20,
            label_formatter = x -> string(Int(x)))

            
for i in 1:2:length(ycols)
    poly!(dot_ax,
    BBox(0, length(intersects) + 1, i-0.5, i+0.5),
    color = :gray95
    )
end

upset_dots!(dot_ax, intersects)

hidexdecorations!(intersection_ax)
hideydecorations!(set_ax)

rowgap!(fig.layout, 0)
linkyaxes!(dot_ax, set_ax)
linkxaxes!(dot_ax, intersection_ax)
hidespines!(intersection_ax, :t, :r, :b)
hidespines!(set_ax, :t, :r, :l)

save("figures/upset_score_stool.pdf", fig)
fig

##


fig, ax, plt = hist(collect(skipmissing(tps.cogScore)))
plt2 = hist!(ax, collect(skipmissing(tps.cogScore[tps.has_stool .|| tps.has_prevstool])))
Legend(fig[1,2], [plt, plt2], ["All scores", "Scores with stool"])
fig

##

tps.under3 = tps.ageMonths .< 36
tps.under5 = tps.ageMonths .< 60

ys = ["cogScore",  "stool", "prev stool", "under 3"]
ycols = [:has_segmentation, :has_stool, :has_prevstool, :under3]

fig = Figure()

intersection_ax = Axis(fig[1,1:2], ylabel="intersection size", yautolimitmargin = (0, 0.15))
dot_ax = Axis(fig[2,1:2], yticklabelpad = 10, yticks = (1:length(ys), ys))
set_ax = Axis(fig[2,3], xlabel="set size", xticklabelrotation=π/4,
                 xautolimitmargin = (0, 0.25),  xgridvisible = false)

intersects = [
    [1],          # 1. score only    
    [2],          # 2. stool only
    [1,2],        # 3. score & stool
    [1,2,4],      # 4. score & stool & under 3
    [1,3,4],      # 4. score & prior stool & under 3
]

barplot!(intersection_ax, 1:length(intersects), [count_set(tps, ycols, i) for i in intersects],
            bar_labels=:y, color = :gray20,
            label_size = 14, label_formatter = x -> string(Int(x)))

barplot!(set_ax, 1:length(ycols), [count(x-> !ismissing(x) && x, tps[!, col]) for col in ycols], 
            direction=:x, bar_labels=:y, label_size = 14, color = :gray20,
            label_formatter = x -> string(Int(x)))

            
for i in 1:2:length(ycols)
    poly!(dot_ax,
    BBox(0, length(intersects) + 1, i-0.5, i+0.5),
    color = :gray95
    )
end

upset_dots!(dot_ax, intersects)

hidexdecorations!(intersection_ax)
hideydecorations!(set_ax)

rowgap!(fig.layout, 0)
linkyaxes!(dot_ax, set_ax)
linkxaxes!(dot_ax, intersection_ax)
hidespines!(intersection_ax, :t, :r, :b)
hidespines!(set_ax, :t, :r, :l)

save("figures/upset_score_stool_under3.pdf", fig)
fig

##

tps.under5 = tps.ageMonths .< 60

ys = ["cogScore",  "stool", "prev stool", "under 5"]
ycols = [:has_segmentation, :has_stool, :has_prevstool, :under5]

fig = Figure()

intersection_ax = Axis(fig[1,1:2], ylabel="intersection size", yautolimitmargin = (0, 0.15))
dot_ax = Axis(fig[2,1:2], yticklabelpad = 10, yticks = (1:length(ys), ys))
set_ax = Axis(fig[2,3], xlabel="set size", xticklabelrotation=π/4,
                 xautolimitmargin = (0, 0.25),  xgridvisible = false)

intersects = [
    [1],          # 1. score only    
    [2],          # 2. stool only
    [1,2],        # 3. score & stool
    [1,2,4],      # 4. score & stool & under 5
    [1,3,4],      # 4. score & prior stool & under 5
]

barplot!(intersection_ax, 1:length(intersects), [count_set(tps, ycols, i) for i in intersects],
            bar_labels=:y, color = :gray20,
            label_size = 14, label_formatter = x -> string(Int(x)))

barplot!(set_ax, 1:length(ycols), [count(x-> !ismissing(x) && x, tps[!, col]) for col in ycols], 
            direction=:x, bar_labels=:y, label_size = 14, color = :gray20,
            label_formatter = x -> string(Int(x)))

            
for i in 1:2:length(ycols)
    poly!(dot_ax,
    BBox(0, length(intersects) + 1, i-0.5, i+0.5),
    color = :gray95
    )
end

upset_dots!(dot_ax, intersects)

hidexdecorations!(intersection_ax)
hideydecorations!(set_ax)

rowgap!(fig.layout, 0)
linkyaxes!(dot_ax, set_ax)
linkxaxes!(dot_ax, intersection_ax)
hidespines!(intersection_ax, :t, :r, :b)
hidespines!(set_ax, :t, :r, :l)

save("figures/upset_score_stool_under5.pdf", fig)
fig

## 

tps.under6mo = tps.ageMonths .<= 6

fig, ax, hst = hist(@rsubset(tps, :under6mo).ageMonths, axis=(; xlabel = "Age (months)"))

ax2, hst2 = hist(fig[1,2], @rsubset(tps, :under6mo, :breastfeeding=="exclusive formula").ageMonths)
hist!(ax2, @rsubset(tps, :under6mo, :breastfeeding=="exclusive breast").ageMonths)
fig

##

(@chain tps begin
    groupby(:subject)
    @combine(:breastfeeding = unique(:breastfeeding))
end).breastfeeding |> countmap