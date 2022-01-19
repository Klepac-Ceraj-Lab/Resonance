using Resonance
using DataFramesMeta

wrangled = CSV.read("data/wrangled.csv", DataFrame)
unique(wrangled, [:subject, :timepoint])

subj = groupby(wrangled, :subject)

transform!(subj, nrow => :n_samples; ungroup=false)

transform!(subj, AsTable(r"subject|breast"i) => (s -> begin
    if any(!ismissing, s.breastFedPercent)
        return fill(true, length(s[1]))
    else
        return fill(false, length(s[1]))
    end
    fill(0, length(s[1]))
end) => :has_bfperc_subj; ungroup=false)

# ## Plotting set intersections
# 
# Load CairoMakie

using CairoMakie

# ### Intersections between stool and brain scans

##

ys = ["scan",  "stool", "prev stool", "breastfeeding"]
ycols = [:has_segmentation, :has_stool, :has_prevstool, :has_bfperc]

fig = Figure()

intersection_ax = Axis(fig[1,1:2], ylabel="intersection size", yautolimitmargin = (0, 0.15))
dot_ax = Axis(fig[2,1:2], yticklabelpad = 10, yticks = (1:4, ys))
set_ax = Axis(fig[2,3], xlabel="set size", xticklabelrotation=π/4,
                 xautolimitmargin = (0, 0.25),  xgridvisible = false)

intersects = [
    [1],        # 1. scan only    
    [2],        # 2. stool only
    [1,2],      # 3. scan & stool
    [1,3],      # 4. scan & prior stool
    [2,3],      # 5. stool & prior stool
    [1,2,3],    # 5. scan & stool & prior stool
    [1,2,4],    # 6. scan & stool & bf
    [1,2,3,4]   # 7. scan & stool & prior stool && bf
]

barplot!(intersection_ax, 1:8, [count_set(wrangled, ycols, i) for i in intersects],
            bar_labels=:y, color = :gray20,
            label_size = 14, label_formatter = x -> string(Int(x)))

barplot!(set_ax, 1:4, [count(x-> !ismissing(x) && x, wrangled[!, col]) for col in ycols], 
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

colldates = sort(wrangled.collectionDate |> skipmissing |> collect)

drange = range(Date(2017,07,01), today(), step=Month(1))

colldatecount = map(drange) do dt
    count(cd-> dt <= cd < dt+Month(1), colldates)
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

hist!(stool_age, collect(skipmissing(subset(wrangled, :has_stool=> ByRow(x-> !ismissing(x) && x)).ageMonths)) ./ 12, color=:gray20)
hist!(scan_age, collect(skipmissing(subset(wrangled, :has_segmentation=> ByRow(x-> !ismissing(x) && x)).ageMonths)) ./ 12, color=:gray20)

hist!(stool_age_young, subset(wrangled, AsTable([:ageMonths, :has_stool]) => ByRow(
                                row -> !ismissing(row.ageMonths) && row.ageMonths < 24 && !ismissing(row.has_stool) && row.has_stool)
                                ).ageMonths |> skipmissing |> collect,
                                color=:gray20)


hist!(scan_age_young, subset(wrangled, AsTable([:ageMonths, :has_segmentation]) => ByRow(
                                row -> !ismissing(row.ageMonths) && row.ageMonths < 24 && !ismissing(row.has_segmentation) && row.has_segmentation)
                                ).ageMonths |> skipmissing |> collect,
                                color=:gray20)

Label(fig[1, 0], text = "All kids", tellheight=false)
Label(fig[2, 1], text = "Young kids", tellheight=false)

stools_bydate = Axis(fig[3, 1:3], title="samples collected by date",
                        xticks=(1:4:length(drange), string.(drange[1:4:end])),
                        xticklabelrotation=π/4)
barplot!(stools_bydate, 1:length(drange), colldatecount, color=:gray20)
vlines!(stools_bydate, [33], linestyle=:dash, color=:gray20)
tightlimits!(stools_bydate, Left(), Right(), Bottom())


save("figures/age_hists.pdf")
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
    df = sort(wrangled, :n_samples)
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
