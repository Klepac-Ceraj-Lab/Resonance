using Resonance
omni, tps = startup([:omni, :tps])

using CairoMakie
using AlgebraOfGraphics

rce = (@chain tps begin
    groupby(:subject)
    combine(:simple_race => (r-> coalesce(r...)) => :race)
end).race


hist(collect(skipmissing(tps.ageMonths)))

hhs = (@chain tps begin
    groupby(:subject)
    combine(:mother_HHS_Education => (r-> coalesce(r...)) => :hhs)
end).hhs


##

labels = [
    "White",
    "Black",
    "Asian",
    "Mixed",
    "Other",
    "Unknown"
]

counts = [
    1142,
    290,
    1,
    343,
    30,
    594
]

colors = [
    :navajowhite3,
    :coral,
    :darkorange3,
    :darkseagreen3,
    :cyan4,
    :dimgray
]

fig, ax, p = pie(counts, color = colors, axis=(;aspect=1))
hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1], [MarkerElement(color=c, marker=:circle) for c in colors],
                labels, "Race",
        orientation=:horizontal,
        tellheight=true, tellwidth=false, nbanks=3
)

fig
##

using CategoricalArrays

ed = skipmissing(unique(tps, [:subject, :mother_HHS_Education]).mother_HHS_Education) |> collect
filter!(!=(-8), ed)
ed = categorical(sort(ed); ordered=true)
ed = recode(ed, 
     2 => "Junior high school",
     3 => "Some high school",
     4 => "High school grad",
     5 => "Some college",
     6 => "College grad",
     7 => "Grad/professional school")

fig, ax, p = pie(map(levels(ed)) do lev
        count(==(lev), ed)
    end;
    color = Makie.to_colormap(:viridis, 6), axis=(;aspect=1)
)
hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1],
    [MarkerElement(color=c, marker=:circle) for c in Makie.to_colormap(:viridis, 6)],
    levels(ed),
    "Maternal Education";
    orientation=:horizontal,
    tellheight=true, tellwidth=false, nbanks=3
)

fig

##

unique!(omni, :sample)

dyads = @chain omni begin
    groupby(:subject)
    combine(:Mother_Child => (mc-> length(unique(mc)) == 2) => :hasdyad)
    @rsubset!(:hasdyad)
end

toplot = @chain omni begin
    @rsubset(:subject in dyads.subject)
    groupby(:subject)
    combine(:Mother_Child=> (mc-> count(==("C"), mc))=> :nchild, :Mother_Child=> (mc -> count(==("M"), mc)) => :nmother)
end

fig, ax1, h1 = hist(toplot.nchild, axis=(; xlabel="number of child samples", ylabel="number of dyads"))
ax2 = Axis(fig[1,2], xlabel="number of mother samples")
h2 = hist!(ax2, toplot.nmother)

linkyaxes!(ax1, ax2)
fig