include("startup.jl")

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

fig, ax, p = pie(demo, color = colors, axis=(;aspect=1))
hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1], [MarkerElement(color=c, marker=:circle) for c in colors],
                labels, "Race",
        orientation=:horizontal,
        tellheight=true, tellwidth=false, nbanks=3
)

fig
##

fig, ax, p = pie(demo, color = Makie.to_colormap(:viridis, 6), axis=(;aspect=1))
hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1], [MarkerElement(color=c, marker=:circle) for c in Makie.to_colormap(:viridis, 6)],
                ["Junior high school",
                 "Some high school",
                 "High school grad",
                 "Some college",
                 "College grad",
                 "Grad/professional school"
                ], "Maternal Education",
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