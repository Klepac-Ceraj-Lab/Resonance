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

Legend(fig[2,1], [MarkerElement(color=c, marker=:circle) for c in colors,
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