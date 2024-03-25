using Resonance
using CairoMakie
using Distances
using MultivariateStats
using Random
using Distributions
using PERMANOVA

comm = Resonance.load_raw_metaphlan()
comm = filter(f-> taxrank(f) === :species, comm)
comm = comm[:, (1:502)[Not([263, 419])]] # 263 and 419 had no bugs

mat = collect(abundances(comm))
notzero = vec(sum(mat; dims=2)) .> 0
mat = mat[notzero, :]

comm2 = CommunityProfile(vcat(mat, reshape([zeros(250); rand(250) .* 10], 1, 500)),
			 [features(comm)[notzero]; [Taxon("diagnosis_bool")]],
			 samples(comm)
)
relativeabundance!(comm2)

comm3 = CommunityProfile(vcat(mat, reshape(range(0,10; length=500), 1, 500)),
			 [features(comm)[notzero]; [Taxon("diagnosis_continuous")]],
			 samples(comm)
)
relativeabundance!(comm3)

corspearman(vec(abundances(comm2[r"diagnosis", :])), repeat([1,2]; inner=250))
permanova(DataFrame(risk = repeat([0,1]; inner=250)),
	  abundances(comm2)',
	  BrayCurtis,
	  @formula(1~risk),
	  1000
)

permanova(DataFrame(risk = collect(1:500)),
	  abundances(comm3)',
	  BrayCurtis,
	  @formula(1~risk),
	  1000
)


#-

fig = Figure(; size = (1200, 1200))

ax1 = Axis(fig[1,1];
	   xlabel = "diagnosis",
	   xticks = ([1,2], ["no", "yes"]),
	   ylabel="bug abundance"
)
ax2 = Axis(fig[2,1];
	   xlabel = "risk",
	   ylabel="bug abundance"
)
ax3 = Axis(fig[1,2];)
ax4 = Axis(fig[2,2];)

scatter!(ax1, repeat([1,2]; inner=250) .+ rand(Normal(0., 0.05), 500),
	 vec(abundances(comm2[r"diagnosis", :]))
)
scatter!(ax2, 1:500, vec(abundances(comm3[r"diagnosis", :])))

p1 = plot_pcoa!(ax3, fit(MDS, Microbiome.braycurtis(comm2); distances=true);
    color=repeat([:blue, :red]; inner=250)
)

p2 = plot_pcoa!(ax4, fit(MDS, Microbiome.braycurtis(comm3); distances=true);
    color=1:500
)

Legend(fig[1,3],
    [MarkerElement(; color=c, marker=:circle) for c in [:blue, :red]],
    ["no", "yes"], "diagnosis"
)
Colorbar(fig[2,3], colorrange=(1,500); label="risk")

Label(fig.layout[1, 1, TopLeft()], "A",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
Label(fig.layout[2, 1, TopLeft()], "B",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
Label(fig.layout[1, 2, TopLeft()], "C",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
Label(fig.layout[1, 2, TopLeft()], "D",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)

save("/home/kevin/Downloads/hypo_bugs.png", fig)


#-
permanova(DataFrame(risk = repeat([0,1]; inner=250)),
	  abundances(comm2[[collect(1:4:20)..., size(comm2, 1)], :])',
	  BrayCurtis,
	  @formula(1~risk),
	  1000
)

permanova(DataFrame(risk = collect(1:500)),
	  abundances(comm3)',
	  BrayCurtis,
	  @formula(1~risk),
	  1000
)

fig2 = Figure(; size=(900,500))
ax1 = Axis(fig2[1,1])
ax2 = Axis(fig2[1,2])

p1 = plot_pcoa!(ax1, fit(MDS, Microbiome.braycurtis(
			relativeabundance(
			    comm2[[collect(20:40)..., size(comm2, 1)], :])
			); distances=true);
    color=repeat([:blue, :red]; inner=250)
)

p2 = plot_pcoa!(ax2, fit(MDS, Microbiome.braycurtis(
			relativeabundance(
			    comm2[[collect(20:40)..., size(comm2, 1)], :])
			); distances=true);
    color=1:500
)

Legend(fig2[1,3],
    [MarkerElement(; color=c, marker=:circle) for c in [:blue, :red]],
    ["no", "yes"], "diagnosis"
)
Colorbar(fig[2,3], colorrange=(1,500); label="risk")

save("/home/kevin/Downloads/hypo_bugs2.png", fig2)

