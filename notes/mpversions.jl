using Resonance
using CairoMakie
using Statistics
using Arrow
using SparseArrays
using VKCComputing
## Data Loading

base = LocalBase()
mdata = Resonance.load(Metadata())
seqs = subset(mdata, "sample"=> ByRow(!ismissing)) 

seqs.sample = [s for s in seqs.sample]
seqs.edfloat = map(x-> ismissing(x) ? missing : Float64(levelcode(x)), seqs.education)
taxa = Resonance.load(TaxonomicProfiles(); timepoint_metadata = seqs) # this can take a bit
taxdf = comm2wide(taxa)

olddf = let 
	tbl = Arrow.Table("input/taxa_old.arrow")
	mdt = Arrow.getmetadata(tbl)
	fs = [taxon(line) for line in eachline(IOBuffer(mdt["features"]))]
	mdt = DataFrame(
           sample = [MicrobiomeSample(line) for line in eachline(IOBuffer(mdt["samples"]))],
           read_depth = map(l-> l=="missing" ? missing : parse(Float64, l), eachline(IOBuffer(mdt["reads"]))),
    )
	transform!(mdt, "sample"=> ByRow(n-> replace(name(n), r"_S\d+"=> ""))=> "sample_base")
	mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
	comm = CommunityProfile(mat, fs, mdt.sample)
	set!(comm, mdt)
	comm2wide(comm)
end

olddf.seqid = map(olddf.sample) do samp
	samp == "FG02471_S89" && return "SEQ02512"
	sb, swell = match(r"(FG\d+)_(S\d+)",samp)
	biosp = base["Biospecimens", String(sb)]
	sid = get(biosp, :seqprep, nothing)
	isnothing(sid) && error("no seqprep for $sb")
	srecs = base[sid]

	srec = filter(sr-> get(sr, :S_well, nothing) == swell, srecs)
	length(srec) != 1 && @info "sample $samp" srecs
	first(srec)[:uid]
end

olddf.sample_base = replace.(olddf.sample, r"_S\d+"=> "")
rename!(olddf, "sample_base" => "biospecimen", "sample"=> "biospecimen_swell", "seqid"=> "sample")

taxdf_new = rename(taxdf, Dict(x => string(x, "_new") for x in names(taxdf, r"^[A-Z][a-z]+_(sp_(CAG_)?[\d_]+|[a-z]+)")))
leftjoin!(olddf, taxdf_new; on="sample", makeunique=true)

#-
function squarecoords(iter; invert = false)
	nside = ceil(Int, √length(iter))
	Iterators.map(enumerate(iter)) do (p, x)
		i = p - nside * floor(Int, (p - 1) / nside)
		j = 1 + floor(Int, (p - 1) / nside)
		return invert ? (x, j, i) : (x, i, j)
	end
end
#-

bugs = [
	("Alistipes_finegoldii", 1),
	("Bacteroides_dorei", 3),
	("Parasutterella_excrementihominis", 3),
	("Roseburia_faecis", 1),
	("Ruminococcus_gnavus", 1),
	("Dorea_longicatena", 1),
	("Blautia_wexlerae", 1),
	("Alistipes_putredinis", 1),
	("Adlercreutzia_equolifaciens", 1),
	("Eubacterium_ramulus", 3),
	("Roseburia_hominis", 1),
	("Gordonibacter_pamelaeae", 1),
	("Eubacterium_eligens", 3),
	("Faecalibacterium_prausnitzii", 3),
	("Asaccharobacter_celatus", 3),
	("Sutterella_wadsworthensis", 2),
	("Barnesiella_intestinihominis", 2),
	("Alistipes_obesi", 2),
	("Eggerthella_sp_CAG_298", 2),
	("Bacteroides_vulgatus", 3)
]

bugcoords = squarecoords(bugs) |> collect

#-
	
fig = Figure(; resolution = (1600, 1600));

for (bug, i, j) in bugcoords
	(bug, imp) = bug
	ax = Axis(fig[i, j]; xlabel = "old", ylabel="new", title=bug)
	newbug = isempty(names(olddf, Regex(string(bug, "_new")))) ?
						fill(0, size(olddf, 1)) :
						olddf[!, string(bug, "_new")]
	bugab = isempty(names(olddf, Regex(string(bug, raw"$")))) ?
						fill(0, size(olddf, 1)) :
						olddf[!, bug]
	scatter!(ax, bugab, newbug; color = (:dodgerblue, :slateblue1, :goldenrod)[imp])
end

Legend(fig[1:ceil(Int,√(length(bugs))), 0],
		[MarkerElement(; marker=:rect, color) for color in (:dodgerblue, :slateblue1, :goldenrod)],
		["Important old", "Important new", "Important both"]
		)

save("manuscript/assets/newmp.png", fig)
fig


