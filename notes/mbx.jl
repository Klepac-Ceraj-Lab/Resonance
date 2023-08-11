using Resonance
using XLSX
using ThreadsX
using Chain

biospecimens = CSV.read("input/Biospecimens-ECHO.csv", DataFrame)
olddb = CSV.read("input/Samples-olddb.csv", DataFrame)

biospdict = Dict{String,String}()
for row in eachrow(biospecimens)
    ismissing(row.aliases) && continue
    if row.aliases isa Vector
        for alias in row.aliases
            biospdict[alias] = row.uid
        end
    else
        biospdict[row.aliases] = row.uid
    end
end



#-

mbxhave = readlines("mbx_samples.txt")
filter!(m-> !startswith(m, "PRE"), mbxhave)

#-

mbxdf = DataFrame(sid = map(mbxhave) do mbx
    startswith(mbx, "FE") && return mbx
    if startswith(mbx, "SD")
        m = match(r"SD(\d+)([A-Z])", mbx)
        isnothing(m) && error("SD sample $mbx doesn't match pattern")
        sid = lpad(m[1], 4, "0")
        tp = findfirst(c-> "$c" == m[2], "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        mbx = string("C", sid, "_", tp, "E", "_1A")
    end
    mbxid = get(biospdict, mbx, nothing)
    isnothing(mbxid) && return mbx
    return mbxid
end)

leftjoin!(mbxdf, select(biospecimens, "uid"=>"sid", "subject", "collection"=> "timepoint"), on=:sid, makeunique=true)



transform!(mbxdf, AsTable(["sid", "subject", "timepoint"])=>ByRow(row-> begin
    !ismissing(row.subject) && return row
    m = match(r"[CM](\d+)_([0-9]+)E_\d+\w", row.sid)
    if isnothing(m)
        sid = row.sid == "FE2790" ? "FE02790" : row.sid

        ix = findfirst(id -> !ismissing(id) && id == sid, olddb.sample)
        isnothing(ix) && return (; sid, subject=missing, timepoint=missing)
        return (; sid, subject= string("resonance-", lpad(olddb.subject[ix], 4, "0")), timepoint=olddb.timepoint[ix])
    else
        return (; sid = row.sid, subject=string("resonance-", lpad(m[1], 4, "0")), timepoint=parse(Int, m[2]))
    end
end)=> ["sid", "subject", "timepoint"])

#- 

transform!(mbxdf, "subject"=> ByRow(s-> parse(Int, replace(s, "resonance-" => "")))=> "subject")

mdata = Resonance.load(Metadata())
unique(mdata.subject) |> length

count(s-> s ∈ mdata.subject, mbxdf.subject)

mdata = leftjoin!(mdata, unique(select(mbxdf, "subject", "timepoint", "sid"=> "mbxsample")), on=["subject", "timepoint"])

subset(mdata, "mbxsample" => ByRow(!ismissing)) |> size
unique(subset(mdata, "mbxsample" => ByRow(!ismissing)), "subject") |> size

subset(mdata, "mbxsample" => ByRow(!ismissing)).ageMonths |> describe

#- OK, you win, Claire 😭

function fixcol(col, name)
    name == "MZ"                ? convert(Vector{Float64}, col) :
    name == "RT"                ? convert(Vector{Float64}, col) :
    name == "Method"            ? fill(first(col) == "HIL-pos" ? "HILIC-pos" : first(col), length(col)) :
    name == "Compound_ID"       ? convert(Vector{String}, col) :
    name == "HMDB_ID"           ? convert(Vector{Union{Missing,String}}, col) :
    name == "Metabolite"        ? convert(Vector{Union{Missing,String}}, col) :
    name == "HMDB_ID_Certainty" ? map(col) do c
        c isa Int && return c
        ismissing(c) && return c
        c == "NA" && return 0
    end : coalesce.(col, 0)
end
#-

mbx_xlsx_files = filter(f-> endswith(f, ".xlsx") && !contains(f, r"centric"i), readdir("input/", join=true))

mbx_long = mapreduce(vcat, mbx_xlsx_files) do file
    xf = XLSX.readxlsx(file)
    sh = xf[first(XLSX.sheetnames(xf))]
    mat = sh[:]
    startrow = findfirst(!ismissing, mat[:,1])
    endcol = let i = findfirst(ismissing, mat[startrow, :])
        isnothing(i) ? size(mat, 2) : i - 1
    end

    df = DataFrame(mat[startrow+1:size(mat, 1), 1:endcol], vec(mat[startrow, 1:endcol]); makeunique=true)
    rename!(df, first(names(df, r"HMDB_ID_")) => "HMDB_ID_Certainty")
    for n in names(df)
        df[!, n] = fixcol(df[!, n], n)
    end

    df = subset(stack(df, Not(["MZ", "RT", "Method", "Compound_ID", "HMDB_ID", "Metabolite", "HMDB_ID_Certainty"]); 
                variable_name = "sample"),
            "sample" => ByRow(s-> !startswith(s, "PRE"))
    )
    transform!(df, "sample"=> ByRow(mbx-> begin
        startswith(mbx, "FE") && return mbx
        if startswith(mbx, "SD")
            m = match(r"SD(\d+)([A-Z])", mbx)
            isnothing(m) && error("SD sample $mbx doesn't match pattern")
            sid = lpad(m[1], 4, "0")
            tp = findfirst(c-> "$c" == m[2], "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            mbx = string("C", sid, "_", tp, "E", "_1A")
        end
        mbxid = get(biospdict, mbx, nothing)
        isnothing(mbxid) && return mbx
        return mbxid
    end) => "sample")
end


@chain mbx_long begin
    groupby(["Method", "sample"])
    transform!("value"=> (v-> all(==(0), v) ? fill(1 / length(v), length(v)) : v ./ sum(v)) => "relab")
    transform!("relab"=> (v-> log.(v .+ (minimum(filter(>(0), v)) / 2))) => "logvalue")
end

metabolites = unique(skipmissing(mbx_long.Metabolite))
metabolites[findall(m-> contains(m, r"trypt"i), metabolites)]
metabolites[findall(m-> contains(m, r"Gamma"i), metabolites)]
metabolites[findall(m-> contains(m, r"acet(ate|ic)"i), metabolites)]
metabolites[findall(m-> contains(m, r"glutam(ate|ic)"i), metabolites)]
metabolites[findall(m-> contains(m, r"quino"i), metabolites)]

#- 

using Arrow

Arrow.write("input/metabolites.arrow", mbx_long)

#-

using Resonance
using VKCComputing
using XLSX
using ThreadsX
using Chain
using Arrow

VKCComputing.set_readonly_pat!()

base = LocalBase()

mbx_long = Arrow.Table("input/metabolites.arrow") |> DataFrame
# mbx_long = subset(mbx_long, "sample"=> ByRow(s-> !startswith(s, "M")))

transform!(groupby(mbx_long, "sample"), "sample"=> (scol-> begin
    len = length(scol)
    sample = first(scol)

    # https://vkc-lab.slack.com/archives/C6A0Z4Z8X/p1691697373467609
    if sample == "C1227_3E_1A"
        sample = "FE01922"
    elseif sample == "C1162_3E_1A"
        sample = "FE01852"
    elseif sample == "FE2790"
        sample = "FE02790"
    elseif sample == "M0932_7E_1A"
        sample = "FE01759"
    elseif sample == "C11680_4E_1A"
        sample = "FE01935"
    end

    default_ret = (; subject = fill(missing, len), timepoint = fill(missing, len), biospecimen = fill(missing, len))
    biosp = get(base["Biospecimens"], sample, nothing)
    if isnothing(biosp)
        alias = get(base["Aliases"], sample, nothing)
        (isnothing(alias) || !haskey(alias, :biospecimens)) && return default_ret

        biosp = base[first(alias[:biospecimens])]
    end

    timepoint = get(biosp, :collection, missing)
    subject = get(get(base, first(get(biosp, :subject, ("",))), Dict()), :uid, missing)
    subject = ismissing(subject) ? missing : parse(Int, replace(subject, "resonance-"=>""))
    biospecimen = biosp[:uid]
    
    return (; subject = fill(subject, len), timepoint = fill(timepoint, len), biospecimen = fill(biospecimen, len))
    
end) => ["subject", "timepoint", "biospecimen"])

#-

newmeta = CSV.read("input/ECHO_metabolomics_preselectedinputs.csv", DataFrame)
subtp = Set([(row.subject,row.timepoint) for row in eachrow(newmeta)])
mbx_long = subset(mbx_long, AsTable(["subject", "timepoint"])=> ByRow(nt-> (nt.subject, nt.timepoint) ∈ subtp))

#-

