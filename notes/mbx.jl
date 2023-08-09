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

#-

fixcol(col, n::String)                  = fixcol(col, Val{Symbol(n)}())
fixcol(col, n::Val)                     = coalesce.(col, 0)
fixcol(col, n::Val{:MZ})                = convert(Vector{Float64}, col)
fixcol(col, n::Val{:RT})                = convert(Vector{Float64}, col)
fixcol(col, n::Val{:Method})            = convert(Vector{String}, col)
fixcol(col, n::Val{:Compound_ID})       = convert(Vector{String}, col)
fixcol(col, n::Val{:HMDB_ID})           = convert(Vector{Union{Missing,String}}, col)
fixcol(col, n::Val{:Metabolite})        = convert(Vector{Union{Missing,String}}, col)
fixcol(col, n::Val{:HMDB_ID_Certainty}) = map(col) do c
    c isa Int && return c
    ismissing(c) && return c
    c == "NA" && return 0
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

append!(mbx_long, let file = mbx_xlsx_files[5]
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
end)

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