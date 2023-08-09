using Resonance
using XLSX
using ThreadsX

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
