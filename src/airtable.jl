"""
    airtable_metadata(key=ENV["AIRTABLE_KEY"])
Get fecal sample metadata table from airtable.
The API `key` comes from https://airtable.com/account.
This is unlikely to work if you're not in the VKC lab,
but published sample metadata is available from OSF.io
using `datadep"sample metadata"`.
"""
function airtable_metadata(key=Airtable.Credential())
    records = []
    req = Airtable.get(key, "/v0/appyRaPsZ5RsY4A1h", "Master"; view="ALL_NO_EDIT")
    append!(records, req.records)
    while haskey(req, :offset)
        @info "Making another request"
        req = Airtable.get(key, "/v0/appyRaPsZ5RsY4A1h/", "Master"; view="ALL_NO_EDIT", offset=req.offset)
        append!(records, req.records)
        sleep(0.250)
    end

    df = DataFrame()
    for record in records 
        append!(df, filter(p -> !(last(p) isa AbstractArray), record.fields), cols=:union)
    end

    rename!(df, "TimePoint"=>"timepoint", "SubjectID"=>"subject")
    subset!(df, AsTable(["subject", "timepoint"])=> ByRow(r-> !any(ismissing, r)))

    transform!(df, "subject"   => ByRow(s-> parse(Int, s)) => "subject",
                   "timepoint" => ByRow(tp-> parse(Int, tp)) => "timepoint",
                   "Mgx_batch" => ByRow(b-> (!ismissing(b) && contains(b, r"Batch (\d+)")) ? parse(Int, match(r"Batch (\d+)", b).captures[1]) : missing) => "Mgx_batch",
                   "16S_batch" => ByRow(b-> (!ismissing(b) && contains(b, r"Batch (\d+)")) ? parse(Int, match(r"Batch (\d+)", b).captures[1]) : missing) => "16S_batch",
                   "Metabolomics_batch" => ByRow(b-> (!ismissing(b) && contains(b, r"Batch (\d+)")) ? parse(Int, match(r"Batch (\d+)", b).captures[1]) : missing) => "Metabolomics_batch")
    return select(df, Cols(:sample, :subject, :timepoint, :))
end