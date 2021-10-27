_letter2number(l) = findfirst(==(only(l)), "abcdefghijklmnopqrstuvwxyz")


function split_tp(id)
    m = match(r"(^\d+)(\w)?$", id)
    isnothing(m) && throw(ArgumentError("`$id` does not match expected pattern: `^\\d+\\w?\$`"))

    sub = parse(Int, m[1])
    tp = isnothing(m[2]) ? 1 : _letter2number(m[2])
    return (subject=sub, timepoint=tp)
end

function brain_ingest(file; label=nothing)
    df = CSV.read(file, DataFrame)
    df = hcat(df, DataFrame(split_tp.(df[!, 1])))

    !isnothing(label) && (df[:, "has_$label"] = fill(true, nrow(df)))
    return df
end