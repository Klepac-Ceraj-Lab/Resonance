
function _read_knead_read_counts(dir=ENV["ANALYSIS_FILES"])
    fileinfo = NamedTuple[]

    paired_paths = Iterators.filter(walkpath(PosixPath(dir))) do file 
        "kneaddata" in file.segments || return false
        contains(filename(file), r"read_counts") || return false
        return true
    end

    df = DataFrame()

    for path in paired_paths
        @info path
        batch = path.segments[findfirst(s-> contains(s, "batch"), path.segments)]
        logs = CSV.read(path, DataFrame; stringtype=String, missingstring=["", "NA"])
        rename!(col-> replace(col, "hg37dec_v0.1"=>"Homo_sapiens"), logs)
        subset!(logs, "raw pair1"=>ByRow(!ismissing))
        foreach(col-> logs[!, col] = coalesce.(logs[!, col], 0), names(logs)[2:end])
        logs[!, :batch] .= batch
        append!(df, logs)
    end

    CSV.write("data/read_counts.csv", df)
end

##

function load_knead(dir=ENV["ANALYSIS_FILES"])
    !isfile("data/read_counts.csv") && _read_knead_read_counts(dir)
    CSV.read("data/read_counts.csv", DataFrame)
end