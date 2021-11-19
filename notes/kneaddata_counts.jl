using Resonance
using CodecZlib
using DataFramesMeta


if isfile("/home/kevin/Desktop/read_counts.csv")
    df = CSV.read("/home/kevin/Desktop/read_counts.csv", DataFrame)
else
    df = DataFrame()

    kneads = readdir("/grace/echo/analysis/biobakery3/links/kneaddata/", join=true)
    for k in kneads[946:end]
        open(k) do io
            stream = GzipDecompressorStream(io)
            c = count(line-> startswith(line, "@"), eachline(stream))
            push!(df, (file=k, count=c))
        end
    end

    @transform!(df, @byrow :sample = first(split(basename(:file), "_")))
    CSV.write("/home/kevin/Desktop/read_counts.csv", df)
end

df = CSV.read("/home/kevin/Desktop/read_counts.csv", DataFrame)
@transform!(df, @byrow :sample = first(split(basename(:file), "_")))
CSV.write("/home/kevin/Desktop/read_counts.csv", df)
