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

##

human = DataFrame() 

for (root, dirs, files) in walkdir("/grace/echo/analysis/biobakery3/")
    contains(root, "kneaddata") || continue
    for f in files
        (contains(f, "hg37") || contains(f, "sapiens")) || continue
        s = first(split(f, "_"))
        try
            open(joinpath(root, f)) do io
                stream = GzipDecompressorStream(io)
                c = count(line-> startswith(line, "@"), eachline(stream))
                push!(human, (sample=s, count=c))
            end
        catch e
            @error e
        end
    end
end

##

using CairoMakie

human = @chain human begin
    groupby(:sample)
    combine(:count=> sum)    
end

hist(human.count_sum .* 150)
using Statistics

median(human.count_sum .* 150)
mean(human.count_sum .* 150)
minimum(human.count_sum .* 150)
quantile(human.count_sum .* 150, [0.10, 0.9])
