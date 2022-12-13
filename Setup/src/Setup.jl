module Setup

using Comonicon
using Downloads
using Random

include("download_locations.jl")

struct Taxa end
struct Unirefs end
struct ECs end
struct KOs end
struct Neuro end
struct Timepoints end

function _downloader(id, dir)
    mkpath(dir)
    tarpath = joinpath(dir, "$(randstring()).tar.gz")
    Downloads.download("https://osf.io/download/" * id, tarpath)
    
    run(Cmd(["tar", "xvzf", tarpath, "--directory", dir]))
    rm(tarpath; force=true)
    return nothing
end

function datadownload(::Taxa; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["taxa.arrow", "taxa_features.txt", "samples.txt"]))
        @info "Downloading taxonomic profiles!"
        _downloader(taxa_osfio, inputdir)
    else
        @info "Taxonomic profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::Unirefs; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["genefamilies.arrow", "genefamilies_features.txt", "samples.txt"]))
        @info "Downloading Uniref profiles!"
        _downloader(genefamilies_osfio, inputdir)
    else
        @info "Uniref profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::KOs; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["kos.arrow", "kos_features.txt", "samples.txt"]))
        @info "Downloading KO profiles!"
        _downloader(kos_osfio, inputdir)
    else
        @info "KO profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::ECs; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["ecs.arrow", "ecs_features.txt", "samples.txt"]))
        @info "Downloading EC profiles!"
        _downloader(ecs_osfio, inputdir)
    else
        @info "EC profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::Neuro; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["brain_normalized.csv"]))
        @info "Downloading normalized brain data!"
        _downloader(brain_osfio, inputdir)
    else
        @info "Normalized brain data is already present, skipping!"
    end
    return nothing
end

function datadownload(::Timepoints; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["timepoint_metadata.csv"]))
        @info "Downloading Timepoint metadata!"
        _downloader(timepoints_osfio, inputdir)
    else
        @info "Timepoint metadata is already present, skipping!"
    end
    return nothing
end



# input
# scratch
# models
# figures
# tables
end # module
