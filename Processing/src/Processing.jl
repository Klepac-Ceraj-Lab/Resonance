module Processing

using Comonicon
using Downloads

include("download_locations.jl")

struct Taxa end
struct Unirefs end
struct ECs end
struct KOs end
struct Neuro end
struct Timepoints end

function _downloader(id, dir)
    mkpath(dir)
    tarpath = joinpath(dir, "tmp.tar.gz")
    Downloads.download("https://drive.google.com/uc?export=download&id=" * id, tarpath)
    
    run(Cmd(["tar", "xvzf", tarpath, "--directory", dir]))
    rm(tarpath; force=true)
    return nothing
end

function datadownload(::Taxa; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["taxa.arrow", "taxa_features.txt", "samples.txt"]))
        _downloader(taxa_gdrive, inputdir)
    else
        @info "Taxonomic profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::Unirefs; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["genefamilies.arrow", "genefamilies_features.txt", "samples.txt"]))
        _downloader(genefamilies_gdrive, inputdir)
    else
        @info "Uniref profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::KOs; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["kos.arrow", "kos_features.txt", "samples.txt"]))
        _downloader(kos_gdrive, inputdir)
    else
        @info "KO profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::ECs; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["ecs.arrow", "ecs_features.txt", "samples.txt"]))
        _downloader(ecs_gdrive, inputdir)
    else
        @info "EC profiles are already present, skipping!"
    end
    return nothing
end

function datadownload(::Neuro; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["brain_normalized.csv"]))
        _downloader(brain_gdrive, inputdir)
    else
        @info "Normalized brain data is already present, skipping!"
    end
    return nothing
end

function datadownload(::Timepoints; basedir=pwd(), overwrite=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    if overwrite || any(!isfile, joinpath.(inputdir, ["timepoint_metadata.csv"]))
        _downloader(timepoints_gdrive, inputdir)
    else
        @info "Timepoint metadata is already present, skipping!"
    end
    return nothing
end


"""
Set up your environment for use with analysis scripts
and machine learning models.

# Options

- `-b, --basedir`: Directory to use as default if other options aren't set.
- `--inputdir`: Directory for downloading files to use as inputs.

# Flags

- `--overwrite`: Overwrite existing files. Default is to skip existing files.

"""
@main function setup(; basedir=pwd(), overwrite::Bool=false, inputdir=get(ENV, "INPUT_FILES", joinpath(basedir, "input")))
    datadownload(Taxa(); basedir, overwrite, inputdir)
    datadownload(Timepoints(); basedir, overwrite, inputdir)
    datadownload(Neuro(); basedir, overwrite, inputdir)
    datadownload(KOs(); basedir, overwrite, inputdir)
    datadownload(ECs(); basedir, overwrite, inputdir)
    datadownload(Unirefs(); basedir, overwrite, inputdir)
    
    return nothing
end


# input
# scratch
# models
# figures
# tables
end # module
