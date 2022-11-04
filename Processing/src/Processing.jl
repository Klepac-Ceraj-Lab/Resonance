module Processing

using Comonicon
using Downloads

include("download_locations.jl")

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
    
    if overwrite || any(!isfile, joinpath.(inputdir, ["taxa.arrow", "taxa_features.txt"]))
        mkpath(inputdir)
        tarpath = joinpath(inputdir, "taxa.tar.gz")
        Downloads.download("https://drive.google.com/uc?export=download&id=" * taxa_gdrive, tarpath)

        run(Cmd(["tar", "xvzf", tarpath, "--directory", inputdir]))
        rm(tarpath; force=true)
    else
        @info "Taxonomic profiles are already present, skipping!"
    end
    return
end


# input
# scratch
# models
# figures
# tables
end # module
