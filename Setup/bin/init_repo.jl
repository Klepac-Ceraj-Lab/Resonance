using Comonicon
using Setup
import Setup: datadownload, Taxa, Timepoints, Neuro, ECs, Unirefs, KOs

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
