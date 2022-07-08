analysisfiles(args...) = joinpath(ENV["ANALYSIS_FILES"], args...)
scratchfiles(args...) = joinpath(ENV["SCRATCH_SPACE"], args...)
datafiles(args...) = joinpath(ENV["DATA_FILES"], args...)