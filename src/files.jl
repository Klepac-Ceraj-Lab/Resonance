inputfiles(args...) = joinpath(get(ENV, "INPUT_FILES", "./input"), args...)
analysisfiles(args...) = joinpath(get(ENV, "ANALYSIS_FILES", "./analysis"), args...)
scratchfiles(args...) = joinpath(get(ENV, "SCRATCH_FILES", "./scratch"), args...)
datafiles(args...) = joinpath(get(ENV, "DATA_FILES", "./data"), args...)
outputfiles(args...) = joinpath(get(ENV, "OUTPUT_FILES", "./output"), args...)
figurefiles(args...) = joinpath(get(ENV, "FIGURE_FILES", "./figures"), args...)
tablefiles(args...) = joinpath(get(ENV, "TABLE_FILES", "./tables"), args...)

@testset "Dir functions" begin
    withenv("ANALYSIS_FILES"=> nothing, "SCRATCH_FILES"=> nothing, "DATA_FILES" => nothing) do
        @test analysisfiles("test", "thing.txt") == "./analysis/test/thing.txt"
        @test scratchfiles("test", "thing.txt") == "./scratch/test/thing.txt"
        @test datafiles("test", "thing.txt") == "./data/test/thing.txt"
    end
end