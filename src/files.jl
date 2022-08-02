analysisfiles(args...) = joinpath(get(ENV, "ANALYSIS_FILES", "./analysis"), args...)
scratchfiles(args...) = joinpath(get(ENV, "SCRATCH_SPACE", "./scratch"), args...)
datafiles(args...) = joinpath(get(ENV, "DATA_FILES", "./data"), args...)

@testset "Dir functions" begin
    @test analysisfiles("test", "thing.txt") == "./analysis/test/thing.txt"    
    @test scratchfiles("test", "thing.txt") == "./data/test/thing.txt"    
    @test datafiles("test", "thing.txt") == "./data/test/thing.txt"    
end