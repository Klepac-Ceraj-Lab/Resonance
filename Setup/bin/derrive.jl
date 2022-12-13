using Resonance
using Chain
using CairoMakie
using AlgebraOfGraphics

mdata = Resonance.load_raw_metadata()
taxa = Resonance.load_raw_metaphlan()
ecs = Resonance.load_raw_humann(; kind = "ecs", names=true, stratified=false)
kos = Resonance.load_raw_humann(; kind = "kos", names=true, stratified=false)
# unirefs = Resonance.load_raw_humann(; kind = "genefamilies", names=false, stratified=false)

#- mgx / mbx

let
    df = @chain mdata begin
        groupby([:subject])
        combine(:omni => (col-> count(!ismissing, col)) => "n_mgx",
                :etoh => (col-> count(!ismissing, col)) => "n_mbx")
        subset("n_mgx"=> ByRow(>(0)))
        groupby(["n_mgx", "n_mbx"])
        combine("n_mgx" => length => "count")
    end
    
    plt = data(df) * mapping(:n_mbx=> "N mbx samples", :count; layout = :n_mgx => nonnumeric) * visual(BarPlot)
    draw(plt; figure=(;title = "number of samples"))
end

#- visits vs samples

hasmgx = @chain subset(mdata, :omni=> ByRow(!ismissing)) begin
    groupby([:subject])
    combine(:omni => length => :n)
end
hasmgx.kind .= "hasmgx"
hasvisit = @chain mdata begin
    groupby([:subject])
    combine(:omni => length => :n)
end
hasvisit.kind .= "hasvisit"
df = vcat(hasmgx, hasvisit)

#-

fig = Figure()
visits = Axis(fig[1,1]; title  = "Visits", xticks=(1.5:12.5, string.(1:12)),
                        xlabel = "number of visits",
                        ylabel = "number of subjects"
)
samples = Axis(fig[1,2]; title  = "Stool samples", xticks=(1.5:12.5, string.(1:12)),
                         xlabel = "number of stool samples",
                         ylabel = "number of subjects"
)
subs = Axis(fig[2,1]; title = "All visits",
                      xlabel = "Age (years)",
                      ylabel = "Count"
)

subws = Axis(fig[2,2]; title = "Visits with stool",
                       xlabel = "Age (years)",
                       ylabel = "Count"
)


hist!(visits, hasvisit.n; bins = length(unique(hasvisit.n)))
hist!(samples, hasmgx.n; bins = length(unique(hasmgx.n)))

linkaxes!(visits, samples)

hist!(subs, collect(skipmissing(mdata.ageMonths)) ./ 12)
hist!(subws, collect(skipmissing(subset(mdata, "omni"=> ByRow(!ismissing)).ageMonths)) ./ 12)


fig