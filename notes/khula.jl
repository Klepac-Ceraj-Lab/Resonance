using Resonance
using Statistics
using CairoMakie

tax = Resonance.load(TaxonomicProfiles())
spec = filter(f-> taxrank(f) == :species, tax)
specdf = DataFrame(spec)
CSV.write("/home/kevin/Desktop/taxa.csv", specdf)

ur = Resonance.load(UnirefProfiles())
ur = filter(f-> !hastaxon(f), ur)

ecs = Resonance.load(ECProfiles())
ecs = filter(f-> !hastaxon(f), ecs)
ecdf = DataFrame(ecs)
CSV.write("/home/kevin/Desktop/ecs.csv", ecdf)


let prev = prevalence(spec)
    @info "species mean prevalence:" mean(prev)
    @info "species median prevalence:" median(prev)
    @info "species minimum prevalence:" minimum(prev)
end

let prev = prevalence(ur)
    @info "uniref mean prevalence:" mean(prev)
    @info "uniref median prevalence:" median(prev)
    @info "uniref minimum prevalence:" minimum(prev)
end

let prev = prevalence(ecs)
    @info "ecs mean prevalence:" mean(prev)
    @info "ecs median prevalence:" median(prev)
    @info "ecs minimum prevalence:" minimum(prev)
end

1 .- vec(abundances(ecs[r"ungrouped"i, :]) ./ sum(abundances(ecs); dims=1)) |> extrema
1 .- vec(abundances(ecs[r"ungrouped"i, :]) ./ sum(abundances(ecs); dims=1)) |> median

hist(vec(prevalence(spec)) .* 100; axis=(; xlabel="prevalence (%)", ylabel = "number of species"), bins=20)
hist(vec(prevalence(ecs)) .* 100; axis=(; xlabel="prevalence (%)", ylabel = "number of ECs"), bins=20)
hist(vec(prevalence(ur)) .* 100; axis=(; xlabel="prevalence (%)", ylabel = "number of UniRefs"), bins=20)

mdf = DataFrame(metadata(spec))

hist(filter(!=(0), vec(prevalence(spec[:, mdf.ageMonths .>= 12]))) .* 100; axis=(; xlabel="prevalence (%)", ylabel = "number of species"), bins=20)
hist!(filter(!=(0), vec(prevalence(spec[:, mdf.ageMonths .< 12]))) .* 100; bins=20)
current_figure()