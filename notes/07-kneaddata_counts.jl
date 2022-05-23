using Resonance

kneadread = Resonance.load_knead()

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
