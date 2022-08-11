using CSV
using DataFrames
using DataFrames.PrettyTables
using GLM
using RCall

dat = CSV.read("/home/kevin/Downloads/dat.csv", DataFrame)
dat=df
modjl = lm(@formula(y ~ x), dat)
modjldf = DataFrame(coeftable(modjl))
pretty_table(select(modjldf, 1:5))

R"library('lme4')"
@rput dat
R"modR <- lm(y ~ x, dat)"
R"summary(modR)"


modjl2 = lm(@formula(y ~ x + c1), dat)
modjl2df = DataFrame(coeftable(modjl2))
pretty_table(select(modjl2df, 1:5))

R"modR2 <- lm(y ~ x + c1, dat)"
R"summary(modR2)"

dat.c3 = dat.c1 ./ 1e6
@rput dat

modjl3 = lm(@formula(y ~ x + c3), dat)
modjl3df = DataFrame(coeftable(modjl3))
pretty_table(select(modjl3df, 1:5))

R"modR3 <- lm(y ~ x + c3, dat)"
R"summary(modR3)"


using Distributions

df = DataFrame(y = rand(Normal(), 200), x = rand(Normal(), 200), c1 = rand(Normal(1e7, 1e6), 200))
df.c2 = df.c1 ./ 1e6

modjl = DataFrame(coeftable(lm(@formula(y ~ x + c1), df; dropcollinear=false)))
modjl2 = DataFrame(coeftable(lm(@formula(y ~ x + c2), df)))

@rput df
R"modR <- lm(y ~ x + c1, df)"
R"modR2 <- lm(y ~ x + c2, df)"

pretty_table(select(modjl, 1:5))
pretty_table(select(modjl2, 1:5))
R"summary(modR)"
R"summary(modR2)"