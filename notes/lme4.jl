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
