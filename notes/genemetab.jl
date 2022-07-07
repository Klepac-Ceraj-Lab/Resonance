using GLMakie
using CSV
using DataFrames

genemetab = CSV.read(datafiles("genemetab.csv"), DataFrame)

scatter(log.(1 .+ genemetab.glutdegr), log.(1 .+ genemetab.glutsynth), genemetab.glutgut,
    axis=(; xlabel="degradataion", ylabel="synthesis", zlabel="concentration"))