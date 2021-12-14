struct Metabolite <: Microbiome.AbstractFeature
    name::String
    commonname::Union{Missing, String}
    mz::Float64
    rt::Float64
end

name(m::Metabolite) = m.name
commonname(m::Metabolite) = m.commonname
masscharge(m::Metabolite) = m.mz
retentiontime(m::Metabolite) = m.rt

function load_metabolites(file; sheet=nothing, header=11)
    df = DataFrame(XLSX.readtable(file, isnothing(sheet) ? 1 : sheet, first_row=header, infer_eltypes=true)...)
end

function pull_row(df, name)
    row = findfirst(x-> !ismissing(x) && x == name, df[!, :Metabolite])
    row = vec(permutedims(Matrix(df[row:row, r"^FE\d+"])))
    [ismissing(x) ? 0 : x for x in row]
end