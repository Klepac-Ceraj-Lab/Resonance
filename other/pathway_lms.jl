using Resonance
using DataFrames.PrettyTables


mdata = Resonance.load(Metadata())
pathways = Resonance.load_raw_humann(; kind="pathabundance", sample_filter=mdata.sample)
relativeabundance!(pathways)
set!(pathways, mdata)
pdf = comm2wide(pathways)

non_spec_cols = [
    "sample", "subject", "timepoint", "ageMonths", "sex", "race", "education", "date",
    "cogScore", "sample_base", "read_depth", "filter_00to120", "filter_00to06", "filter_18to120"]


plm_00to06 = runlms(pdf[pdf.filter_00to06, :], scratchfiles("pathways"), names(pdf, Not(non_spec_cols)))
sort!(plm_00to06, :qvalue)
pretty_table(select(plm_00to06, "feature", "coef", "pvalue", "qvalue"))


plm_18to120 = runlms(pdf[pdf.filter_18to120, :], scratchfiles("pathways"), names(pdf, Not(non_spec_cols)))
sort!(plm_18to120, :qvalue)
pretty_table(select(plm_18to120, "feature", "coef", "pvalue", "qvalue"))