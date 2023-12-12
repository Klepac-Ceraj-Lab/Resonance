using Resonance
using DataFrames
using CSV
using Chain

subject_centric_path = "/brewster/kevin/resonance_data/resonance_fmp/Participant_Centric_081023.xlsx"


request = DataFrame(id = string.(getindex.(split.(readlines("data/oif_cohort_ids.txt"), "."), 2)))
fmp_subject = DataFrame(XLSX.readtable(subject_centric_path, "Sheet1", infer_eltypes=true))

leftjoin!(request, select(fmp_subject, "ECHOProtocolDBioMom"=> "id", Cols(:)); on="id")
transform!(request, "ECHOProtocolDChild"=> ByRow(cid-> ismissing(cid) ? missing : last(split(cid, "-")))=> "child_id")

CSV.write("/home/kevin/Downloads/tiange_request_ids.csv", select(request, "studyID"=> "subject", "child_id", "childGender", "birthType", r"BeforePregnancy", r"Pregnancy::antibiotics", r"Prenatal"))

#-

base = LocalBase()

sample_metadata = DataFrame(mapreduce(vcat, base["Subjects"][map(s-> string("resonance-", lpad(s, 4, '0')), skipmissing(request.studyID))]) do rec
    biosp = base[rec[:Biospecimens]]
    subject = replace(rec[:uid], "resonance-"=> "")
    seqprep = mapreduce(vcat, biosp) do b
        haskey(b, :seqprep) || return []
        biospecimen = b[:uid]
        collection = get(b, :collection, missing)
    
        [(; subject, biospecimen, collection, seqprep=seq[:uid]) for seq in base[b[:seqprep]]]
    end
    seqprep
end)

timepoint_centric_path = "/brewster/kevin/resonance_data/resonance_fmp/Timepoint_Centric_081023.xlsx"
fmp_timepoint = DataFrame(XLSX.readtable(timepoint_centric_path, "Sheet1", infer_eltypes=true))
transform!(fmp_timepoint,
    "studyID"=> ByRow(s-> lpad(s, 4, '0')) => "subject", 
    AsTable(["assessmentAgeDays", "assessmentAgeMonths", "scanAgeDays", "scanAgeMonths"]) => ByRow(row-> begin
    aad, aam, sad, sam = row
    calcage = coalesce(aad / 365 * 12 + aam, sad / 365 * 12 + sam)

    return calcage
    end) => "ageMonths"
)

sample_metadata = leftjoin(subset(sample_metadata, "subject"=> ByRow(!ismissing), "collection"=> ByRow(!ismissing)), 
          select(subset(fmp_timepoint, "subject"=> ByRow(!ismissing), "timepoint"=> ByRow(!ismissing)),"subject", "timepoint"=> "collection", "ECHOTPCoded", "ageMonths");
          on=["subject", "collection"]
)

sample_metadata
CSV.write("/home/kevin/Downloads/tiange_request_samples.csv", sample_metadata)

CSV.write("/home/kevin/Downloads/need_ga.csv",
    select(sort(subset(sample_metadata, "ECHOTPCoded"=> ByRow(tp-> startswith(tp, "Pre"))), ["subject", "collection"]), "subject", "collection"=>"timepoint", "ECHOTPCoded")
)

