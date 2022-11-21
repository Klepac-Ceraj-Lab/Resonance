# SRA and DAC metadata for uploads

## SRA

### Get samples not already uploaded

```julia
using Resonance
using Airtable

atmdata = airtable_metadata()
samplemap = Dict(row.sample => (; subject = row.subject, 
                                  timepoint = row.timepoint,
                                  DOC = row.DOC,
                                  childAgeMonths = row.childAgeMonths,
                                  Mother_Child = row.Mother_Child) for row in eachrow(atmdata)
)
                    
mdata = CSV.read(datafiles("wrangled", "timepoints.csv"), DataFrame)

base = AirBase("appSWOVVdqAi5aT5u")
biosamples = Airtable.query(AirTable("NCBI Biosample", base))
recs = Iterators.flatten([get(rec, :Samples, []) for rec in biosamples]) |> collect
samples = filter(rec -> rec.id in recs, Airtable.query(AirTable("Samples", base)))
completed = DataFrame(map(s-> (; id = s.id, sample = get(s, :sample, missing)), samples))
subset!(atmdata, "sample" => ByRow(s-> startswith(s, "FG") && !(s in completed.sample)), 
                 "Mgx_batch" => ByRow(!ismissing)
)


```

### Get files

This step should be run on `hopper`

```julia
files = filter(f-> startswith(f, "FG") && contains(f, r"paired_[12]\.fastq\.gz"), 
              readdir("/grace/sequencing/processed/mgx/kneaddata")
)

filter!(f-> first(split(f, "_")) in atmdata.sample, files)
biosample = DataFrame(
    "Sample Name" => unique(map(f-> first(split(f, "_")), files)),
    "file1"       => unique(f-> first(split(f, "_")), filter(f-> contains(f, "paired_1"), files)),
)
biosample.file1 = map(f-> replace(f, "_kneaddata_paired_1"=> "_L001_R1_001"), biosample.file1)
biosample.file2 = map(f-> replace(f, "_L001_R1_001"=> "_L001_R2_001"), biosample.file1)
biosample.file3 = map(f-> replace(f, "_L001_R1_001"=> "_L002_R1_001"), biosample.file1)
biosample.file4 = map(f-> replace(f, "_L001_R1_001"=> "_L002_R2_001"), biosample.file1)
biosample.file5 = map(f-> replace(f, "_L001_R1_001"=> "_L003_R1_001"), biosample.file1)
biosample.file6 = map(f-> replace(f, "_L001_R1_001"=> "_L003_R2_001"), biosample.file1)
biosample.file7 = map(f-> replace(f, "_L001_R1_001"=> "_L004_R1_001"), biosample.file1)
biosample.file8 = map(f-> replace(f, "_L001_R1_001"=> "_L004_R2_001"), biosample.file1)


transform!(biosample, "Sample Name" => ByRow(s-> samplemap[s].subject) => "subject", 
                      "Sample Name" => ByRow(s-> samplemap[s].timepoint) => "timepoint"
)
```

### Global traits for SRA submission

```julia
biosample[!, "organism"]            .= "human gut metagenome" # txid408170
biosample[!, "host"]                .= "Homo sapiens"
biosample[!, "investigation_type"]  .= "metagenome"
biosample[!, "project_name"]        .= "Human gut metagenomes of healthy mothers and their children, Jan 19 '21"
biosample[!, "lat_lon"]             .= "41.8786 N 71.3831 W" # location of Pawtucket, RI
biosample[!, "geo_loc_name"]        .= "USA: Rhode Island"
biosample[!, "env_biome"]           .= "not aplicable"
biosample[!, "env_feature"]         .= "not aplicable"
biosample[!, "env_material"]        .= "not aplicable"
biosample[!, "env_package"]         .= "human-gut"
biosample[!, "rel_to_oxygen"]       .= "anaerobe"
biosample[!, "samp_collect_device"] .= "Omnigene OMR-200 tube"
biosample[!, "seq_meth"]            .= "illumina NextSeq 550"
biosample[!, "seq_quality_check"]   .= "none"
biosample[!, "host_disease_stat"]   .= "none"
biosample[!, "host_body_product"]   .= "fma64183"
biosample[!, "samp_store_temp"]     .= "-80 degree celcius"
biosample[!, "nucl_acid_ext"]       .= "https://www.qiagen.com/us/resources/resourcedetail?id=84c1f2e7-8db6-4957-a504-92bf9f82dd84"
```

And now for sample-specific attributes:

```julia
# make DataFrame with same rows as SRA
samplemeta = leftjoin(rename(biosample, "Sample Name" => "sample"), unique(mdata, ["subject", "timepoint"]), on=["subject", "timepoint"])

## biosample[!, "adapters"]           .=
## biosample[!, "mid"]                .= {Multiplex ID}
rename!(biosample, "subject" => "subject_id")
transform!(biosample, "Sample Name" => ByRow(s-> samplemap[s].DOC) => "collection_date")
biosample[!, "host_age"]        = map(a-> ismissing(a) ? missing : string.(a ./ 12 .* 365) .* " days", samplemeta.ageMonths)
biosample[!, "host_sex"]        = map(row-> samplemap[row."sample"] == "M" ? "Female" : row.childGender, eachrow(samplemeta))

# Generate output for BioSample Attributes.
rename!(biosample, "subject" => "subject_id")
CSV.write(datafiles("biosample_mgx_attributes.csv"), biosample[!, Not([:file1, :file2, :file3, :file4, :file5, :file6, :file7, :file8])], delim='\t')
```

```julia
sra = select(biosample, ["Sample Name", "file1", "file2", "file3", "file4", "file5", "file6", "file7", "file8"])

sra.title = map(row-> "Shotgun metagenomic sequence of stool sample: subject $(row.subject_id), timepoint $(row.timepoint)", eachrow(biosample))
sra.library_ID = map(s-> "$s-mgx", sra."Sample Name")
sra[!, "library_strategy"] .= "WGS"
sra[!, "library_source"] .= "METAGENOMIC"
sra[!, "library_selection"] .= "RANDOM"
sra[!, "library_layout"] .= "paired"
sra[!, "platform"] .= "ILLUMINA"
sra[!, "instrument_model"] .= "NextSeq 550"
sra[!, "design_description"] .= replace("""
    Libraries were prepared with Illumina Nextera Flex Kit
    for MiSeq and NextSeq from 1 ng of each sample.
    Samples were then pooled onto a plate and sequenced
    on the Illumina NextSeq 550 platform
    using 150+150 bp paired-end “high output” chemistry""",
    '\n'=> " ")
sra[!, "filetype"] .= "fastq"
rename!(sra, "file1"=> "filename", 
             "file2"=> "filename2",
             "file3"=> "filename3",
             "file4"=> "filename4",
             "file5"=> "filename5",
             "file6"=> "filename6",
             "file7"=> "filename7",
             "file8"=> "filename8",
             "Sample Name"=>"sample_name"
)

# Generate output for SRA Attributes.

CSV.write(datafiles("sra_mgx_attributes.csv"), sra[!, :], delim='\t')
open(datafiles("srafiles.txt"), "w") do io
    for row in eachrow(sra)
        println.(io, values(row[["filename", "filename2", "filename3", "filename4", "filename5", "filename6", "filename7", "filename8"]]))
    end
end
```

## ECHO DAC Uploads

Using SOP-EAW-016-Extant-Microbiome protocol.

### A1. Demultiplexed Sequence Data Files

Raw sequence files in fastq format
that were generated from the Illu-mina (or other platform) sequenced
paired-end (and/or sequencing line)
fastq files that have been split
(or “demultiplexed”) by sample

For whole metagenome data:

- Forward lane 1:
  `CohortID_ParticipantID_ SampleNumber_SetID_SequencingType_L001_R1_DMF_YYYY_MMDD.fastq`
- Reverse lane 1:
  `CohortID_ParticipantID_SampleNumber_SetID_SequencingType_L001_R2_DMF_YYYY_MMDD.fastq`
- Forward lane 2:
  `CohortID_ParticipantID_SampleNumber_SetID_SequencingType_L002_R1_DMF_YYYY_MMDD.fastq`
- Reverse lane 2:
  `CohortID_ParticipantID_SampleNumber_SetID_SequencingType_L002_R2_DMF_YYYY_MMDD.fastq`

<!-- Code to rename files -->

### B. Participant-visit level data (1 folder with four required data files)

Submission type: Participant-visit level (four files per cohort stored in a single folder)

There will be 4 types of files (Panel, Batch, Specimen, and Analysis files) 
requested as part (B)“Participant-visit level data” of this transfer of data.
These files contain information about the samples with data being transferred.
**Each row represents one sample**. 
It is possible that an ECHO participant has multiple samples
(e.g., collected and measured at different times or in different tissues).  
The ECHO participant ID should have multiple rows in the file
if the participant has measurements from multiple samples.
The columns of this file will contain information about the specific sample.

Details specifying column names, column order, and a description of the data format for these files are
provided in Tables below. Please note:

- It is critical that the order of the columns and the column names appear exactly as listed
  in Tables B1. Panel File, B2. Batch File, B3. Specimen File, and B4. Analysis file and the data dictionaries, including capitalization.
- Please do not use any spaces or special characters (e.g. &, %, #) for character strings;
  the only allowable special characters are underscores.

The columns are required unless otherwise noted. Please include all column names into the file, even if
data was not collected for that variable.

#### File Nomenclature

Individual file nomenclature:

- `CohortID_SetID_SampleInformation_MicroPanelFile_YYYY_MMDD.csv`
- `CohortID_SetID_SampleInformation_MicroBatchFile_YYYY_MMDD.csv`
- `CohortID_SetID_SampleInformation_MicroSpecimenFile_YYYY_MMDD.csv`
- `CohortID_SetID_SampleInformation_MicroAnalysisFile_YYYY_MMDD.csv`

Folder nomenclature:

Please combine the 4 Table files (Panel, Batch, Specimen, and Analysis file)
within each set into a compressed folder named as follows:

`CohortID_SetID_SampleInformation_YYYY_MMDD.tar.gz`

or

`CohortID_SetID_SampleInformation_YYYY_MMDD.tgz`


```julia
using Resonance
using XLSX
using CSV
using DataFrames
using Airtable
using Dates
using CategoricalArrays

atmdata = airtable_metadata()
extant = DataFrame(XLSX.readtable(datafiles("Extant_Microbiome_092022.xlsx"), "Sheet1", infer_eltypes=true))

subset!(extant, "Participation Level"=> ByRow(p-> parse(Int, p) > 0))

extant.fastq_file_name = map(eachrow(extant)) do row
    join([
        row.cohortid,
        row.participantid,
        row.collectionNum,
        "1", #setID
        "WMS_L001_R1_DMF",
        "2022_1001.fastq.gz"
    ], "_"
    )
end

```

### B1. Panel File

```julia
b1_common = (;
  PanelID = "resonance_upload_1",
  specimen_type = 1, # Stool
  specimen_type_oth_sp = missing,
  panel_lab_name = "Integrated Microbiome Resource",
  panel_lab_country = "CAN", # Canada
  panel_lab_state = "NS",
  panel_lab_city = "Halifax",
  panel_lab_meth = 4, # metagenomic
  panel_lab_meth_oth = missing,
  panel_lab_instr = 9, # other
  panel_lab_instr_oth = "Illumina NextSeq 550",
  panel_data_public = 1, # subset of data

  panel_data_srs = 1,
  panel_data_link = "https://dataview.ncbi.nlm.nih.gov/object/PRJNA695570"
)

b1_panel = copy(extant)
select!(b1_panel,
    "siteid"          => identity => "SiteID",
    "cohortid"        => identity => "CohortID",
    "participantid"   => identity => "ParticipantID",
    "collectionNum"   => identity => "sample_number",
    "fastq_file_name" => identity => "fastq_file_name",
    "SampleID"        => identity => "sample",
    "studyID"         => identity => "subject",
    "collectionNum"   => identity => "timepoint"
)

for key in keys(b1_common)
    b1_panel[!, key] .= b1_common[key]
end


```

### B2 Batch File

```julia
b2_batch = copy(b1_panel)
select!(b2_batch,
    "PanelID",
    "SiteID",
    "CohortID",
    "ParticipantID",
    "sample_number",
    "fastq_file_name",
    "sample",
    "subject",
    "timepoint"
)

atgrp = groupby(atmdata, "sample")

base = AirBase("appSWOVVdqAi5aT5u")
batches = DataFrame([rec.fields for rec in Airtable.query(AirTable("MGX Batches", base))])
batchdict = Dict(row.Name => Date(row."Date Shipped") for row in eachrow(batches))

b2_batch = hcat(b2_batch, DataFrame(map(eachrow(b2_batch)) do row
    misrow = (; batch_num=missing, batch_date_yr=missing, batch_date_mn=missing, batch_date_day=missing)
    key = (; sample = row.sample)
    haskey(atgrp, key) || return misrow
    grp = first(atgrp[key])
    batch_num = grp.Mgx_batch
    ismissing(batch_num) && return misrow
    dt = batchdict[batch_num]
    batch_date_yr = year(dt)
    batch_date_mn = lpad(month(dt), 2, "0")
    batch_date_day = lpad(day(dt), 2, "0")
    
    return (; batch_num, batch_date_yr, batch_date_mn, batch_date_day)
end))

subset!(b2_batch, :batch_num => ByRow(!ismissing))
let keep = Set(b2_batch.sample)
    subset!(b1_panel, "sample"=> ByRow(s-> s in keep))
end
```

### B3 Specimen File

```julia
b3_common = (;
    spec_coll_comp = -6, # Not Applicable
    spec_type      =  1, # Stool
    spec_stl_meth  =  1, # Fecal sample
    spec_coll_meth =  2, # collected by trained staff
    spec_stl_kit   = -6, # Not Applicable
    spec_ship_meth =  5, # Other
    spec_ship_meth_oth = "In box with dry ice",
    spec_coll_prss =  2, # Specimen collected with preservatives
    
)

function lifestage(stage)
    startswith(stage, "P")   && return 2
    startswith(stage, "Inf") && return 3
    startswith(stage, "EC")  && return 4
    startswith(stage, "MC")  && return 5
    startswith(stage, "Ado") && return 6
    return missing
end

b3_specimen = copy(extant)
b3_specimen.Sex = categorical(b3_specimen.Sex; levels = ["Male", "Female"])
select!(b3_specimen,
    "siteid"          => identity => "SiteID",
    "cohortid"        => identity => "CohortID",
    "participantid"   => identity => "ParticipantID",
    "collectionNum"   => identity => "sample_number",
    "fastq_file_name" => identity => "fastq_file_name",
    "SampleID"        => identity => "sample",
    "studyID"         => identity => "subject",
    "collectionNum"   => identity => "timepoint",
    "AgeatCollection_Years" => identity => "spec_age_yrs",
    "AgeOn7012020"    => identity => "spec_age_fixed",
    "echotpcoded"     => ByRow(lifestage) => "spec_coll_ls",
    "Sex"             => ByRow(levelcode) => "spec_sex"
)


for key in keys(b3_common)
    b3_specimen[!, key] .= b3_common[key]
end

let keep = Set(b2_batch.sample)
    subset!(b3_specimen, "sample"=> ByRow(s-> s in keep))
end

```

### B4 Analysis File

```julia
b4_common_mgx = (;
    analysis_dna_extr = 1, # yes
    an_dna_extr_meth  = "extracted",
    an_dna_extr_kit   = 1, # MoBio PowerSoil kit
    analysis_rna_extr = 2, # no
    an_hmgz           = 2, # no
    an_shtgn          = 1, # yes
    an_shtgn_chem_sp  = "Illumina Nextera Flex kit for MiSeq+NextSeq",
    an_shtgn_std_sp   = "None",
    an_amp            = 2
)

b4_common_amp = (;
    an_amp = 1,                                   # Amplicon sequencing conducted  1,Yes | 2,No            
    an_amp_std_sp = "N/A",                        # Standards used
    an_amp_prod = 2,                              # Amplicon products stored       1,Yes | 2,No
    an_amp_chem = "V3",                           # Sequencing chemistry
    an_amp_prim_fwd_seq = "GTGYCAGCMGCCGCGGTAA",  # The forward primer
    an_amp_prim_rev_seq = "CCGYCAATTYMTTTRAGTTT", # The reverse primer
    an_amp_prim_trim = 1,                         # Were the primers trimmed       1,Yes | 2,No
    an_amp_pcrcycle = "",                         # Number of cycles
    an_amp_pairseq = 1,                           # Paired-end sequencing used     1,Yes | 2,No
    an_amp_read = 300,                            # Read length       
)

b4_analysis = copy(b2_batch)
select!(b4_analysis,
    "batch_num" => identity => "BatchID",
    "SiteID",
    "CohortID",
    "ParticipantID",
    "sample_number",
    "fastq_file_name",
    "sample",
    "subject",
    "timepoint",
)


for key in keys(b4_common_mgx)
    b4_analysis[!, key] .= b4_common_mgx[key]
end

```


### Outputs


```julia
CSV.write(datafiles("10701_1_SampleInformation_MicroPanelFile_2022_1001.csv"),
    unique(select(b1_panel, Not(["sample", "subject", "timepoint"])), "fastq_file_name")
)

CSV.write(datafiles("10701_1_SampleInformation_MicroBatchFile_2022_1001.csv"),
    unique(select(b2_batch, Not(["sample", "subject", "timepoint"])), "fastq_file_name")
)

CSV.write(datafiles("10701_1_SampleInformation_MicroSpecimenFile_2022_1001.csv"),
    unique(select(b3_specimen, Not(["sample", "subject", "timepoint"])), "fastq_file_name")
)

CSV.write(datafiles("10701_1_SampleInformation_MicroAnalysisFile_2022_1001.csv"),
    unique(select(b3_specimen, Not(["sample", "subject", "timepoint"])), "fastq_file_name")
)
```


### File copying

```julia

files = [readdir("/lovelace/sequencing/raw/mgx/fastq"; join=true);
         readdir("/lovelace/sequencing/raw/mgx/sra_uploads"; join=true)]
newdir = "/lovelace/sequencing/raw/mgx/echo_dac"

using ThreadsX

for row in eachrow(unique(b2_batch, "fastq_file_name"))
    sample = row.sample
    sfs = filter(f-> occursin(sample, basename(f)), files)
    isempty(sfs) && continue
    length(sfs) == 16 && (sfs = first(sfs, 8))

    @assert length(sfs) == 8
    bits = split(row.fastq_file_name, "_")
    @assert bits[[6,7]] == ["L001", "R1"]
    @info sample

    ThreadsX.foreach(sfs) do sf
        newf = copy(bits)
        lane, fr = match(r"FG\d+_S\d+_(L00[1-4])_(R[12])_001", sf).captures
        newf[[6,7]] = [lane, fr]
        dest = joinpath(newdir, joinpath(newdir, join(newf, "_")))
        !isfile(dest) && cp(sf, dest)
    end
end
```
