
```sh
❯ exa -Tl data
drwxrwx---    - kevin 25 Apr 12:39 data
drwxrwxr-x    - kevin 13 Jan 14:58 ├── brain
.rw------- 237k kevin 13 Jan 14:39 │  ├── freesurfer_curvature_leftHemi_oct2021.csv
.rw------- 237k kevin 13 Jan 14:39 │  ├── freesurfer_curvature_rightHemi_oct2021.csv
.rw------- 245k kevin 13 Jan 14:39 │  ├── freesurfer_thickness_leftHemi_oct2021.csv
.rw------- 245k kevin 13 Jan 14:39 │  ├── freesurfer_thickness_rightHemi_oct2021.csv
.rw------- 585k kevin 13 Jan 14:58 │  └── segmentationVolumeMeasurements_oct2021.csv
.rw-rw-r--  63k kevin 22 Nov  2021 ├── echo_id_map.csv
.rw-------  65k kevin  1 Apr 14:06 ├── echoids.csv
.rw-rw-r--  11k kevin 13 Dec  2021 ├── gbm.txt
.rw-------  22k kevin  4 Jan 21:46 ├── genemetab.csv
drwxrwxr-x    - kevin 27 Feb 21:45 ├── lm_results
.rw-------  51k kevin 27 Feb 13:32 │  ├── old_species_predict_cogscore.csv
.rw------- 134k kevin 27 Feb 13:32 │  ├── over12m_cogScores_taxa.csv
.rw------- 154k kevin 27 Feb 13:31 │  ├── species_cogscore.csv
.rw-------  53k kevin 27 Feb 13:32 │  ├── under12m_cogScores_taxa.csv
.rw-------  33k kevin 27 Feb 13:32 │  └── young_species_predict_cogscore.csv
.rw-------  18M kevin 13 Dec  2021 ├── map_ko_uniref90.txt.gz
drwxrwxr-x    - kevin 13 Jan 14:36 ├── metabolites
.rw-rw-r--  25M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_C8-pos_results.xlsx
.rw-r--r--  41M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_C18-neg_results.xlsx
.rw-r--r--  17M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_HILIC-neg_results.xlsx
.rw-r--r--  30M kevin 13 Jan 14:36 │  └── 21_0924_VKC_HILIC-pos_results.xlsx
.rw------- 3.5k kevin 22 Nov  2021 ├── methylation_files.txt
.rw------- 146k kevin 19 Nov  2021 ├── read_counts.csv
.rw------- 4.7k kevin 25 Apr 12:39 ├── README.md
drwxrwxr-x    - kevin  1 Apr 13:17 ├── resonance_fmp
.rw------- 233k kevin  2 Feb 13:11 │  ├── breast.csv
.rw-rw-r--  28k kevin  2 Feb 11:05 │  ├── Child_Medications_012422.xlsx
.rw-rw-r--  43k kevin  1 Apr 12:57 │  ├── COVID_Fecal_040122.xlsx
.rw-------  17k kevin 13 Jan 14:38 │  ├── COVID_Fecal_10252021.xlsx
.rw-rw-r-- 123k kevin  2 Feb 11:06 │  ├── Family_Medical_History_012422.xlsx
.rw-rw-r-- 128k kevin  1 Apr 13:17 │  ├── Fecal_All_033022.xlsx
.rw------- 122k kevin 13 Jan 14:37 │  ├── Sample_Centric_10252021.xlsx
.rw-rw-r-- 507k kevin  2 Feb 11:06 │  ├── Subject_Centric_012522.xlsx
.rw-rw-r-- 444k kevin  1 Apr 12:57 │  ├── Subject_Centric_040122.xlsx
.rw-rw-r-- 1.6M kevin  2 Feb 11:06 │  ├── Timepoint_Centric_012522.xlsx
.rw-rw-r-- 1.9M kevin  7 Feb 13:35 │  ├── Timepoint_Centric_020422.xlsx
.rw-rw-r-- 1.7M kevin  1 Apr 12:57 │  └── Timepoint_Centric_033122.xlsx
.rw------- 4.0M kevin 29 Apr 11:41 ├── species.csv
.rw-r--r--  14k kevin 22 Nov  2021 ├── targets_passed.csv
.rw-rw-r-- 3.8k kevin  1 Oct  2021 ├── VKC.DATA.MI.csv
.rwxrwxrw- 1.3M kevin 15 Jul  2021 ├── VKC_ASVTAB.csv
.rwxrwxrw-  52k kevin 15 Jul  2021 ├── VKC_METATAB.csv
.rwxrwxrw-  81k kevin 15 Jul  2021 ├── VKC_TAXTAB.csv
drwxrwxr-x    - kevin  1 Apr 11:41 └── wrangled
.rw------- 271k kevin  8 Apr 16:43    ├── covid.csv
.rw------- 113k kevin  8 Apr 16:43    ├── etohsamples.csv
.rw------- 301M kevin  7 Feb 16:13    ├── metabolites.csv
.rw------- 138k kevin  8 Apr 16:43    ├── omnisamples.csv
.rw------- 3.8M kevin 13 Jan 15:26    ├── species.csv
.rw-------  81k kevin  8 Apr 16:45    ├── tidy_subjects.csv
.rw-------  31M kevin 23 Feb 19:48    ├── tidy_timepoints.csv
.rw-------  31M kevin  8 Apr 16:45    ├── tidy_timepoints_with_brain.csv
.rw-------  20M kevin  8 Apr 16:42    └── timepoints.csv
```

## Metadata issues identified

1. for `ageLabel`, the only values are "1 and under", "prenatal", or missing.
   Let's just leave that column out of future exports
2. "simple_race" still has mix of `missing`, "Unknown", and "Decline to answer",
   also both "Mixed" and "Mixed Race"
3. `everBreastFed` is totally wrong. Only 18 non-missing values, all 0

## gbm.txt

Supplementary dataset 1 from Valles-Colomer et. al. (2019)

https://doi.org/10.1038/s41564-018-0337-x

https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-018-0337-x/MediaObjects/41564_2018_337_MOESM3_ESM.txt

## Files transferred to Guilherme

- `metadata/` 
  - `covid.csv`: covid-specific sample metadata
  - `etohsamples.csv`: metadata for samples collected in ethanol,
    mostly needed for metabolomics.
    These are the samples that start with `FE`
  - `metabolites.csv`: normalized metabolomics data
  - `omnisamples.csv`: metadata for samples collected in genotek omnigene,
    which is mostly what's used for metagenomics.
    These are the samples that start with `FG`
  - `species.csv`: combined taxonomic profiles (see `notes/basic/species.jl`)
  - `timepoints.csv`: time-point specific metadata for all subjects
  - `volumes.csv`: combined, but not normalized brain volume data
- `batchXXX/output`
  - `humann/`
    - `main/`
      - `XXX_genefamilies.tsv`: stratified gene function profiles (UniRef90)
      - `XXX_pathabundance.tsv`: metabolic pathways relative abundance
      - `XXX_pathcoverage.tsv`: metabolic pathway coverage (see `humann` docs for explanation)
    - `regroup/`: regrouping genefamiles into `ko`, `pfam` etc
    - `rename/`: adding human-readable names to contents of `regroup/`
  - `metaphlan`
    - `XXX_profile.tsv`: taxonomic profile
- `links/`: symlinks to `_genefamilies.tsv` and `_profile.tsv` files for all batches


## Relevant REDCAP coding

- `Redcap_Ess_CHB_IFP::ifp_b01` : ever breastmilk
  - 1: yes
  - 2: no
  - -8: don't know
- `Redcap_Ess_CHB_IFP::ifp_b08_1`: past 7 days mom's breastmilk
  - 1: yes
  - 2: no
  - -8: don't know
- `Redcap_Ess_CHB_IFP::ifp_b09_1`: past 7 days other's breastmilk
  - 1: yes
  - 2: no
  - -8: don't know
