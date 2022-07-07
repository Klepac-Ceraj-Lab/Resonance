
```sh
❯ pwd
/brewster/kevin/resonance_data

❯ exa -lTa ./
drwxrwxr-x    - kevin  7 Jul 13:11 .
.rw-rw-r-- 386k kevin  6 Jun 11:42 ├── 85_snps_20220401.csv
.rw-rw-r-- 1.1k kevin  4 May 14:27 ├── aditi_samples.txt
drwxrwxr-x    - kevin 13 Jan 14:58 ├── brain
.rw-r--r-- 237k kevin 13 Jan 14:39 │  ├── freesurfer_curvature_leftHemi_oct2021.csv
.rw-r--r-- 237k kevin 13 Jan 14:39 │  ├── freesurfer_curvature_rightHemi_oct2021.csv
.rw-r--r-- 245k kevin 13 Jan 14:39 │  ├── freesurfer_thickness_leftHemi_oct2021.csv
.rw-r--r-- 245k kevin 13 Jan 14:39 │  ├── freesurfer_thickness_rightHemi_oct2021.csv
.rw-r--r-- 585k kevin 13 Jan 14:58 │  └── segmentationVolumeMeasurements_oct2021.csv
.rw-rw-r--  63k kevin 22 Nov  2021 ├── echo_id_map.csv
.rw-r--r--  65k kevin  1 Apr 14:06 ├── echoids.csv
.rw-rw-r-- 3.8k kevin 29 Jun 13:24 ├── fsea_all.csv
.rw-rw-r--  11k kevin 13 Dec  2021 ├── gbm.txt
.rw-r--r--  22k kevin  4 Jan 21:46 ├── genemetab.csv
.rw-r--r--  18M kevin 13 Dec  2021 ├── map_ko_uniref90.txt.gz
drwxrwxr-x    - kevin 13 Jan 14:36 ├── metabolites
.rw-rw-r--  25M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_C8-pos_results.xlsx
.rw-r--r--  41M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_C18-neg_results.xlsx
.rw-r--r--  17M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_HILIC-neg_results.xlsx
.rw-r--r--  30M kevin 13 Jan 14:36 │  └── 21_0924_VKC_HILIC-pos_results.xlsx
.rw-r--r-- 3.5k kevin 22 Nov  2021 ├── methylation_files.txt
.rw-r--r-- 258k kevin  2 Jun 15:28 ├── read_counts.csv
drwxrwxr-x    - kevin  6 Jul 23:20 ├── resonance_fmp
.rw-r--r-- 233k kevin  2 Feb 13:11 │  ├── breast.csv
.rw-rw-r--  28k kevin  2 Feb 11:05 │  ├── Child_Medications_012422.xlsx
.rw-rw-r--  43k kevin  1 Apr 12:57 │  ├── COVID_Fecal_040122.xlsx
.rw-r--r--  17k kevin 13 Jan 14:38 │  ├── COVID_Fecal_10252021.xlsx
.rw-rw-r-- 123k kevin  2 Feb 11:06 │  ├── Family_Medical_History_012422.xlsx
.rw-rw-r-- 128k kevin  1 Apr 13:17 │  ├── Fecal_All_033022.xlsx
.rw-r--r-- 122k kevin 13 Jan 14:37 │  ├── Sample_Centric_10252021.xlsx
.rw-rw-r-- 507k kevin  2 Feb 11:06 │  ├── Subject_Centric_012522.xlsx
.rw-rw-r-- 444k kevin  1 Apr 12:57 │  ├── Subject_Centric_040122.xlsx
.rw-rw-r-- 523k kevin  2 Jun 13:11 │  ├── Subject_Centric_060222.xlsx
.rw-rw-r-- 1.6M kevin  2 Feb 11:06 │  ├── Timepoint_Centric_012522.xlsx
.rw-rw-r-- 1.9M kevin  7 Feb 13:35 │  ├── Timepoint_Centric_020422.xlsx
.rw-rw-r-- 1.7M kevin  1 Apr 12:57 │  └── Timepoint_Centric_033122.xlsx
.rw-r--r-- 214k kevin  6 Jun 15:02 ├── snps_cog.csv
.rw-r--r-- 4.0M kevin  2 Jun 14:03 ├── species.csv
.rw-r--r--  14k kevin 22 Nov  2021 ├── targets_passed.csv
.rw-rw-r-- 3.8k kevin  1 Oct  2021 ├── VKC.DATA.MI.csv
.rwxrwxrwx 1.3M kevin 15 Jul  2021 ├── VKC_ASVTAB.csv
.rwxrwxrwx  52k kevin 15 Jul  2021 ├── VKC_METATAB.csv
.rwxrwxrwx  81k kevin 15 Jul  2021 ├── VKC_TAXTAB.csv
drwxrwxr-x    - kevin  1 Apr 11:41 └── wrangled
.rw-r--r-- 272k kevin  2 Jun 14:33    ├── covid.csv
.rw-r--r-- 113k kevin  2 Jun 14:33    ├── etohsamples.csv
.rw-r--r-- 301M kevin  7 Feb 16:13    ├── metabolites.csv
.rw-r--r-- 138k kevin  2 Jun 14:33    ├── omnisamples.csv
.rw-r--r-- 3.8M kevin 13 Jan 15:26    ├── species.csv
.rw-r--r--  81k kevin  8 Apr 16:45    ├── tidy_subjects.csv
.rw-r--r--  31M kevin 23 Feb 19:48    ├── tidy_timepoints.csv
.rw-r--r--  31M kevin  8 Apr 16:45    ├── tidy_timepoints_with_brain.csv
.rw-r--r--  20M kevin  2 Jun 14:33    └── timepoints.csv
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
