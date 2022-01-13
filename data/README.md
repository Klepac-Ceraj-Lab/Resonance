
```sh
❯ ❯ exa -Tl data
drwxrwx---    - kevin 13 Jan 15:21 data
drwxrwxr-x    - kevin 13 Jan 14:58 ├── brain
.rw------- 237k kevin 13 Jan 14:39 │  ├── freesurfer_curvature_leftHemi_oct2021.csv
.rw------- 237k kevin 13 Jan 14:39 │  ├── freesurfer_curvature_rightHemi_oct2021.csv
.rw------- 245k kevin 13 Jan 14:39 │  ├── freesurfer_thickness_leftHemi_oct2021.csv
.rw------- 245k kevin 13 Jan 14:39 │  ├── freesurfer_thickness_rightHemi_oct2021.csv
.rw------- 585k kevin 13 Jan 14:58 │  └── segmentationVolumeMeasurements_oct2021.csv
.rw-rw-r--  63k kevin 22 Nov  2021 ├── echo_id_map.csv
.rw-rw-r--  11k kevin 13 Dec  2021 ├── gbm.txt
.rw-------  22k kevin  4 Jan 21:46 ├── genemetab.csv
.rw-------  18M kevin 13 Dec  2021 ├── map_ko_uniref90.txt.gz
drwxrwxr-x    - kevin 13 Jan 14:36 ├── metabolites
.rw-rw-r--  25M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_C8-pos_results.xlsx
.rw-r--r--  41M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_C18-neg_results.xlsx
.rw-r--r--  17M kevin 13 Jan 14:36 │  ├── 21_0924_VKC_HILIC-neg_results.xlsx
.rw-r--r--  30M kevin 13 Jan 14:36 │  └── 21_0924_VKC_HILIC-pos_results.xlsx
.rw------- 301M kevin 12 Dec  2021 ├── metabolites.csv
.rw------- 3.5k kevin 22 Nov  2021 ├── methylation_files.txt
.rw------- 146k kevin 19 Nov  2021 ├── read_counts.csv
.rw-rw-r-- 2.7k kevin 13 Dec  2021 ├── README.md
drwxrwxr-x    - kevin 13 Jan 14:58 ├── resonance_fmp
.rw-------  17k kevin 13 Jan 14:38 │  ├── COVID_Fecal_10252021.xlsx
.rw------- 122k kevin 13 Jan 14:37 │  ├── Sample_Centric_10252021.xlsx
.rw-rw-r-- 426k kevin 13 Jan 14:37 │  ├── Subject_Centric_10262021.xlsx
.rw-rw-r-- 1.1M kevin 13 Jan 14:37 │  └── Timepoint_Centric_122321.xlsx
.rw------- 3.8M kevin 13 Dec  2021 ├── species.csv
.rw-r--r--  14k kevin 22 Nov  2021 ├── targets_passed.csv
.rw-rw-r-- 3.8k kevin  1 Oct  2021 ├── VKC.DATA.MI.csv
.rwxrwxrw- 1.3M kevin 15 Jul  2021 ├── VKC_ASVTAB.csv
.rwxrwxrw-  52k kevin 15 Jul  2021 ├── VKC_METATAB.csv
.rwxrwxrw-  81k kevin 15 Jul  2021 ├── VKC_TAXTAB.csv
.rw------- 110k kevin 19 Nov  2021 ├── volumes.csv
.rw-------  13M kevin 13 Jan 15:13 └── wrangled.csv
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