
```sh
❯ ls -al data/
.rw-rw-r-- 9.1M kevin 22 Nov 15:15 21_0924_VKC_C8-pos_results.xlsx
.rw-------  17k kevin 25 Oct 17:09 COVID_Fecal_10252021.xlsx
.rw-rw-r--  63k kevin 22 Nov 11:54 echo_id_map.csv
.rw------- 237k kevin 26 Oct 09:23 freesurfer_curvature_leftHemi_oct2021.csv
.rw------- 237k kevin 26 Oct 09:23 freesurfer_curvature_rightHemi_oct2021.csv
.rw------- 245k kevin 26 Oct 09:23 freesurfer_thickness_leftHemi_oct2021.csv
.rw------- 245k kevin 26 Oct 09:23 freesurfer_thickness_rightHemi_oct2021.csv
.rw------- 679k kevin 22 Nov 15:00 halla_species.tsv
.rw-------  57k kevin 22 Nov 15:00 halla_volumes.tsv
.rw------- 3.5k kevin 22 Nov 11:27 methylation_files.txt
.rw------- 146k kevin 19 Nov 12:38 read_counts.csv
.rw-rw-r-- 1.4k kevin 19 Nov 18:29 README.md
.rw------- 122k kevin 25 Oct 17:09 Sample_Centric_10252021.xlsx
.rw------- 585k kevin 26 Oct 09:23 segmentationVolumeMeasurements_oct2021.csv
.rw------- 227k kevin 25 Oct 17:09 Subject_Centric_10252021.xlsx
.rw-rw-r-- 426k kevin 26 Oct 13:50 Subject_Centric_10262021.xlsx
.rw-r--r--  14k kevin 22 Nov 11:12 targets_passed.csv
.rw------- 903k kevin 25 Oct 17:09 Timepoint_Centric_10252021.xlsx
.rw-rw-r-- 3.8k kevin  1 Oct 14:42 VKC.DATA.MI.csv
.rwxrwxrw- 1.3M kevin 15 Jul 14:23 VKC_ASVTAB.csv
.rwxrwxrw-  52k kevin 15 Jul 14:29 VKC_METATAB.csv
.rwxrwxrw-  81k kevin 15 Jul 14:23 VKC_TAXTAB.csv
.rw------- 110k kevin 19 Nov 13:07 volumes.csv
.rw------- 5.9M kevin 17 Nov 13:10 wrangled.csv
```

## Metadata issues identified

1. for `ageLabel`, the only values are "1 and under", "prenatal", or missing.
   Let's just leave that column out of future exports
2. "simple_race" still has mix of `missing`, "Unknown", and "Decline to answer",
   also both "Mixed" and "Mixed Race"
3. `everBreastFed` is totally wrong. Only 18 non-missing values, all 0
