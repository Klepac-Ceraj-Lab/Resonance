
```sh
$ ❯ ls -al data/
.rw-------  17k kevin 25 Oct 17:09 COVID_Fecal_10252021.xlsx
.rw------- 237k kevin 26 Oct 09:23 freesurfer_curvature_leftHemi_oct2021.csv
.rw------- 237k kevin 26 Oct 09:23 freesurfer_curvature_rightHemi_oct2021.csv
.rw------- 245k kevin 26 Oct 09:23 freesurfer_thickness_leftHemi_oct2021.csv
.rw------- 245k kevin 26 Oct 09:23 freesurfer_thickness_rightHemi_oct2021.csv
.rw-rw-r--  327 kevin 25 Oct 17:10 README.md
.rw------- 122k kevin 25 Oct 17:09 Sample_Centric_10252021.xlsx
.rw------- 585k kevin 26 Oct 09:23 segmentationVolumeMeasurements_oct2021.csv
.rw------- 227k kevin 25 Oct 17:09 Subject_Centric_10252021.xlsx
.rw------- 903k kevin 25 Oct 17:09 Timepoint_Centric_10252021.xlsx
```

## Metadata issues identified

1. for `ageLabel`, the only values are "1 and under", "prenatal", or missing.
   Let's just leave that column out of future exports
2. "simple_race" still has mix of `missing`, "Unknown", and "Decline to answer",
   also both "Mixed" and "Mixed Race"
3. `everBreastFed` is totally wrong. Only 18 non-missing values, all 0
