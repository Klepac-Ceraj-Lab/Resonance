# Resonance Analysis

Analysis code for the RESONANCE microbiome cohort
of the NIH **E**nvironmental influencences on **C**hild **H**ealth **O**utcomes (ECHO) project.

Grant #  NIH UG3 OD023313 (VK-C).

[![zenodo DOI](https://zenodo.org/badge/448058945.svg)](https://zenodo.org/badge/latestdoi/448058945)
[![OSF.io DOI](https://img.shields.io/badge/OSF.io-10.17605%2FOSF.IO%2FYBS32-informational)](https://doi.org/10.17605/OSF.IO/YBS32)

## Repo organization

- This repository is a julia package, with reusable code in the [`src/`](src/) directory.
  For more details, see below.
- Before running analysis code, your environment must be set up.
  Please see instructions in [`setup/README.md`](setup/README.md)
- Individual analyis "notebooks" are found in the [`manuscript/`](manuscript/) directory
- The manuscript draft and references can be also be found in the [`manuscript/`](manuscript/) directory.

## Running the code

Helper scripts for downloading data and running code necessary for reproducing these analyses
are available in the `Setup/` directory.
See [`Setup/README.md`](Setup/README.md) for additional information
on how to install julia and this repository.
