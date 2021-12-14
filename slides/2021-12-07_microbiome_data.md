---
title: Microbiome data types
author: Kevin Bonham, PhD
date: "2021-12-07"
notes: "For presenation to Khula team meeting"
institute: "Wellesley College"
theme: "Frankfurt"
colortheme: "beaver"
fonttheme: "professionalfonts"
mainfont: "Hack Nerd Font"
fontsize: 12pt
urlcolor: red
linkstyle: bold
aspectratio: 169
# titlegraphic: img/aleph0.png
# logo: img/aleph0-small.png
date: "2021-12-17"
lang: en-US
section-titles: false
toc: true
---

# Design and primary data

## Sample collection design

![](https://i.imgur.com/ESDSLTw.png){ height=90% width=90% }

## Samples
- Stool samples (~500mg) put in buffer as quickly as possible

. . .

- 2 (or 3??) collection types:
  - Zymo DNA/RNA protect - proprietary buffer to stabilize nucleic acids (for sequencing)
  - Ethanol - for metabolomics
  - Direct freeze - for culturing

## Primary Data types

- Shotgun metagenomic sequencing: FASTQ files (sequences + quality scores)
  - paired-end reads, 2x150 bp
  - ~10M reads / sample

. . .

- Metabolomics (LCMS)
  - 4 column types that target different molecule types
  - chromatograph with peaks with m / z & retention time

## Reads per sample

![](https://i.imgur.com/TsCPuUg.png){ height=90% width=90% }


# Derived data types

## Feature Profiles

- Shotgun metagenomics
  - Taxonomic profiles: "Who's there?" - relative abundance of taxa (eg species, genera) in each sample
  - Functional profiles: "What can they do?" - relative abundance of genes (some stratified by species)

. . .

- Metabolomic profiles: "What have they (and we) done?"
  - relative abundance of metabolites, ~5% known

## Pipeline

![](slides/assets/pipeline.png){ height=90% width=90% }

## Shotgun metagenomics profiles

- Reads are aligned to reference database to identify "marker genes" for taxa
- Reduced gene database for identified taxa is generated
  - Reads aligned to reduced database
  - unexplained reads are aligned to all-gene database (translated search)

## Shotgun metagenomics profiles

![](slides/assets/comm_profile.png){ height=90% width=90% }

## Expected taxonomic diversity

![](https://i.imgur.com/UQqFd4C.png){ height=90% width=90% }


## Expected taxonomic diversity - by age

![](https://i.imgur.com/uN9336C.png){ height=90% width=90% }

## Expected taxonomic diversity - by age

![](https://i.imgur.com/PgA5D73.png){ height=90% width=90% }

## Expected functional diversity - by age

![](https://i.imgur.com/V2HqYEY.png){ height=90% width=90% }

## Metabolomics

![](slides/assets/lcms.png){ height=90% width=90% }

<!-- https://en.wikipedia.org/wiki/Liquid_chromatography%E2%80%93mass_spectrometry -->

## Table of integrated areas

![](https://i.imgur.com/kehk1hD.png){ height=90% width=90% }

- ~500 named metabolites
- ~70,000 unnamed metabolites

## Limitations

- Sparsity: most features are not in most samples, most samples don't have most features

. . .

- Heteroskedasticity: variances of features are not uniform

. . .

- Compositionality: features sum to 1 (so all features are dependent on all others)

. . .

- High dimensionality
  - hundreds to thousands in taxonomic profiles
  - tens of thousands to millions of genes in functional profiles
  - tens of thousands of metabolites

## Extra limitations in kids

![](https://i.imgur.com/1A4KEFU.png){ height=90% width=90% }