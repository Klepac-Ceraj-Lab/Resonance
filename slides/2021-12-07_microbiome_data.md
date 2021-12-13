---
title: Microbiome data types
author: Kevin Bonham, PhD
date: "2021-12-07"
notes: "For presenation to Khula team meeting"
---

## Sample collection design

<p class="stretch"><img src="../assets/study-design.png"></p>

## Samples and primary data

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
  - chromatograph with peaks with $m / z$ & retention time

## Reads per sample

<p class="stretch"><img src="../assets/kneaddatacounts.png"></p>


## Derived data types

- Shotgun metagenomics
  - Taxonomic profiles: "Who's there?" - relative abundance of taxa (eg species, genera) in each sample
  - Functional profiles: "What can they do?" - relative abundance of genes (some stratified by species)

. . .

- Metabolomic profiles: "What have they (and we) done?"
  - relative abundance of metabolites, ~5% known

## Pipeline

<p class="stretch"><img src="../assets/pipeline.png"></p>

## Shotgun metagenomics profiles

- Reads are aligned to reference database to identify "marker genes" for taxa
- Reduced gene database for identified taxa is generated
  - Reads aligned to reduced database
  - unexplained reads are aligned to all-gene database (translated search)

## Shotgun metagenomics profiles

<p class="stretch"><img src="../assets/comm_profile.png"></p>

## Expected taxonomic diversity

<p class="stretch"><img src="../assets/echo_richness.png"></p>


## Expected taxonomic diversity - by age

<p class="stretch"><img src="../assets/echo_richness_byage.png"></p>

## Expected taxonomic diversity - by age

<p class="stretch"><img src="https://i.imgur.com/PgA5D73.png"></p>

## Expected functional diversity - by age

<p class="stretch"><img src="https://i.imgur.com/V2HqYEY.png"></p>

## Metabolomics

<p class="stretch"><img src="https://upload.wikimedia.org/wikipedia/en/thumb/7/71/Liquid_chromatography_MS_spectrum_3D_analysis.png/1920px-Liquid_chromatography_MS_spectrum_3D_analysis.png" alt="https://en.wikipedia.org/wiki/Liquid_chromatography%E2%80%93mass_spectrometry"></p>

## Table of integrated areas

<p class="stretch"><img src="https://i.imgur.com/kehk1hD.png"></p>

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

<p class="stretch"><img src="../assets/age_taxon_ratios.png"></p>