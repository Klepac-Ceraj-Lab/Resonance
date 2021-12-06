---
title: Microbiome data types and analyses
author: Kevin Bonham, PhD
date: "2021-12-07"
---

## Sample collections

- Stool samples (~500mg) put in buffer as quickly as possible
- 2 (or 3??) collection types:
  - Zymo DNA/RNA protect - proprietary buffer to stabilize nucleic acids (for sequencing)
  - Ethanol - for metabolomics
  - Direct freeze - for culturing

## Primary Data types

- Shotgun metagenomic sequencing
  - paired-end reads, 2x150 bp
  - ~10M reads / sample
  - FASTQ files (sequences + quality scores)
- Metabolomics (LCMS)
  - 4 column types that target different molecule types
  - chromatograph with peaks with Mz / ellution time

## Derived data types

- Shotgun metagenomics
  - Taxonomic profiles: relative abundance of taxa (eg species, genera) in each sample
  - Functional profiles: relatie abundance of genes (some stratified by species)
- Metabolomic profile
  - relative abundance of metabolites, ~5% known

## Shotgun metagenomics profiles

<!-- TODO: insert from Microbiome.jl -->

### Expected taxonomic diversity

<!-- TODO: some charts about taxonomy from ECHO -->

### Expected functional diversity

<!-- TODO: some charts about functions from ECHO -->

## Metabolomics

<!-- TODO: some charts about metabolomics from ECHO -->

## Limitations

- Sparsity
- Heteroschedasticity
- Compositionality
- High dimensionality