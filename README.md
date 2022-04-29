# Resonance Analysis

Analysis code for the RESONANCE microbiome cohort
of the NIH **E**nvironmental influencences on **C**hild **H**ealth **O**utcomes (ECHO) project.

Grant #  NIH UG3 OD023313 (VK-C).

## Repo organization

- Data is stored in the `data/` directory,
  though mostly not tracked via git.
  See [`data/README.md`](data/README.md) for more details.
- This repository is a julia package, with reusable code in the [`src/`](src/) directory.
  For more details, see below.
- Individual analyis "notebooks" are found in the [`notes/`](notes/) directory
  - TODO: Make TOC / build scripts for running code that generates figures
  - The `notes/startup.jl` has some code that is useful to run at the beginning of many notebooks.
- The manuscript draft and references can be found in the [`manuscript/`](manuscript/) directory.

## Running julia code

This code is known compatible with julia v1.7.
Project dependencies are provided in the `Project.toml` file.

TODO: More details on code and data re-use