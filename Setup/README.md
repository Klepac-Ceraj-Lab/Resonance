# Repository setup

Scripts in the directory check for the presence of necessary data files,
and if they don't exit, downloads them to an appropriate location.

## Requirements

### Julia

The code in this repository should work with julia v1.8 or higher.
There are several ways to install julia,
we recommend using [`juliaup`](https://github.com/JuliaLang/juliaup),
or downloading and installing from the [julia website](https://julialang.org/downloads/).

You should be able to open a terminal and type `julia` to open a julia REPL:

```
❯ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.8.2 (2022-09-29)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

type `exit()` and press <kbd>Enter</kbd> to return to your terminal

### Disk space

<!-- TODO: add disk spaces requirements -->

The minimal data to run the repository code
requires ~4 Gb of space.
Running all of the code generates an additional ~5 GB.

By default, scripts will place all data
into your current working directory in the following directories:

- `input/`: This is where the initial input files will be downloaded and read from.
- `scratch/`: Intermediate processed files that are generated from inputs are placed here.
- `models/`: Serialized ML models trained on input data will be place here.
- `figures/`: Main text and supplementary figures from the manuscript will be placed here.
- `tables/`: Main text and supplementary tables from the manuscript will be placed here.

Each of these locations can also be specified using environmental variables.
For example, if you have a large hard drive that would be a better location
for downloads and ML models, you can set

```
INPUT_FILES=/some/big/harddrive/input
MODEL_FILES=/some/big/harddrive/models
```

If you wish to change the locations we recommend using `direnv` or some other method
of reproducibly setting environmental variables,
since if the scripts do not find specific files,
they will be re-downloaded or re-generated.

## Using this julia project

This code is known compatible with julia v1.8 and up.
Project dependencies are provided in the `Project.toml` file.

Once julia is installed, you can download this repository
using `git`, or by selecting the most recent release
from the [releases page](https://github.com/Klepac-Ceraj-Lab/Resonance/releases).
Then, run julia, and [activate this project](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).
If julia is installed in your system `PATH`,
you can activate it directly from the command line using `julia --project=@.`. 
For example, if you downloaded and unpacked the code in your `Documents` directory,
you can run:

```
$ cd Documents/Resonance

$ julia --project=@.
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.2 (2023-07-05)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

Alternatively, you can change directories from the julia REPL,
and then activate the project from the REPL:

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.2 (2023-07-05)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> cd(expanduser("~/Documents/Resonance))

julia> # press ] to enter the Pkg REPL

(@v1.8) pkg> activate .
Activating project at `~/Documents/Resonance`

(Resonance) pkg>
```

Either way, once you have the project activated,
run `instantiate` from the `Pkg` REPL

```
(Resonance) pkg> instantiate
  Installed SIMDDualNumbers ────────────────── v0.1.1
  Installed StatsFuns ──────────────────────── v1.1.1
  Installed CategoricalDistributions ───────── v0.1.9
  Installed JpegTurbo_jll ──────────────────── v2.1.2+0
  #... etc
```


## Repo Setup

Before code is run for the first time,
you should run the scripts found in the `setup/bin` directory.
These will download required data
and build the data structures and models needed
to run the analysis and figure-generating notebooks.

1. `init_repo.jl`: downloads all of the raw data necessary for replicating the analyses.
   Those data are stored on [`OSF.io`](https://doi.org/10.17605/OSF.IO/YBS32)
2. `derive.jl`: uses the data downloaded in step (1) and processes it to abe a bit nicer for downstream analysis.
