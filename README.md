# microclimfPara

Fast, mechanistic microclimate modelling above, below, or within vegetation
canopies on gridded spatial data — built on
[`microclimf`](https://github.com/ilyamaclean/microclimf) (Maclean 2026),
with parallelisation and a number of bug fixes and memory-conservation
improvements layered on top.

## Installation

```r
require(devtools)
install_github("BSpen0820/microclimfPara")
```

## Quick start

```r
library(microclimfPara)
library(terra)

# Run the point microclimate model with the package's bundled example data
micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)

# Subset to representative hours (e.g. hottest/coldest day of each month)
micropoint_mx <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")

# Run the grid model 5cm above ground (add parallel = TRUE, ncores = N to parallelise)
mout_mx <- runmicro(micropoint_mx, reqhgt = 0.05, vegp, soilc, dtmcaerth)
```

See `vignette("running-microclimf")` for the full walkthrough (input data
formats, vegetation/soil parameters, the snow model, running over large
areas, and bioclim variables), or the docs (`vignettes/`) generated on
install.

## Key functions

| Function | Purpose |
|---|---|
| `runpointmodel()` / `runpointmodela()` | Point microclimate model (data.frame / array weather input) |
| `subsetpointmodel()` / `subsetpointmodela()` | Subset point-model output to representative hours before running the grid model |
| `runmicro()` / `runmicro_big()` | Gridded microclimate model |
| `runsnowmodel()` / `subsetsnowmodel()` | Snow accumulation/melt model and its point-model subsetting |
| `vegpfromhab()` | Derive vegetation parameters from habitat classification |
| `checkinputs()` | Validate weather/vegetation/soil/terrain inputs before running a model |
| `writetonc()` | Write gridded output to NetCDF |

Pass `parallel = TRUE, ncores = N` to `runmicro()`/`runsnowmodel()` to use
the TBB-backed parallel grid/snow runners instead of the serial C++ path.

## How this compares to the original microclimf

This fork's output is validated directly against the original
`ilyamaclean/microclimf` package — not just internally. **`test/`
contains the validation suite and, most usefully for a new user or
reviewer, `test/Validation_Results.csv` and
`test/Validation_Results_byvariable.csv`: checked-in results from the
latest such comparison, including a `note` column that calls out every
variable with a confirmed, investigated cause of divergence, with a
worked numerical example and a verdict on which package's value is
physically reasonable in that case.** Start with `test/README.md` for the
full rundown — what's validated, how to re-run it yourself, and why some
values differ from the original (short version: the differences found so
far trace to real, already-fixed bugs in the original package, not this
fork).

## License

GPL (>= 2) — see `LICENSE.md`.
