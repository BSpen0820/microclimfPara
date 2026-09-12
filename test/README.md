# test/

Validation scripts for `microclimfPara`. There are two independent tracks —
run whichever one answers the question you actually have.

## Track 1: does this fork match the original microclimf?

This is what most people want: confirmation that `microclimfPara`'s output
(serial and parallel) matches `ilyamaclean/microclimf`.

1. **`01_original_vs_fork.R`** — installs the original `ilyamaclean/microclimf`,
   runs the point/grid/snow models, then installs this fork
   (`BSpen0820/microclimfPara`, from GitHub) and runs the same inputs both
   serial and parallel. Saves all three sets of `.rds` outputs to
   `Base_Results/`, `MicroPar_Ser/`, `MicroPar_Par/`.
2. **`02_compare_results.R`** — hashes every `.rds` produced by step 1 and
   compares each against the matching `Base_Results/` file (exact MD5 match,
   falling back to `all.equal(tolerance = 1e-10)`), writing
   `Validation_Results.csv`. It also breaks the comparison for `mout_nosnow`,
   `mout_snow` and `smod_snow` down **by individual output variable** (`Tz`,
   `tleaf`, `Tc`, `Tg`, ...), writing `Validation_Results_byvariable.csv` —
   useful for seeing exactly which variables differ from the original and by
   how much, rather than a single pass/fail for the whole object.

Both scripts expect to be run **from the `test/` directory** (they read/write
`Base_Results/`, `MicroPar_Ser/`, `MicroPar_Par/` relative to the working
directory):

```r
setwd("test")
source("01_original_vs_fork.R")
source("02_compare_results.R")
```

`01_original_vs_fork.R` installs microclimfPara from GitHub, not the local
working copy — it validates whatever's currently pushed, not uncommitted
local changes.

Pre-generated `Validation_Results.csv` and `Validation_Results_byvariable.csv`
are checked into this directory so you can see the outcome without running
either script yourself. Re-run both scripts and overwrite them if you want to
confirm the comparison against the current `main`.

### Why some values differ from the original

`Validation_Results_byvariable.csv` reports, per variable, the mean/median/
max absolute difference plus a **median** relative difference (`median_rel_diff_pct`
— median, not mean, because a mean is dominated by the rare near-zero-value
cells where any relative-difference metric blows up; see the comment above
`compare_variable()` in `02_compare_results.R`).

It shows the no-snow grid variables (`mout_nosnow`) matching the original
closely — `soilm` is an exact match, `Rdirdown` has the largest mean absolute
difference (1.56 W/m², max ~1009 W/m², see below) but its median relative
difference is still 0%, meaning most cells/hours are essentially unaffected.
The snow-covered variables (`mout_snow`, `smod_snow`) diverge more and more
broadly — including variables (wind speed, all four radiation components)
that aren't obviously related to ground temperature at all — though even
there the median relative difference is often 0% (`Rdirdown`, `Rdifdown`,
`Rswup`) or small (`windspeed` 2.6%, `Tz` 21%), with the larger *mean*
differences (radiation fields: 9–46 W/m², max in the hundreds) reflecting a
smaller tail of more strongly affected cells/hours rather than a uniform
shift. `tleaf`'s 50% median relative difference is the standout exception.
Three known, expected sources explain this pattern:

- **A snow-depth chunk carry-over bug, since fixed (commit `8e32e4b`,
  "Fix snow-model chunk carry-over and floating-point day classification").**
  This fork processes a full run in 5-day chunks (for memory efficiency);
  the original `microclimf` doesn't chunk at all. A copy-paste bug meant the
  carried-over ground snow depth state (`other$isnowdg`) was assigned to the
  wrong variable and immediately clobbered by the next line, so it was never
  actually updated between chunks — ground snow depth collapsed back to ~0
  at every 5-day boundary regardless of how much had actually accumulated.
  This is fixed in the validated build, but its existence is exactly why the
  divergence isn't limited to ground temperature: this model dynamically
  adjusts canopy height, wind profile, and radiation transmission as snow
  accumulates, so anything that changes the snow depth/SWE trajectory
  cascades into wind speed and every radiation component for the whole
  snow-covered grid, not just the below-ground temperature fields. Even with
  the bug fixed, this fork's chunked accumulation and the original's
  unchunked accumulation aren't guaranteed to track identically, which is
  the main reason `mout_snow`'s `windspeed`/`Rdirdown`/`Rdifdown`/`Rlwdown`/
  `Rswup`/`Rlwup` show a real tail of large-magnitude differences (mean
  absolute differences of 9–46 W/m² / 0.28 m/s, max values in the hundreds)
  from the original that the no-snow grid doesn't, even though most
  individual cells/hours still match closely.
- **The snow ground-temperature smoothing fix (see Track 2 below).** It
  deliberately changes below-ground temperature under snow cover — that's
  the whole point of the fix (it corrects a staircase artifact and an
  unblended bare-soil/snow regime switch present in the original). Shows up
  mainly in `smod_snow$Tc`/`Tg` and `mout_snow$Tz`.
- **Floating-point noise from parallelization.** `microclimfParallel.cpp`
  and `microclimfCpp.cpp` are separate translation units, so the compiler
  makes different inlining/register-allocation decisions for each — parallel
  output can differ from serial by up to a few ULP (~2–4×10⁻¹⁶), which is
  physically meaningless but means an exact hash match isn't guaranteed run
  to run (see CLAUDE.md, "Parallel vs serial floating-point tolerance").
  This is a separate, much smaller effect than the two sources above, and
  only affects `MicroPar_Par` vs. `MicroPar_Ser`, not the comparison against
  `Base_Results`.

## Track 2: snow ground-temperature smoothing correctness

### The problem

`belowpointsnow()` computes below-ground temperature under snow cover using
a smoothed proxy of snow-surface temperature (`Tzd`) at the query depth. The
original implementation built `Tzd` as a **day-block-constant mean** — one
flat value per calendar day. Two real bugs fell out of that:

- **Daily "staircase" artifact.** Because `Tzd` was constant within a day and
  jumped to a new constant at midnight, below-ground temperature showed an
  artificial discontinuity every day boundary, confirmed against real
  station data (Grand Targhee, Tetons, WY, USA): a ~9°C staircase jump on
  Jan 8–9, 2018 that has no physical basis — the ground doesn't reset at
  midnight.
- **Unblended bare-soil/snow-covered regime switch.** The model hard-switched
  from bare-soil to snow-covered physics the instant *any* local snow
  appeared, however thin, with no allowance for thermal inertia — producing
  an unrealistic single-hour cliff (e.g. a 9°C jump at the May 23, 2018
  station snow-onset event) instead of a gradual transition.

### The fix

The snow-covered ground temperature now uses *the same* depth-dependent
trailing rolling window formula the no-snow model's own `Tbelowgroundv()`
"complete" case already uses (`n = round(-118.35 * reqhgt / meanD)`,
`manCpp(series, n)`), instead of a bespoke smoother — so below-ground
temperature is computed consistently whether or not the pixel is
snow-covered. Doing that safely required two pieces:

- **`.build_snow_Tgref()`** (`R/internal.R`) assembles a full-year, NA-free
  reference series *before* smoothing: it splices snow-surface temperature
  (`Tg`) onto snowdays/mixed-day hours and the no-snow model's own
  below-ground prediction (`moutn$Tz`) onto pure-nosnow-day hours, then fills
  any remaining gap by linear interpolation from that pixel's own nearest
  valid hours. This is what makes it safe to hand the series to `manCpp()`
  directly — it's genuinely calendar-complete, so there's no
  season-end/season-start wraparound for a circular moving-average window to
  smear (`manCpp()`'s window is circular, so an incomplete series would wrap
  into itself).
- **`snowdaymov()`** (`src/microclimfCpp.cpp`) applies that rolling window to
  the assembled `Tgref` series, then slices the smoothed result back down to
  just the requested snowdays hours.

Two edge cases surfaced during review and were fixed as part of this design:
a pixel with fewer than 2 valid hours all year can't be interpolated
(`stats::approx()` requires ≥2 points) and is blanked entirely rather than
crashing or silently passing through a single stray value; and a
`tsteps %% 24 == 0` day-alignment guard ensures a violated invariant errors
loudly instead of corrupting tail hours with `manCpp()`'s `n > 48` day-block
path.

### The scripts

Each script below is self-contained and can be run independently, in
increasing order of scope:

3. **`03_snow_tgref_build.R`** — unit checks for `.build_snow_Tgref()`:
   splicing `Tg`/`moutn$Tz` onto the right calendar days (including a
   genuine "mixed day" flagged as both a snowday and a nosnowday, which must
   use `Tg`, never `moutn$Tz`), filling isolated gaps by interpolation, and
   blanking pixels with too few valid hours to interpolate safely.
4. **`04_snow_daymov.R`** — unit checks for `snowdaymov()`: confirms its
   output matches `manCpp(Tgref, n)` directly (the same formula
   `Tbelowgroundv()` uses), sliced to the requested snowdays hours; plus the
   masked-pixel, non-positive-`meanD`, and `n >= tsteps` fallback edge cases.
5. **`05_snow_wiring_consistency.R`** — end-to-end check that this plumbing
   (steps 3 and 4, threaded through `runmicro()`'s snow path via
   `.runmicrosnow1()` → `gridmicrosnow1()`/`gridmicrosnow1Par()`) gives
   identical results serial vs. parallel, using the package's bundled
   example data.

Run these from the **package root**:

```r
devtools::load_all()
source("test/03_snow_tgref_build.R")
source("test/04_snow_daymov.R")
source("test/05_snow_wiring_consistency.R")
```

or `Rscript test/0N_....R` from the package root.
