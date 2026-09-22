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

### Is the largest `max_abs_diff` a bug in this fork?

No — investigated case by case (2026-09-22), the largest outliers are
confirmed bugs in the **original** `microclimf`, already fixed in this fork.
`microclimfPara`'s values are the physically reasonable ones:

- **Radiation/temperature (`mout_nosnow`/`mout_snow`, up to `Tz` 21–56 degC,
  `tleaf` 34–51 degC, `relhum` 59–67 pts, radiation fields up to ~1009–1352
  W/m²).** At the worst `mout_nosnow$Tz` cell/hour (2017-10-06 07:00,
  near sunrise), total incoming shortwave (`swdown`) was only 62.5 W/m², yet
  the *original* reports `Rdirdown` (direct beam alone) of 743–816 W/m² —
  over 12x the total energy that arrived, which violates energy conservation
  and is not physically possible. `microclimfPara` correctly reports
  `Rdirdown = 0` there (the already-documented near-horizon
  `cos(zenith) < 0.065` amplification bug — see "Direct-beam radiation fix"
  in `CLAUDE.md`), and its `Tz` of ~9.49 degC sits almost exactly on the raw
  ambient air temperature for that hour (9.448 degC) — exactly what's
  expected right after sunrise with genuinely weak sun. The affected cells
  cluster in one interior patch of the terrain (rows 2–13, columns 27–34 of
  the 50x50 grid — not the domain edge), consistent with a slope/aspect/
  horizon-shading combination that's more prone to the near-horizon
  division-by-~0 that causes the bug.
- **Snow ground temperature (`smod_snow$Tg` up to 76 degC, `Tc` up to
  14 degC, `totalSWE`/`groundsnowdepth`).** At the worst `Tg` cell/hour
  (2017-05-24 noon, full sun), the *original*'s `groundsnowdepth` reads
  exactly `0.0000` while its own `totalSWE` still shows ~5.4 units of snow
  water equivalent at the same timestep — an internal contradiction, since
  `groundsnowdepth` and `totalSWE` are tracked as independent state in the
  snow model and can drift apart. That triggers a spurious full bare-ground
  regime switch, and the *original*'s ground temperature spikes to a
  physically implausible 76 degC in bright May sun before decaying back
  over several hours. `microclimfPara`'s depth-proportional regime blending
  (commit `c070816`, "Fix snow-covered below-ground staircase and unblended
  regime switch") avoids this; its `Tg` stays a sane, flat `0.000` degC
  throughout.

**How rare is this?** Out of 20,778,720 compared `mout_nosnow$Tz` cell-hours:

| Percentile | Abs. difference |
|---|---|
| 50% (median) | ~0 (machine precision) |
| 90% | 0.000005 degC |
| 99% | 0.0003 degC |
| 99.9% | 0.03 degC |
| 99.99% | 0.81 degC |
| 99.999% | 5.7 degC |
| 100% (max) | 21.2 degC |

Only 70 of those 20,778,720 cell-hours (0.00034%) differ by more than 10
degC, all tied to the same near-horizon radiation bug. 99.99% of the grid
agrees with the original to within a degree.

`Validation_Results_byvariable.csv`'s `note` column flags every variable
row where a specific cause has been investigated and confirmed, so the raw
`max_abs_diff`/`mean_abs_diff` numbers aren't read in isolation.

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
