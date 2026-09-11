# Validation Baseline Report

**Date produced:** 2026-06-10  
**Package version:** microclimf 2.0.1  
**Script:** `Validation_test/Validation_base.R`  
**Purpose:** Establishes numeric baselines for `runsnowmodel()` and `runmicro()` against which future modifications must be compared.

---

## Setup

| Parameter | Value |
|-----------|-------|
| Dataset | `climdata`, `vegp`, `soilc`, `dtmcaerth` (bundled package data, Caerthillian Cove UK) |
| Temperature offset | `climdata$temp - 12` (shifted cold enough to trigger sustained snow) |
| Spatial grid | 50 × 50 cells, dtm unpacked from `PackedSpatRaster` |
| Coarse DTM | `dtmc` = `aggregate(dtm, 10)` (5 × 5 cells) |
| Time series | 8,760 hourly steps (2017, full year) |
| Request height | `reqhgt = 0.05 m` |
| Alt correction | `altcorrect = 0` |

Climate arrays are spatially uniform (each cell repeats the point time series). The spatial variation in outputs is driven entirely by DTM-derived terrain parameters.

---

## Baseline Objects

All four objects are saved to `Validation_test/base/` and should be treated as read-only reference files.

| File | Object | Description |
|------|--------|-------------|
| `micropointa.rds` | `micropointa` | Point model output (class `micropoint`) used as input to both snow and grid models |
| `smod_slow.rds` | `smod_slow` | Snow model output, `method = "slow"` |
| `smod_fast.rds` | `smod_fast` | Snow model output, `method = "fast"` |
| `mout_slow_snow.rds` | `mout_slow_snow` | Grid microclimate output using slow snow model |
| `mout_fast_snow.rds` | `mout_fast_snow` | Grid microclimate output using fast snow model |

---

## Object Structures

### `micropointa` — class `micropoint`, list of 9

```
$ weather : data.frame 8760 × 10
    temp, relhum, pres, swdown, difrad, lwdown, windspeed, precip, winddir, obs_time
$ dfo     : data.frame 8760 × 10
    umu, kp, muGp, DDp, T0p, dtrp, G, soilm, Tg, Tc
$ Tbz     : NA  (reqhgt > 0, no below-ground output)
$ lat     : 50
$ long    : -5.21
$ zref    : 2
$ subs    : int[1:8760]
$ tmeorig : POSIXlt[1:8760]
$ matemp  : -0.639
```

Key point model values at t=1:
- `Tc` = -3.40 °C, `Tg` = -4.34 °C, `soilm` = 0.419

### `smod_slow` / `smod_fast` — list of 6, dim [50, 50, 8760]

```
$ Tc             : num [50, 50, 8760]
$ Tg             : num [50, 50, 8760]
$ groundsnowdepth: num [50, 50, 8760]   (metres)
$ totalSWE       : num [50, 50, 8760]   (mm water equivalent)
$ snowden        : num [50, 50, 8760]   (kg/m³)
$ umu            : num [50, 50, 8760]
```

### `mout_slow_snow` / `mout_fast_snow` — list of 11, dim [50, 50, 8760]

```
$ Tz       : num [50, 50, 8760]   (°C, air temp at reqhgt)
$ tleaf    : num [50, 50, 8760]   (°C)
$ relhum   : num [50, 50, 8760]   (%)
$ soilm    : num [50, 50, 8760]
$ windspeed: num [50, 50, 8760]   (m/s)
$ Rdirdown : num [50, 50, 8760]   (W/m²)
$ Rdifdown : num [50, 50, 8760]   (W/m²)
$ Rlwdown  : num [50, 50, 8760]   (W/m²)
$ Rswup    : num [50, 50, 8760]   (W/m²)
$ Rlwup    : num [50, 50, 8760]   (W/m²)
$ tme      : POSIXct[1:8760]
```

---

## Baseline Summary Statistics

### `smod_slow` and `smod_fast` (snow model outputs)

| Variable | Min | Mean | Max |
|----------|-----|------|-----|
| `Tc` (°C) | -15.1041 | -2.0128 | 8.2096 |
| `Tg` (°C) | -15.1041 | -1.2986 | 76.1575 |
| `groundsnowdepth` (m) | 0.0000 | 0.0897 | 1.6397 |
| `totalSWE` (mm) | -0.0290 | 36.4738 | 367.0626 |
| `snowden` (kg/m³) | 217.0 | 217.0 | 217.0 |
| `umu` | 0.5261 | 1.0841 | 2.4372 |

### `mout_slow_snow` and `mout_fast_snow` (grid microclimate outputs)

| Variable | Min | Mean | Max |
|----------|-----|------|-----|
| `Tz` (°C) | -15.1175 | -0.0730 | 56.3220 |
| `tleaf` (°C) | -15.1041 | 0.2259 | 54.8900 |
| `relhum` (%) | 7.7029 | 97.3586 | 100.0000 |
| `soilm` | 0.1101 | 0.3879 | 0.4190 |
| `windspeed` (m/s) | -4.4350 | 1.3523 | 11.2357 |
| `Rdirdown` (W/m²) | 0.0000 | 72.4055 | 1352.0000 |
| `Rdifdown` (W/m²) | 0.0000 | 46.5085 | 484.6595 |
| `Rlwdown` (W/m²) | 0.0000 | 222.1520 | 600.2587 |
| `Rswup` (W/m²) | 0.0000 | 60.0147 | 868.2266 |
| `Rlwup` (W/m²) | 0.0000 | 209.0360 | 781.5854 |

---

## slow vs fast Comparison

The `"slow"` and `"fast"` methods of `runsnowmodel()` produce **bit-for-bit identical outputs** for this dataset — confirmed by matching MD5 hashes (see below). Any future modification that causes a hash difference between methods should be investigated.

---

## Known Baseline Behaviours

- **`totalSWE` min = -0.029 mm:** A small negative SWE arises from the numerical snow melt scheme in edge cases (near-zero accumulation cells). Pre-existing behaviour.
- **`windspeed` min = -4.435 m/s:** Negative wind speed can occur at below-canopy heights under very stable conditions where `windtiCpp` returns a counter-gradient value. Expected behaviour.
- **`Rlwdown` / `Rlwup` min = 0:** Occurs at timesteps where snow depth exceeds `reqhgt`, placing the output height below the snow surface; the model returns zero radiation in this case.
- **`snowden` constant at 217 kg/m³:** The Taiga snow environment uses fixed density parameters that produce a single density value. Expected.

---

## Baseline Hashes

Hashes are computed with `digest::digest(obj, algo = "md5")` — this hashes the serialized R object directly, bypassing gzip timestamp non-determinism in `.rds` files. The hash will match if and only if the object is numerically identical to the baseline.

| Object | MD5 hash |
|--------|----------|
| `micropointa` | `9375c48eb328e54b65182db0dff686e7` |
| `smod_slow` | `924717ec3f7a08924aa037617ac60116` |
| `smod_fast` | `924717ec3f7a08924aa037617ac60116` |
| `mout_slow_snow` | `f8e2d10ac5efe38ef4952aa348dc1a9c` |
| `mout_fast_snow` | `f8e2d10ac5efe38ef4952aa348dc1a9c` |

`smod_slow` = `smod_fast` and `mout_slow_snow` = `mout_fast_snow` (identical hashes confirm the two methods are exactly equivalent on this dataset).

---

## How to Use for Regression Validation

After making changes to `runsnowmodel()` or `runmicro()`, re-run the validation script and check hashes. **All hashes must match the table above — this is a hard requirement.**

```r
library(digest)
devtools::load_all()
library(terra)

baseline_hashes <- c(
  micropointa    = "9375c48eb328e54b65182db0dff686e7",
  smod_slow      = "924717ec3f7a08924aa037617ac60116",
  smod_fast      = "924717ec3f7a08924aa037617ac60116",
  mout_slow_snow = "f8e2d10ac5efe38ef4952aa348dc1a9c",
  mout_fast_snow = "f8e2d10ac5efe38ef4952aa348dc1a9c"
)

check_hash <- function(obj, name) {
  h <- digest(obj, algo = "md5")
  expected <- baseline_hashes[name]
  status <- if (h == expected) "PASS" else "FAIL <-- outputs changed"
  cat(sprintf("  %-20s  %s  [%s]\n", name, h, status))
}

# Re-run (source Validation_base.R or reproduce manually)
source("Validation_test/Validation_base.R")

cat("\n=== Hash Verification ===\n")
check_hash(micropointa,    "micropointa")
check_hash(smod_slow,      "smod_slow")
check_hash(smod_fast,      "smod_fast")
check_hash(mout_slow_snow, "mout_slow_snow")
check_hash(mout_fast_snow, "mout_fast_snow")
```
