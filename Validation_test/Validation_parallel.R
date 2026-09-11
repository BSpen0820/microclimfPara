# Validation_parallel.R
#
# Two-part check:
#   Part 1 — Serial path still matches the stored baseline MD5 hashes.
#   Part 2 — Parallel outputs are numerically identical to serial outputs.
#
# Run from the package root:
#   Rscript Validation_test/Validation_parallel.R
# or interactively after devtools::load_all().

devtools::load_all()
library(terra)
library(digest)

# ---- helpers ----------------------------------------------------------------

pass <- function(label) cat(sprintf("  [PASS] %s\n", label))
fail <- function(label, detail = "") {
  cat(sprintf("  [FAIL] %s%s\n", label, if (nzchar(detail)) paste0(" — ", detail) else ""))
}

check_hash <- function(obj, expected, label) {
  got <- digest(obj, algo = "md5")
  if (identical(got, expected)) pass(label) else fail(label, paste("expected", expected, "got", got))
}

check_identical <- function(a, b, label) {
  if (isTRUE(all.equal(a, b))) pass(label) else fail(label, "parallel and serial outputs differ")
}

# ---- shared data ------------------------------------------------------------

cat("\nPreparing data...\n")

data(climdata)
data(vegp)
data(soilc)
data(dtmcaerth)

.rast <- function(m, tem) { r <- rast(m); ext(r) <- ext(tem); crs(r) <- crs(tem); r }
.ta   <- function(x, dtm, xdim = 5, ydim = 5) {
  a <- array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x)))
  .rast(a, dtm)
}

dtm  <- rast(dtmcaerth)
dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)

climdata$temp <- climdata$temp - 12   # cold enough for snow

climarrayr <- list(
  temp      = .ta(climdata$temp,      dtm),
  relhum    = .ta(climdata$relhum,    dtm),
  pres      = .ta(climdata$pres,      dtm),
  swdown    = .ta(climdata$swdown,    dtm),
  difrad    = .ta(climdata$difrad,    dtm),
  lwdown    = .ta(climdata$lwdown,    dtm),
  windspeed = .ta(climdata$windspeed, dtm),
  winddir   = .ta(climdata$winddir,   dtm),
  precip    = .ta(climdata$precip,    dtm)
)

tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")

# Re-use the stored micropointa baseline (avoids re-running the point model)
micropointa <- readRDS("Validation_test/base/micropointa.rds")

cat("Data ready.\n\n")

# ============================================================================
# Part 1: Serial path matches baseline hashes
# ============================================================================

cat("=== Part 1: serial baseline hashes ===\n")

# 1a. micropointa (point model — not affected by parallelisation, but sanity check)
check_hash(micropointa, "9375c48eb328e54b65182db0dff686e7",
           "micropointa hash unchanged")

# 1b. runsnowmodel() serial
cat("  Running runsnowmodel() serial (this takes ~1 min)...\n")
smod_serial <- runsnowmodel(climarrayr, micropointa, vegp, soilc, dtm, dtmc,
                             tme = tme, altcorrect = 0, method = "slow",
                             parallel = FALSE)
check_hash(smod_serial, "924717ec3f7a08924aa037617ac60116",
           "runsnowmodel() serial hash matches baseline")

# 1c. runmicro() serial with snow
cat("  Running runmicro() serial with snow...\n")
mout_serial <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                         altcorrect = 0, snow = TRUE, snowmod = smod_serial,
                         parallel = FALSE)
check_hash(mout_serial, "f8e2d10ac5efe38ef4952aa348dc1a9c",
           "runmicro() serial hash matches baseline")

# 1d. runmicro() serial without snow (data.frame climate — exercises runmicro1Cpp)
micropoint_df <- subsetpointmodel(
  runpointmodel(climdata, 0.05, dtmcaerth, vegp, soilc),
  tstep = "month", what = "tmax"
)
cat("  Running runmicro() serial no-snow (data.frame climate)...\n")
mout_nosnow_serial <- runmicro(micropoint_df, reqhgt = 0.05, vegp, soilc,
                                dtmcaerth, parallel = FALSE)

cat("\n")

# ============================================================================
# Part 2: Parallel outputs numerically identical to serial
# ============================================================================

cat("=== Part 2: parallel == serial ===\n")

# 2a. runsnowmodel() parallel
cat("  Running runsnowmodel() parallel (ncores = 2)...\n")
smod_par <- runsnowmodel(climarrayr, micropointa, vegp, soilc, dtm, dtmc,
                          tme = tme, altcorrect = 0, method = "slow",
                          parallel = TRUE, ncores = 2)
check_identical(smod_serial, smod_par,
                "runsnowmodel() parallel == serial")

# 2b. runmicro() parallel with snow — use smod_serial so we isolate runmicro parallelisation
#     from snow model parallelisation (tests are independent even if 2a fails)
cat("  Running runmicro() parallel with snow (ncores = 2)...\n")
mout_par <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                      altcorrect = 0, snow = TRUE, snowmod = smod_par,
                      parallel = TRUE, ncores = 2)
check_identical(mout_serial, mout_par,
                "runmicro() parallel with snow == serial")

# 2c. runmicro() parallel no-snow (data.frame — exercises runmicro1Par)
cat("  Running runmicro() parallel no-snow, data.frame climate (ncores = 2)...\n")
mout_nosnow_par <- runmicro(micropoint_df, reqhgt = 0.05, vegp, soilc,
                             dtmcaerth, parallel = TRUE, ncores = 2)
check_identical(mout_nosnow_serial, mout_nosnow_par,
                "runmicro() parallel no-snow (df) == serial")

# 2d. runmicro() parallel no-snow (array climate — exercises runmicro2Par)
micropointa_sub <- subsetpointmodela(micropointa, tstep = "month", what = "tmax")
cat("  Running runmicro() parallel no-snow, array climate (ncores = 2)...\n")
mout_arr_serial <- runmicro(micropointa_sub, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                             altcorrect = 0, parallel = FALSE)
mout_arr_par    <- runmicro(micropointa_sub, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                             altcorrect = 0, parallel = TRUE, ncores = 2)
check_identical(mout_arr_serial, mout_arr_par,
                "runmicro() parallel no-snow (array) == serial")

cat("\nDone.\n")
