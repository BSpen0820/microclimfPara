# Validation_snow_wiring_consistency.R
# End-to-end serial vs parallel consistency check for the new Tgref/snowdays
# plumbing threaded through .runmicrosnow1() (data.frame climate) ->
# gridmicrosnow1()/gridmicrosnow1Par(). Uses the package's own bundled
# example datasets (small, fast, no saved fixture files needed) rather than
# a stored fixture - avoids relying on stale fixture provenance (a fixture
# built with a mismatched reqhgt caused a hard-to-diagnose crash earlier in
# this investigation).
#
# NOTE: runsnowmodel() segfaults on a small *cropped* dtmcaerth (confirmed
# directly while designing this test - a pre-existing bug unrelated to this
# change, out of scope for this plan). Use the full, uncropped 50x50
# dtmcaerth grid, which is confirmed to run in a few seconds.
# Run from the package root: Rscript test/Validation_snow_wiring_consistency.R

library(microclimfPara)
library(terra)

data(dtmcaerth); data(climdata); data(vegp); data(soilc)
dtm <- terra::unwrap(dtmcaerth)

# 20 contiguous days (480 hours), cooled by 8 degC per the package's own
# documented pattern (see ?subsetsnowmodel) to guarantee a snow/no-snow mix.
weather <- climdata[1:480, ]
weather$temp <- weather$temp - 8

reqhgt <- -0.15
micropoint <- runpointmodel(weather, reqhgt = reqhgt, dtm = dtm, vegp = vegp,
                             soilc = soilc, runchecks = FALSE)
smod <- runsnowmodel(weather, micropoint, vegp, soilc, dtm)
stopifnot(any(smod$totalSWE > 0, na.rm = TRUE))   # confirms a snow-covered day exists
stopifnot(any(smod$totalSWE == 0, na.rm = TRUE))  # confirms a snow-free day exists

mout_serial <- runmicro(micropoint, reqhgt, vegp, soilc, dtm, snow = TRUE,
                         snowmod = smod, runchecks = FALSE, parallel = FALSE)
mout_parallel <- runmicro(micropoint, reqhgt, vegp, soilc, dtm, snow = TRUE,
                           snowmod = smod, runchecks = FALSE, parallel = TRUE, ncores = 2)

na_serial <- is.na(mout_serial$Tz)
na_parallel <- is.na(mout_parallel$Tz)
stopifnot(identical(na_serial, na_parallel))
cat("PASS: serial and parallel produce NA at exactly the same positions.\n")

maxdiff <- max(abs(mout_serial$Tz - mout_parallel$Tz), na.rm = TRUE)
cat(sprintf("max |serial - parallel| Tz difference (non-NA cells): %.10f\n", maxdiff))
stopifnot(maxdiff < 1e-8)
cat("PASS: serial and parallel Tgref/snowdays plumbing agree on Tz.\n")

cat("\nPASS: end-to-end snow-model wiring (Tgref assembly -> gridmicrosnow1/1Par) is serial/parallel consistent.\n")
