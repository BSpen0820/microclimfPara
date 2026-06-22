# =============================================================================
# Validation: memory-reduction changes preserve model results
# =============================================================================
# Confirms that the recent memory work does NOT change microclimate outputs:
#   * output-array gating in the parallel C++ runners
#     (runmicro1-4Par, gridmicrosnow1/2Par)
#   * rm() + silent gc() added through runmicro() / runsnowmodel() internals
#
# Checks performed (bundled data only - no network installs):
#   A. runmicro no-snow, data.frame climate : serial == parallel
#   B. output gating : a restricted `out=` returns values identical to the full
#      run for the requested layer, and omits the unrequested layers
#   C. runsnowmodel + snow runmicro, data.frame : serial == parallel (fast)
#   D. runsnowmodel, data.frame : serial == parallel (slow)
#   E. runmicro no-snow, array climate : serial == parallel
#   F. runsnowmodel + snow runmicro, array climate : serial == parallel
#
# Comparison tolerance is 1e-10. Parallel C++ output differs from serial by
# <2 ULP because the two translation units inline differently; this is
# physically meaningless (see CLAUDE.md, "Parallel vs serial floating-point
# tolerance"). Use all.equal(tolerance=1e-10), never identical().
#
# Usage:  Rscript test/Validation_memfix.R     (run from the package root)
#         or devtools::load_all() then source this file.
# =============================================================================

suppressMessages({
  if (requireNamespace("pkgload", quietly = TRUE) &&
      !"microclimfPara" %in% loadedNamespaces()) {
    # Prefer the working-tree sources so the script tests local edits.
    try(pkgload::load_all(".", quiet = TRUE), silent = TRUE)
  }
  if (!"microclimfPara" %in% loadedNamespaces()) library(microclimfPara)
  library(terra)
})

data(climdata)
data(vegp)
data(soilc)
data(dtmcaerth)

NCORES <- 2L
TOL    <- 1e-10

# ---- result collection ------------------------------------------------------
results <- data.frame(check = character(), description = character(),
                      pass = logical(), detail = character(),
                      stringsAsFactors = FALSE)

add_result <- function(check, description, comparison, detail = "") {
  ok <- isTRUE(comparison)
  if (!ok && is.character(comparison)) {          # all.equal returns a message
    detail <- paste(comparison, collapse = "; ")
  }
  results <<- rbind(results,
                    data.frame(check = check, description = description,
                               pass = ok, detail = detail,
                               stringsAsFactors = FALSE))
  cat(sprintf("[%s] %-58s %s\n", check,
              substr(description, 1, 58),
              if (ok) "PASS" else "FAIL"))
  invisible(ok)
}

# all.equal returns TRUE or a character vector describing the difference.
eq <- function(a, b) all.equal(a, b, tolerance = TOL)

# =============================================================================
# Shared point models (data.frame climate)
# =============================================================================
cat("\n--- Preparing data.frame point models ---\n")
mp_nosnow <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)
mp_nosnow <- subsetpointmodel(mp_nosnow, tstep = "month", what = "tmax")

mp_snow <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)
mp_snow <- subsetpointmodel(mp_snow, tstep = "month", what = "tmin")

# =============================================================================
# A. runmicro, no snow, data.frame  -- serial vs parallel
#    Exercises .runmicronosnow -> .runmodel1Cpp / runmicro1Par
# =============================================================================
cat("\n--- A: runmicro no-snow (data.frame) ---\n")
A_ser <- runmicro(mp_nosnow, reqhgt = 0.05, vegp, soilc, dtmcaerth)
A_par <- runmicro(mp_nosnow, reqhgt = 0.05, vegp, soilc, dtmcaerth,
                  parallel = TRUE, ncores = NCORES)
add_result("A", "no-snow runmicro: parallel == serial", eq(A_ser, A_par))

# =============================================================================
# B. Output gating  -- the core test for the parallel allocation fix
#    A restricted `out=` must (i) return the requested layer bit-for-bit equal
#    to the full run, and (ii) omit the unrequested layers.
# =============================================================================
cat("\n--- B: output gating (parallel) ---\n")
out_tz <- c(TRUE, rep(FALSE, 9))   # request Tz only
B_tz <- runmicro(mp_nosnow, reqhgt = 0.05, vegp, soilc, dtmcaerth,
                 parallel = TRUE, ncores = NCORES, out = out_tz)
add_result("B1", "gated Tz identical to full-run Tz", eq(A_par$Tz, B_tz$Tz))
# Unrequested layers must be absent (memory not allocated/returned).
omitted <- !any(c("tleaf", "relhum", "windspeed", "Rdirdown", "Rlwup") %in%
                names(B_tz))
add_result("B2", "unrequested output layers omitted", omitted,
           detail = paste("names:", paste(names(B_tz), collapse = ",")))

# Serial gating should behave the same way (sanity check on parity).
B_tz_ser <- runmicro(mp_nosnow, reqhgt = 0.05, vegp, soilc, dtmcaerth,
                     out = out_tz)
add_result("B3", "gated Tz: parallel == serial", eq(B_tz_ser$Tz, B_tz$Tz))

# =============================================================================
# C. Snow, data.frame, method "fast" -- serial vs parallel
#    Exercises .snowmodelq1 / gridmodelsnow1Par and
#    .runmicrosnow1 / gridmicrosnow1Par (and their rm()/gc()).
# =============================================================================
cat("\n--- C: snow fast (data.frame) ---\n")
C_smod_ser <- runsnowmodel(climdata, mp_snow, vegp, soilc, dtmcaerth,
                           method = "fast")
C_smod_par <- runsnowmodel(climdata, mp_snow, vegp, soilc, dtmcaerth,
                           method = "fast", parallel = TRUE, ncores = NCORES)
add_result("C1", "runsnowmodel fast: parallel == serial",
           eq(C_smod_ser, C_smod_par))

C_mout_ser <- runmicro(mp_snow, reqhgt = 0.05, vegp, soilc, dtmcaerth,
                       snow = TRUE, snowmod = C_smod_ser)
C_mout_par <- runmicro(mp_snow, reqhgt = 0.05, vegp, soilc, dtmcaerth,
                       snow = TRUE, snowmod = C_smod_par,
                       parallel = TRUE, ncores = NCORES)
add_result("C2", "snow runmicro: parallel == serial",
           eq(C_mout_ser, C_mout_par))

# =============================================================================
# D. Snow, data.frame, method "slow" -- serial vs parallel
#    Exercises .snowmodel1 (full grid + subset) and its per-chunk rm()/gc().
# =============================================================================
cat("\n--- D: snow slow (data.frame) ---\n")
D_smod_ser <- runsnowmodel(climdata, mp_snow, vegp, soilc, dtmcaerth,
                           method = "slow")
D_smod_par <- runsnowmodel(climdata, mp_snow, vegp, soilc, dtmcaerth,
                           method = "slow", parallel = TRUE, ncores = NCORES)
add_result("D", "runsnowmodel slow: parallel == serial",
           eq(D_smod_ser, D_smod_par))

# =============================================================================
# Array climate setup (covers the *2 / *4 internal paths)
# =============================================================================
cat("\n--- Preparing array climate inputs ---\n")
.rast <- function(m, tem) {
  r <- rast(m); ext(r) <- ext(tem); crs(r) <- crs(tem); r
}
.ta <- function(x, dtm, xdim = 5, ydim = 5) {
  a <- array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x)))
  .rast(a, dtm)
}
dtm  <- rast(dtmcaerth)
dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)
tme  <- as.POSIXlt(climdata$obs_time, tz = "UTC")

mk_array <- function(df) list(
  temp = .ta(df$temp, dtm),       relhum = .ta(df$relhum, dtm),
  pres = .ta(df$pres, dtm),       swdown = .ta(df$swdown, dtm),
  difrad = .ta(df$difrad, dtm),   lwdown = .ta(df$lwdown, dtm),
  windspeed = .ta(df$windspeed, dtm), winddir = .ta(df$winddir, dtm),
  precip = .ta(df$precip, dtm))

climarr      <- mk_array(climdata)
climdata_cld <- climdata; climdata_cld$temp <- climdata$temp - 12  # cold enough for snow
climarr_snow <- mk_array(climdata_cld)

# runpointmodela is left serial (parallel path there was intentionally untouched)
mpa_nosnow <- runpointmodela(climarr,      tme, reqhgt = 0.05, dtm, vegp, soilc)
mpa_snow   <- runpointmodela(climarr_snow, tme, reqhgt = 0.05, dtm, vegp, soilc)

# =============================================================================
# E. runmicro, no snow, array  -- serial vs parallel
#    Exercises .runmodel2Cpp / runmicro2Par
# =============================================================================
cat("\n--- E: runmicro no-snow (array) ---\n")
E_ser <- runmicro(mpa_nosnow, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0)
E_par <- runmicro(mpa_nosnow, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0,
                  parallel = TRUE, ncores = NCORES)
add_result("E", "no-snow runmicro (array): parallel == serial", eq(E_ser, E_par))

# =============================================================================
# F. Snow, array  -- serial vs parallel
#    Exercises .snowmodel2 / gridmodelsnow2Par and
#    .runmicrosnow2 / gridmicrosnow2Par (largest per-chunk objects).
# =============================================================================
cat("\n--- F: snow (array) ---\n")
F_smod_ser <- runsnowmodel(climarr_snow, mpa_snow, vegp, soilc, dtm, dtmc,
                           tme = tme, altcorrect = 0, method = "slow")
F_smod_par <- runsnowmodel(climarr_snow, mpa_snow, vegp, soilc, dtm, dtmc,
                           tme = tme, altcorrect = 0, method = "slow",
                           parallel = TRUE, ncores = NCORES)
add_result("F1", "runsnowmodel (array): parallel == serial",
           eq(F_smod_ser, F_smod_par))

F_mout_ser <- runmicro(mpa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                       altcorrect = 0, snow = TRUE, snowmod = F_smod_ser)
F_mout_par <- runmicro(mpa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                       altcorrect = 0, snow = TRUE, snowmod = F_smod_par,
                       parallel = TRUE, ncores = NCORES)
add_result("F2", "snow runmicro (array): parallel == serial",
           eq(F_mout_ser, F_mout_par))

# =============================================================================
# Summary
# =============================================================================
cat("\n=============================== SUMMARY ===============================\n")
print(results[, c("check", "description", "pass")], row.names = FALSE)
try(write.csv(results, "Validation_memfix_results.csv", row.names = FALSE),
    silent = TRUE)

n_fail <- sum(!results$pass)
cat(sprintf("\n%d/%d checks passed.\n", sum(results$pass), nrow(results)))
if (n_fail > 0) {
  cat("FAILURES:\n")
  print(results[!results$pass, c("check", "description", "detail")], row.names = FALSE)
  quit(status = 1, save = "no")
} else {
  cat("All checks passed: memory changes preserve results.\n")
}
