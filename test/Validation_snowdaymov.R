# Validation_snowdaymov.R
# Correctness check for the redesigned snowdaymov(): confirms it applies the
# same depth-dependent trailing rolling window as the no-snow model's own
# Tbelowgroundv() "complete" case (n = round(nb), manCpp(Tg, n) - see
# src/microclimfCpp.cpp), operating on a full-year, NA-free reference series
# built by .build_snow_Tgref() (R/internal.R, see
# test/Validation_build_snow_Tgref.R for that function's own tests), then
# slices the smoothed result down to just the requested snowdays hours.
# Run from the package root: Rscript test/Validation_snowdaymov.R

library(microclimfPara)

ndays <- 60
tsteps <- ndays * 24
set.seed(1)
# A synthetic full-year-like series: a slow seasonal sine wave (so a wide
# trailing window has something non-trivial to smooth) plus a small amount
# of hourly noise.
day_seq <- seq_len(ndays)
seasonal <- 5 + 10 * sin(2 * pi * day_seq / ndays)
Tgref_vec <- rep(seasonal, each = 24) + rnorm(tsteps, sd = 0.5)
Tgref <- array(Tgref_vec, dim = c(1, 1, tsteps))

# snowdays: every day from day 20 to day 40 (inclusive) - a contiguous
# mid-series block, so the output slice is easy to reason about by hand.
snowdays <- 20:40

reqhgt <- -0.15
meanD <- matrix(0.35, nrow = 1, ncol = 1)

Tzd <- microclimfPara:::snowdaymov(Tgref, snowdays, meanD, reqhgt)

# --- (1) Output shape: rows x cols x length(snowdays)*24
stopifnot(identical(dim(Tzd), c(1L, 1L, length(snowdays) * 24L)))
cat("PASS: output shape is rows x cols x length(snowdays)*24.\n")

# --- (2) Output matches manCpp(Tgref, n) directly, sliced to snowdays hours.
nb <- -118.35 * reqhgt / meanD[1, 1]
n <- round(nb)
cat(sprintf("nb = %.2f, n = %d (n > 48 uses manCpp's day-block-then-hourly path)\n", nb, n))
reference_full <- microclimfPara:::manCpp(as.vector(Tgref), n)
snow_hours <- as.vector(sapply(snowdays, function(d) ((d - 1) * 24 + 1):(d * 24)))
reference_sliced <- reference_full[snow_hours]
stopifnot(max(abs(as.vector(Tzd) - reference_sliced)) < 1e-9)
cat("PASS: snowdaymov() output matches manCpp(Tgref, n) sliced to snowdays hours (matches Tbelowgroundv()'s own complete-case formula).\n")

# --- (3) A masked pixel (all-NA Tgref) produces all-NA output, doesn't crash.
Tgref_masked <- array(NA_real_, dim = c(1, 1, tsteps))
Tzd_masked <- microclimfPara:::snowdaymov(Tgref_masked, snowdays, meanD, reqhgt)
stopifnot(all(is.na(as.vector(Tzd_masked))))
cat("PASS: fully-masked pixel (all-NA Tgref) produces all-NA output without error.\n")

# --- (4) A pixel with non-positive meanD is skipped (all-NA output for that pixel).
meanD_bad <- matrix(0, nrow = 1, ncol = 1)
Tzd_bad <- microclimfPara:::snowdaymov(Tgref, snowdays, meanD_bad, reqhgt)
stopifnot(all(is.na(as.vector(Tzd_bad))))
cat("PASS: non-positive meanD produces all-NA output for that pixel.\n")

# --- (5) n >= tsteps falls back to a flat pixel-mean (matches Tbelowgroundv()'s
# own guard for this case) rather than relying on manCpp's degenerate behavior
# at n close to/exceeding the series length.
short_tsteps <- 48  # 2 days
Tgref_short <- array(rep(c(0, 10), each = 24), dim = c(1, 1, short_tsteps))
snowdays_short <- 1:2
meanD_tiny <- matrix(0.001, nrow = 1, ncol = 1)  # huge nb -> n >= tsteps
Tzd_short <- microclimfPara:::snowdaymov(Tgref_short, snowdays_short, meanD_tiny, reqhgt)
expected_flat <- mean(as.vector(Tgref_short))
stopifnot(max(abs(as.vector(Tzd_short) - expected_flat)) < 1e-9)
cat(sprintf("PASS: n >= tsteps falls back to the flat pixel mean (%.4f), matching Tbelowgroundv()'s own guard.\n", expected_flat))

cat("\nPASS: all snowdaymov() correctness checks succeeded.\n")
