# 02_ncores1_test.R
#
# PURPOSE: Diagnostic script — only run if 01_quantify_differences.R shows
# LARGE differences (> 1e-3). Tests parallel with ncores=1 and a 1x1 grid
# to isolate whether the bug is in per-cell computation or TBB scheduling.
#
# Run from: D:/Code/PhD/R_Packages/microclimfPar/
# setwd("D:/Code/PhD/R_Packages/microclimfPar/")
#
# Requires microclimfPara to be installed (serial build is fine for this test).

library(microclimfPara)
library(terra)
library(digest)

data(climdata)
data(vegp)
data(soilc)
data(dtmcaerth)

dtm <- rast(dtmcaerth)

.ta <- function(x, dtm, xdim = 5, ydim = 5) {
    a <- array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x)))
    r <- rast(a)
    ext(r) <- ext(dtm)
    crs(r) <- crs(dtm)
    r
}

climarrayr <- list(
    temp      = .ta(climdata$temp, dtm),
    relhum    = .ta(climdata$relhum, dtm),
    pres      = .ta(climdata$pres, dtm),
    swdown    = .ta(climdata$swdown, dtm),
    difrad    = .ta(climdata$difrad, dtm),
    lwdown    = .ta(climdata$lwdown, dtm),
    windspeed = .ta(climdata$windspeed, dtm),
    winddir   = .ta(climdata$winddir, dtm),
    precip    = .ta(climdata$precip, dtm)
)

tme   <- as.POSIXlt(climdata$obs_time, tz = "UTC")
dtmc  <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)

dir.create("./ClaudeTest_Validation/ncores1_results", showWarnings = FALSE)

cat("Running point model...\n")
micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)

# ---- Test 1: ncores=1 parallel vs serial ----
cat("Running serial grid model...\n")
mout_serial   <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                          altcorrect = 0, parallel = FALSE)

cat("Running parallel ncores=1...\n")
mout_par1     <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                          altcorrect = 0, parallel = TRUE, ncores = 1)

cat("Running parallel ncores=4...\n")
mout_par4     <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                          altcorrect = 0, parallel = TRUE, ncores = 4)

saveRDS(mout_serial, "./ClaudeTest_Validation/ncores1_results/mout_serial.rds",   compress = FALSE)
saveRDS(mout_par1,   "./ClaudeTest_Validation/ncores1_results/mout_par1.rds",     compress = FALSE)
saveRDS(mout_par4,   "./ClaudeTest_Validation/ncores1_results/mout_par4.rds",     compress = FALSE)

cat("\n=== ncores=1 test results ===\n")
cat("serial identical to par(ncores=1):", identical(mout_serial, mout_par1), "\n")
cat("serial identical to par(ncores=4):", identical(mout_serial, mout_par4), "\n")
cat("par(ncores=1) identical to par(ncores=4):", identical(mout_par1, mout_par4), "\n")

# If serial != par(ncores=1), the issue is per-cell computation (compiler codegen)
# If serial == par(ncores=1) != par(ncores=4), the issue is TBB scheduling / race

# ---- Test 2: 1x1 grid ----
cat("\n--- 1x1 grid test ---\n")
dtm1    <- dtm[1, 1, drop = FALSE]
dtmc1   <- aggregate(dtm1, 1, fun = "mean", na.rm = TRUE)

# Build 1x1 climate array (reuse same time series)
.ta1 <- function(x, dtm) {
    a <- array(x, dim = c(1, 1, length(x)))
    r <- rast(a)
    ext(r) <- ext(dtm)
    crs(r) <- crs(dtm)
    r
}
climarrayr1 <- lapply(names(climarrayr), function(nm) .ta1(climdata[[nm]], dtm1))
names(climarrayr1) <- names(climarrayr)
# Fix field name mismatches if needed
if ("swdown" %in% names(climdata)) climarrayr1$swdown <- .ta1(climdata$swdown, dtm1)

micropointa1 <- runpointmodela(climarrayr1, tme, reqhgt = 0.05, dtm1, vegp, soilc)
mout1_ser <- runmicro(micropointa1, reqhgt = 0.05, vegp, soilc, dtm1, dtmc1,
                      altcorrect = 0, parallel = FALSE)
mout1_par <- runmicro(micropointa1, reqhgt = 0.05, vegp, soilc, dtm1, dtmc1,
                      altcorrect = 0, parallel = TRUE, ncores = 1)

cat("1x1 serial vs parallel(ncores=1) identical:", identical(mout1_ser, mout1_par), "\n")

ae <- all.equal(mout1_ser, mout1_par, tolerance = 0)
if (!isTRUE(ae)) {
    cat("Differences in 1x1 grid:\n")
    for (m in head(ae, 10)) cat("  ", m, "\n")
}
