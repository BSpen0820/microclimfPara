# 03_revalidate_after_fix.R
#
# PURPOSE: Full revalidation after rebuilding the package with -ffp-contract=off.
# Compares base (original microclimf), serial (microclimfPara no-parallel),
# and parallel (microclimfPara parallel=TRUE) using both identical() and all.equal().
#
# Prerequisites:
#   1. Rebuilt microclimfPara with updated Makevars (devtools::install() or
#      devtools::clean_dll(); devtools::load_all())
#   2. Original test/ RDS files from test/Base_Results/ are present
#
# Run from: D:/Code/PhD/R_Packages/microclimfPar/
# setwd("D:/Code/PhD/R_Packages/microclimfPar/")

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

climdata_snow <- climdata
climdata_snow$temp <- climdata$temp - 12

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

climarrayr_snow <- list(
    temp      = .ta(climdata_snow$temp, dtm),
    relhum    = .ta(climdata_snow$relhum, dtm),
    pres      = .ta(climdata_snow$pres, dtm),
    swdown    = .ta(climdata_snow$swdown, dtm),
    difrad    = .ta(climdata_snow$difrad, dtm),
    lwdown    = .ta(climdata_snow$lwdown, dtm),
    windspeed = .ta(climdata_snow$windspeed, dtm),
    winddir   = .ta(climdata_snow$winddir, dtm),
    precip    = .ta(climdata_snow$precip, dtm)
)

tme  <- as.POSIXlt(climdata$obs_time, tz = "UTC")
dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)

dir.create("./ClaudeTest_Validation/revalidation_results", showWarnings = FALSE)

cat("Running microclimfPara models (post-fix)...\n\n")

# ---- No-snow ----
cat("[1/5] Point model (no snow)...\n")
micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa, "./ClaudeTest_Validation/revalidation_results/micropointa_nosnow.rds",
        compress = FALSE)

cat("[2/5] Grid model serial (no snow)...\n")
mout_ser <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                     altcorrect = 0, parallel = FALSE)
saveRDS(mout_ser, "./ClaudeTest_Validation/revalidation_results/mout_nosnow_serial.rds",
        compress = FALSE)

cat("[3/5] Grid model parallel (no snow)...\n")
mout_par <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                     altcorrect = 0, parallel = TRUE, ncores = 4)
saveRDS(mout_par, "./ClaudeTest_Validation/revalidation_results/mout_nosnow_parallel.rds",
        compress = FALSE)

# ---- Snow ----
cat("[4/5] Snow model + grid model (snow)...\n")
micropointa_snow <- runpointmodela(climarrayr_snow, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa_snow,
        "./ClaudeTest_Validation/revalidation_results/micropointa_snow.rds",
        compress = FALSE)

smod_par <- runsnowmodel(climarrayr_snow, micropointa_snow, vegp, soilc, dtm, dtmc,
                         tme = tme, altcorrect = 0, method = "slow",
                         parallel = TRUE, ncores = 4)
saveRDS(smod_par, "./ClaudeTest_Validation/revalidation_results/smod_snow_parallel.rds",
        compress = FALSE)

mout_snow_par <- runmicro(micropointa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                          altcorrect = 0, parallel = TRUE, ncores = 4,
                          snow = TRUE, snowmod = smod_par)
saveRDS(mout_snow_par,
        "./ClaudeTest_Validation/revalidation_results/mout_snow_parallel.rds",
        compress = FALSE)

cat("[5/5] Comparing against base results...\n\n")

# ---- Comparison ----
compare <- function(nm_base, obj_new, label) {
    base <- readRDS(file.path("./test/Base_Results", paste0(nm_base, ".rds")))
    is_id <- identical(base, obj_new)
    ae    <- all.equal(base, obj_new, tolerance = 1e-10)
    ae_strict <- all.equal(base, obj_new, tolerance = 0)

    cat("=====", label, "=====\n")
    cat("  identical():              ", is_id, "\n")
    cat("  all.equal(tol=1e-10):     ", isTRUE(ae), "\n")
    cat("  all.equal(tol=0):         ", isTRUE(ae_strict), "\n")

    if (!isTRUE(ae_strict) && !is_id) {
        msgs <- if (is.character(ae_strict)) ae_strict else character(0)
        for (m in head(msgs, 5)) cat("  diff:", m, "\n")
    }

    if (is_id) {
        cat("  RESULT: PASS (bit-exact match)\n")
    } else if (isTRUE(ae)) {
        cat("  RESULT: PASS (within 1e-10 tolerance — acceptable numerical equivalence)\n")
    } else {
        cat("  RESULT: FAIL — differences exceed 1e-10\n")
    }
    cat("\n")

    list(name = label, identical = is_id, pass_1e10 = isTRUE(ae))
}

results <- list(
    compare("micropointa_nosnow", micropointa,       "micropointa_nosnow"),
    compare("mout_nosnow",        mout_ser,           "mout_nosnow (serial)"),
    compare("mout_nosnow",        mout_par,           "mout_nosnow (parallel)"),
    compare("micropointa_snow",   micropointa_snow,   "micropointa_snow"),
    compare("smod_snow",          smod_par,           "smod_snow (parallel)"),
    compare("mout_snow",          mout_snow_par,      "mout_snow (parallel)")
)

# Summary table
cat("============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")
cat(sprintf("%-35s %-10s %-10s\n", "Output", "identical", "pass_1e10"))
cat(strrep("-", 58), "\n")
for (r in results) {
    cat(sprintf("%-35s %-10s %-10s\n",
                r$name,
                as.character(r$identical),
                as.character(r$pass_1e10)))
}

# Save summary
summary_df <- do.call(rbind, lapply(results, as.data.frame))
write.csv(summary_df,
          "./ClaudeTest_Validation/revalidation_results/summary.csv",
          row.names = FALSE)
cat("\nSummary written to ClaudeTest_Validation/revalidation_results/summary.csv\n")
