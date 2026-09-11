# 04_parallel_only_revalidate.R
# Lean revalidation: load existing base/serial results one at a time, run only
# the parallel models, and compare immediately before freeing memory.
# Does NOT re-run the serial grid model (avoids OOM).
# Run this AFTER rebuilding the package with -ffp-contract=off.

library(microclimfPara)
library(terra)
library(digest)

cat("microclimfPara version:", as.character(packageVersion("microclimfPara")), "\n")
cat("Loaded at:", format(Sys.time()), "\n\n")

# ── Helpers ──────────────────────────────────────────────────────────────────
.rast <- function(m, tem) { r <- rast(m); ext(r) <- ext(tem); crs(r) <- crs(tem); r }
.ta   <- function(x, dtm, xdim = 5, ydim = 5) {
  a <- array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x)))
  .rast(a, dtm)
}

compare_pair <- function(label, ref, par, tol = 1e-10) {
  ae <- all.equal(ref, par, tolerance = tol)
  id <- isTRUE(ae)
  detail <- if (id) "" else paste(head(ae, 2), collapse = "; ")
  cat(sprintf("  %-32s | %s\n", label, if (id) "PASS" else paste("FAIL:", detail)))
  data.frame(label = label, pass = id, detail = detail, stringsAsFactors = FALSE)
}

results <- list()

# ── Minimal data needed throughout ───────────────────────────────────────────
data(vegp); data(soilc); data(dtmcaerth); data(climdata)
dtm  <- rast(dtmcaerth)
dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)
tme  <- as.POSIXlt(climdata$obs_time, tz = "UTC")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 1: mout_nosnow
# runmicro() only needs micropointa — no climate array required
# ─────────────────────────────────────────────────────────────────────────────
cat("=== Block 1: mout_nosnow ===\n")

micropointa_nosnow <- readRDS("test/MicroPar_Ser/micropointa_nosnow.rds")

cat("  Running parallel mout_nosnow...\n")
t1 <- proc.time()
par_mout_nosnow <- runmicro(micropointa_nosnow, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                             altcorrect = 0, parallel = TRUE, ncores = 2)
cat("  Done in", round((proc.time() - t1)[3], 1), "s\n")

rm(micropointa_nosnow); gc()

ref <- readRDS("test/Base_Results/mout_nosnow.rds")
results[["mout_nosnow_base"]] <- compare_pair("mout_nosnow base vs par", ref, par_mout_nosnow)
rm(ref); gc()

ref <- readRDS("test/MicroPar_Ser/mout_nosnow.rds")
results[["mout_nosnow_ser"]]  <- compare_pair("mout_nosnow ser  vs par", ref, par_mout_nosnow)
rm(ref); gc()

cat("  hash:", digest(par_mout_nosnow, algo = "md5"), "\n")
rm(par_mout_nosnow); gc()

# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 2: smod_snow
# runsnowmodel() needs climarrayr_snow — build it here, free after use
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Block 2: smod_snow ===\n")

climdata_snow       <- climdata
climdata_snow$temp  <- climdata$temp - 12

climarrayr_snow <- list(
  temp      = .ta(climdata_snow$temp,      dtm),
  relhum    = .ta(climdata_snow$relhum,    dtm),
  pres      = .ta(climdata_snow$pres,      dtm),
  swdown    = .ta(climdata_snow$swdown,    dtm),
  difrad    = .ta(climdata_snow$difrad,    dtm),
  lwdown    = .ta(climdata_snow$lwdown,    dtm),
  windspeed = .ta(climdata_snow$windspeed, dtm),
  winddir   = .ta(climdata_snow$winddir,   dtm),
  precip    = .ta(climdata_snow$precip,    dtm)
)
rm(climdata_snow, climdata); gc()

micropointa_snow <- readRDS("test/MicroPar_Ser/micropointa_snow.rds")

cat("  Running parallel smod_snow...\n")
t1 <- proc.time()
par_smod_snow <- runsnowmodel(climarrayr_snow, micropointa_snow, vegp, soilc, dtm, dtmc,
                               tme = tme, altcorrect = 0, method = "slow",
                               parallel = TRUE, ncores = 2)
cat("  Done in", round((proc.time() - t1)[3], 1), "s\n")

rm(climarrayr_snow); gc()

ref <- readRDS("test/Base_Results/smod_snow.rds")
results[["smod_snow_base"]] <- compare_pair("smod_snow   base vs par", ref, par_smod_snow)
rm(ref); gc()

ref <- readRDS("test/MicroPar_Ser/smod_snow.rds")
results[["smod_snow_ser"]]  <- compare_pair("smod_snow   ser  vs par", ref, par_smod_snow)
rm(ref); gc()

cat("  hash:", digest(par_smod_snow, algo = "md5"), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 3: mout_snow
# needs micropointa_snow (still in memory) and par_smod_snow (still in memory)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Block 3: mout_snow ===\n")

cat("  Running parallel mout_snow...\n")
t1 <- proc.time()
par_mout_snow <- runmicro(micropointa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                           altcorrect = 0, parallel = TRUE, ncores = 2,
                           snow = TRUE, snowmod = par_smod_snow)
cat("  Done in", round((proc.time() - t1)[3], 1), "s\n")

rm(micropointa_snow, par_smod_snow); gc()

ref <- readRDS("test/Base_Results/mout_snow.rds")
results[["mout_snow_base"]] <- compare_pair("mout_snow   base vs par", ref, par_mout_snow)
rm(ref); gc()

ref <- readRDS("test/MicroPar_Ser/mout_snow.rds")
results[["mout_snow_ser"]]  <- compare_pair("mout_snow   ser  vs par", ref, par_mout_snow)
rm(ref); gc()

cat("  hash:", digest(par_mout_snow, algo = "md5"), "\n")
rm(par_mout_snow); gc()

# ── Summary ───────────────────────────────────────────────────────────────────
summary_df <- do.call(rbind, results)
dir.create("ClaudeTest_Validation/revalidation_results", showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, "ClaudeTest_Validation/revalidation_results/summary.csv", row.names = FALSE)

cat("\n=== Final Summary ===\n")
print(summary_df[, c("label", "pass", "detail")])

base_rows <- grep("base vs par", summary_df$label)
all_pass  <- all(summary_df$pass[base_rows])
cat("\n", if (all_pass) "ALL PARALLEL OUTPUTS MATCH BASE — FIX CONFIRMED"
         else "SOME PARALLEL OUTPUTS STILL DIFFER — FURTHER INVESTIGATION NEEDED", "\n\n")

cat("Summary written to ClaudeTest_Validation/revalidation_results/summary.csv\n")
