# 05_snowmicro_par_validate.R
# Validates gridmicrosnow1Par / gridmicrosnow2Par by comparing parallel output
# against the existing serial mout_snow result.
# Loads one reference at a time; rm()/gc() after each comparison to keep memory lean.

library(microclimfPara)
library(terra)
library(digest)

cat("microclimfPara version:", as.character(packageVersion("microclimfPara")), "\n")
cat("Loaded at:", format(Sys.time()), "\n\n")

# ── Minimal shared inputs ─────────────────────────────────────────────────────
data(vegp); data(soilc); data(dtmcaerth)
dtm  <- rast(dtmcaerth)
dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)

results <- list()

# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 1: mout_snow (snow-day microclimate, parallel vs serial)
# Uses existing micropointa_snow and smod_snow — no need to re-run those models.
# ─────────────────────────────────────────────────────────────────────────────
cat("=== Block 1: mout_snow parallel vs serial ===\n")

micropointa_snow <- readRDS("test/MicroPar_Ser/micropointa_snow.rds")
smod_snow        <- readRDS("test/MicroPar_Ser/smod_snow.rds")

cat("  Running parallel mout_snow...\n")
t1 <- proc.time()
par_mout_snow <- runmicro(micropointa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                           altcorrect = 0, parallel = TRUE, ncores = 2,
                           snow = TRUE, snowmod = smod_snow)
cat("  Done in", round((proc.time() - t1)[3], 1), "s\n")

rm(micropointa_snow, smod_snow); gc()

# Compare against serial result
ref <- readRDS("test/MicroPar_Ser/mout_snow.rds")
ae  <- all.equal(ref, par_mout_snow, tolerance = 1e-10)
pass <- isTRUE(ae)
detail <- if (pass) "" else paste(head(ae, 3), collapse = "; ")
cat(sprintf("  %-40s | %s\n", "mout_snow ser vs par", if (pass) "PASS" else paste("FAIL:", detail)))
results[["mout_snow_ser_vs_par"]] <- data.frame(label = "mout_snow ser vs par",
                                                 pass = pass, detail = detail,
                                                 stringsAsFactors = FALSE)
rm(ref); gc()

cat("  par hash:", digest(par_mout_snow, algo = "md5"), "\n")
rm(par_mout_snow); gc()

# ── Summary ───────────────────────────────────────────────────────────────────
summary_df <- do.call(rbind, results)
dir.create("ClaudeTest_Validation/revalidation_results", showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, "ClaudeTest_Validation/revalidation_results/05_summary.csv", row.names = FALSE)

cat("\n=== Final Summary ===\n")
print(summary_df[, c("label", "pass", "detail")])

all_pass <- all(summary_df$pass)
cat("\n", if (all_pass) "ALL SNOW MICROCLIMATE PARALLEL OUTPUTS MATCH SERIAL — IMPLEMENTATION CONFIRMED"
         else "SNOW MICROCLIMATE PARALLEL OUTPUTS DIFFER — INVESTIGATION NEEDED", "\n\n")
