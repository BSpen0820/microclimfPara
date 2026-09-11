# 06_pointmodela_par_validate.R
# Validates runpointmodela(parallel=TRUE) by comparing against existing serial
# reference outputs. Loads one reference at a time; rm()/gc() after each block.

library(microclimfPara)
library(terra)
library(digest)

cat("microclimfPara version:", as.character(packageVersion("microclimfPara")), "\n")
cat("Loaded at:", format(Sys.time()), "\n\n")

.ta <- function(x, dtm, xdim = 5, ydim = 5) {
  a <- array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x)))
  r <- rast(a)
  ext(r) <- ext(dtm)
  crs(r) <- crs(dtm)
  r
}

data(vegp); data(soilc); data(dtmcaerth); data(climdata)
dtm <- rast(dtmcaerth)
tme <- as.POSIXlt(climdata$obs_time, tz = "UTC")

results <- list()

# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 1: no-snow case
# ─────────────────────────────────────────────────────────────────────────────
cat("=== Block 1: runpointmodela no-snow (parallel vs serial) ===\n")

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

cat("  Running parallel runpointmodela (no-snow)...\n")
t1 <- proc.time()
par_micropointa_nosnow <- runpointmodela(climarrayr, tme, reqhgt = 0.05,
                                          dtm, vegp, soilc,
                                          parallel = TRUE, ncores = 2)
cat("  Done in", round((proc.time() - t1)[3], 1), "s\n")

rm(climarrayr); gc()

ref <- readRDS("test/MicroPar_Ser/micropointa_nosnow.rds")
ae  <- all.equal(ref, par_micropointa_nosnow, tolerance = 1e-10)
pass <- isTRUE(ae)
detail <- if (pass) "" else paste(head(ae, 3), collapse = "; ")
cat(sprintf("  %-44s | %s\n", "micropointa_nosnow ser vs par",
            if (pass) "PASS" else paste("FAIL:", detail)))
results[["nosnow"]] <- data.frame(label = "micropointa_nosnow ser vs par",
                                   pass = pass, detail = detail,
                                   stringsAsFactors = FALSE)
rm(ref); gc()

cat("  par hash:", digest(par_micropointa_nosnow, algo = "md5"), "\n")
rm(par_micropointa_nosnow); gc()

# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 2: snow case (temp - 12)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Block 2: runpointmodela snow (parallel vs serial) ===\n")

climdata_snow      <- climdata
climdata_snow$temp <- climdata$temp - 12

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
rm(climdata_snow); gc()

cat("  Running parallel runpointmodela (snow)...\n")
t1 <- proc.time()
par_micropointa_snow <- runpointmodela(climarrayr_snow, tme, reqhgt = 0.05,
                                        dtm, vegp, soilc,
                                        parallel = TRUE, ncores = 2)
cat("  Done in", round((proc.time() - t1)[3], 1), "s\n")

rm(climarrayr_snow); gc()

ref <- readRDS("test/MicroPar_Ser/micropointa_snow.rds")
ae  <- all.equal(ref, par_micropointa_snow, tolerance = 1e-10)
pass <- isTRUE(ae)
detail <- if (pass) "" else paste(head(ae, 3), collapse = "; ")
cat(sprintf("  %-44s | %s\n", "micropointa_snow ser vs par",
            if (pass) "PASS" else paste("FAIL:", detail)))
results[["snow"]] <- data.frame(label = "micropointa_snow ser vs par",
                                 pass = pass, detail = detail,
                                 stringsAsFactors = FALSE)
rm(ref); gc()

cat("  par hash:", digest(par_micropointa_snow, algo = "md5"), "\n")
rm(par_micropointa_snow); gc()

# ── Summary ───────────────────────────────────────────────────────────────────
summary_df <- do.call(rbind, results)
dir.create("ClaudeTest_Validation/revalidation_results", showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, "ClaudeTest_Validation/revalidation_results/06_summary.csv", row.names = FALSE)

cat("\n=== Final Summary ===\n")
print(summary_df[, c("label", "pass", "detail")])

all_pass <- all(summary_df$pass)
cat("\nOverall:", if (all_pass) "ALL PASS" else "FAILURES PRESENT", "\n")
if (!all_pass) stop("Validation failed — see detail above.")
