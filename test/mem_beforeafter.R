# Peak-memory before/after harness for the runmicro snow path.
# Measures R's high-water mark (captures R intermediates AND C++ NumericVector
# outputs, both of which use R's allocator) via gc(reset=TRUE).
# A single cached smod (unaffected by the changes) is reused so the comparison
# isolates runmicro(snow=TRUE), which exercises Fix 1 (.runmodel2Cpp intermediate
# freeing) and Fix 2 (.runmicrosnow2 in-place merge).
suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages(library(terra)); suppressMessages(library(digest))

label <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(label)) label <- "RUN"
outspec <- commandArgs(trailingOnly = TRUE)[2]
out_v <- if (!is.na(outspec) && outspec == "tz") c(TRUE, rep(FALSE, 9)) else rep(TRUE, 10)

data(climdata); data(vegp); data(soilc); data(dtmcaerth)
.rast <- function(m, tem) { r <- rast(m); ext(r) <- ext(tem); crs(r) <- crs(tem); r }
.ta   <- function(x, dtm, xdim = 5, ydim = 5) {
  a <- array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x))); .rast(a, dtm)
}
dtm <- rast(dtmcaerth)
climdata$temp <- climdata$temp - 12
climarrayr <- list(
  temp = .ta(climdata$temp, dtm), relhum = .ta(climdata$relhum, dtm),
  pres = .ta(climdata$pres, dtm), swdown = .ta(climdata$swdown, dtm),
  difrad = .ta(climdata$difrad, dtm), lwdown = .ta(climdata$lwdown, dtm),
  windspeed = .ta(climdata$windspeed, dtm), winddir = .ta(climdata$winddir, dtm),
  precip = .ta(climdata$precip, dtm))
tme  <- as.POSIXlt(climdata$obs_time, tz = "UTC")
dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)
micropointa <- readRDS("Validation_test/base/micropointa.rds")

# Reuse one smod across both runs (runsnowmodel is unaffected by the changes)
smod_cache <- "test/.smod_cache.rds"
if (file.exists(smod_cache)) {
  smod_fast <- readRDS(smod_cache)
} else {
  smod_fast <- runsnowmodel(climarrayr, micropointa, vegp, soilc, dtm, dtmc,
                            tme = tme, altcorrect = 0, method = "fast")
  saveRDS(smod_fast, smod_cache)
}

use_snow <- !(!is.na(outspec) && outspec == "nosnow")
invisible(gc(reset = TRUE, full = TRUE))
if (use_snow) {
  mout <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                   altcorrect = 0, snow = TRUE, snowmod = smod_fast, out = out_v)
} else {
  mout <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                   altcorrect = 0, snow = FALSE)
}
g <- gc(full = TRUE)
peak_vec <- g["Vcells", ncol(g)]         # R gc high-water mark (gc-timing dependent)
# True OS peak RSS for the whole process lifetime (Windows tracks this natively,
# independent of gc timing) -- the analogue of the cluster's MaxRSS.
pid <- Sys.getpid()
pk  <- system(sprintf(
  'powershell -NoProfile -Command "(Get-Process -Id %d).PeakWorkingSet64"', pid),
  intern = TRUE)
peak_rss_mb <- as.numeric(pk) / 1024 / 1024
cat(sprintf("\n[%s] PEAK_RSS_Mb=%.1f  gc_vcells_Mb=%.1f  hash=%s\n",
            label, peak_rss_mb, peak_vec, digest(mout, "md5")))
