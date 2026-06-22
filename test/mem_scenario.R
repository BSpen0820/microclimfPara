# =============================================================================
# Memory measurement harness  (one scenario per fresh R process)
# =============================================================================
# Usage:
#   Rscript test/mem_scenario.R <pkgdir> setup        <fixturedir>
#   Rscript test/mem_scenario.R <pkgdir> <scenario>   <fixturedir> <outcsv> <label> [ncores]
#
# `pkgdir`  is the package source dir to load (HEAD for "after", a git worktree
#           at the parent commit for "before"); pkgload::load_all() compiles it.
# Peak memory is reported two ways:
#   * peak_r_mb  : R's own high-water allocation (gc max-used, reset per scenario)
#   * peak_rss_mb: process peak working set (ps::peak_wset) - true OS peak,
#                  captures the transient C++ output arrays that R may free
#                  before the next gc samples them.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
pkgdir     <- args[1]
scenario   <- args[2]
fixturedir <- args[3]

suppressMessages({
  library(pkgload)
  pkgload::load_all(pkgdir, quiet = TRUE, recompile = FALSE)
  library(terra)
})

data(climdata); data(vegp); data(soilc); data(dtmcaerth)
dtm <- rast(dtmcaerth)

# ---- memory helpers ---------------------------------------------------------
peak_reset <- function() invisible(gc(reset = TRUE, full = TRUE))
peak_r_mb  <- function() sum(gc(full = TRUE)[, 6])            # max used (Ncells+Vcells) Mb
rss_now_mb <- function() {
  mi <- tryCatch(ps::ps_memory_info(), error = function(e) NULL)
  if (!is.null(mi) && "wset" %in% names(mi)) return(unname(mi[["wset"]]) / 1024^2)
  NA_real_
}
rss_peak_mb <- function() {
  mi <- tryCatch(ps::ps_memory_info(), error = function(e) NULL)
  if (!is.null(mi) && "peak_wset" %in% names(mi)) return(unname(mi[["peak_wset"]]) / 1024^2)
  NA_real_
}

# =============================================================================
# Setup mode: build shared fixtures (point models + a snow model) ONCE.
# Point/snow outputs are identical across versions (validated to 1e-10), so the
# same fixtures are reused for before/after to keep the comparison fair & fast.
# =============================================================================
if (scenario == "setup") {
  dir.create(fixturedir, showWarnings = FALSE, recursive = TRUE)
  mp_ns <- runpointmodel(climdata, reqhgt = 0.05, dtm, vegp, soilc)        # full year (no subset)
  saveRDS(mp_ns, file.path(fixturedir, "micropoint_nosnow.rds"))
  mp_sn <- runpointmodel(climdata, reqhgt = 0.05, dtm, vegp, soilc)
  mp_sn <- subsetpointmodel(mp_sn, tstep = "month", what = "tmin")
  saveRDS(mp_sn, file.path(fixturedir, "micropoint_snow.rds"))
  smod  <- runsnowmodel(climdata, mp_sn, vegp, soilc, dtm, method = "fast")
  saveRDS(smod, file.path(fixturedir, "smod_fast.rds"))
  cat("fixtures written to", fixturedir, "\n")
  quit(save = "no")
}

outcsv <- args[4]
label  <- args[5]
ncores <- if (length(args) >= 6) as.integer(args[6]) else 2L

# =============================================================================
# Scenarios. All use parallel = TRUE: the output-gating fix lives in the
# parallel C++ runners, and the rm()/gc() edits sit on the same R paths.
# =============================================================================
run_scenario <- function() {
  if (scenario == "micro_tz_par") {
    # no-snow, full year, request Tz only -> gating fix: 1 vs 10 output arrays
    mp <- readRDS(file.path(fixturedir, "micropoint_nosnow.rds"))
    rss0 <- rss_now_mb(); peak_reset()
    el <- system.time(
      mout <- runmicro(mp, reqhgt = 0.05, vegp, soilc, dtm,
                       parallel = TRUE, ncores = ncores,
                       out = c(TRUE, rep(FALSE, 9)))
    )[["elapsed"]]
    pr <- peak_r_mb(); pk <- rss_peak_mb(); rm(mout)

  } else if (scenario == "micro_full_par") {
    # no-snow, full year, all 10 outputs -> gating is a no-op here; isolates rm()/gc()
    mp <- readRDS(file.path(fixturedir, "micropoint_nosnow.rds"))
    rss0 <- rss_now_mb(); peak_reset()
    el <- system.time(
      mout <- runmicro(mp, reqhgt = 0.05, vegp, soilc, dtm,
                       parallel = TRUE, ncores = ncores)
    )[["elapsed"]]
    pr <- peak_r_mb(); pk <- rss_peak_mb(); rm(mout)

  } else if (scenario == "runsnow_par") {
    # snow model (chunk-loop rm()/gc())
    mp <- readRDS(file.path(fixturedir, "micropoint_snow.rds"))
    rss0 <- rss_now_mb(); peak_reset()
    el <- system.time(
      smod <- runsnowmodel(climdata, mp, vegp, soilc, dtm,
                           method = "fast", parallel = TRUE, ncores = ncores)
    )[["elapsed"]]
    pr <- peak_r_mb(); pk <- rss_peak_mb(); rm(smod)

  } else if (scenario == "micro_snow_par") {
    # snow microclimate (gating in gridmicrosnow*Par + rm(snowin) in .runmicrosnow1)
    mp   <- readRDS(file.path(fixturedir, "micropoint_snow.rds"))
    smod <- readRDS(file.path(fixturedir, "smod_fast.rds"))
    rss0 <- rss_now_mb(); peak_reset()
    el <- system.time(
      mout <- runmicro(mp, reqhgt = 0.05, vegp, soilc, dtm,
                       snow = TRUE, snowmod = smod,
                       parallel = TRUE, ncores = ncores)
    )[["elapsed"]]
    pr <- peak_r_mb(); pk <- rss_peak_mb(); rm(mout)

  } else stop("unknown scenario: ", scenario)

  data.frame(version = label, scenario = scenario,
             peak_r_mb = round(pr, 1), peak_rss_mb = round(pk, 1),
             rss_pre_mb = round(rss0, 1), elapsed_s = round(el, 1))
}

row <- run_scenario()
write.table(row, outcsv, sep = ",", append = file.exists(outcsv),
            col.names = !file.exists(outcsv), row.names = FALSE)
cat(sprintf("%-6s %-15s  peakR=%7.1f Mb   peakRSS=%7.1f Mb   (%.0fs)\n",
            row$version, row$scenario, row$peak_r_mb, row$peak_rss_mb, row$elapsed_s))
