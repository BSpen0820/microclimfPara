devtools::load_all()
library(terra)
library(digest)

# ============================================================================ #
# Validation: regression test — default behaviour must match baseline hashes   #
#                                                                               #
# Any change to runsnowmodel() or runmicro() must produce outputs with hashes  #
# matching the table below. This is a hard requirement.                         #
#                                                                               #
# NOTE: Terrain params (slr, apr, hor, wsa, svf) were NOT added to             #
# runsnowmodel() because the active dispatch path for a full-year dataset is   #
# .snowmodel2, which recomputes terrain per iteration from the snow-adjusted   #
# DTM (dtms). Any fixed terrain passed externally would differ from the        #
# per-iteration computation and produce non-matching hashes.                   #
# ============================================================================ #

baseline_hashes <- c(
  micropointa    = "9375c48eb328e54b65182db0dff686e7",
  smod_fast      = "924717ec3f7a08924aa037617ac60116",
  mout_fast_snow = "f8e2d10ac5efe38ef4952aa348dc1a9c"
)

check_hash <- function(obj, name) {
  h <- digest(obj, algo = "md5")
  expected <- baseline_hashes[name]
  status <- if (!is.na(expected) && h == expected) "PASS" else
            if (is.na(expected)) "INFO (no baseline)" else "FAIL <-- outputs changed"
  cat(sprintf("  %-22s  %s  [%s]\n", name, h, status))
  invisible(h == expected || is.na(expected))
}

# ============================================================================ #
# Data prep — mirrors Validation_base.R exactly                                 #
# ============================================================================ #
data(climdata); data(vegp); data(soilc); data(dtmcaerth)

.rast <- function(m, tem) { r <- rast(m); ext(r) <- ext(tem); crs(r) <- crs(tem); r }
.ta   <- function(x, dtm, xdim = 5, ydim = 5) {
  a <- array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x)))
  .rast(a, dtm)
}

dtm <- rast(dtmcaerth)
climdata$temp <- climdata$temp - 12

climarrayr <- list(
  temp      = .ta(climdata$temp,      dtm),
  relhum    = .ta(climdata$relhum,    dtm),
  pres      = .ta(climdata$pres,      dtm),
  swdown    = .ta(climdata$swdown,    dtm),
  difrad    = .ta(climdata$difrad,    dtm),
  lwdown    = .ta(climdata$lwdown,    dtm),
  windspeed = .ta(climdata$windspeed, dtm),
  winddir   = .ta(climdata$winddir,   dtm),
  precip    = .ta(climdata$precip,    dtm))

tme  <- as.POSIXlt(climdata$obs_time, tz = "UTC")
dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)

micropointa <- readRDS("Validation_test/base/micropointa.rds")

# ============================================================================ #
# Regression test                                                                #
# ============================================================================ #
cat("\n=== Regression test (must match baseline hashes) ===\n")

smod_fast <- runsnowmodel(climarrayr, micropointa, vegp, soilc, dtm, dtmc,
                          tme = tme, altcorrect = 0, method = "fast")

mout_fast_snow <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc,
                           altcorrect = 0, snow = TRUE, snowmod = smod_fast)

cat("\n")
p1 <- check_hash(micropointa,   "micropointa")
p2 <- check_hash(smod_fast,     "smod_fast")
p3 <- check_hash(mout_fast_snow,"mout_fast_snow")

if (all(c(p1, p2, p3))) {
  cat("\n=== VALIDATION PASSED: all hashes match baseline ===\n")
} else {
  cat("\n=== VALIDATION FAILED: one or more hashes differ from baseline ===\n")
}
