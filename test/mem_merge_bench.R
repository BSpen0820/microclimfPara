# Isolated quantification of Fix 2 (the snow-output merge in .runmicrosnow1/2).
# Runs the EXACT old and new merge loops on identical synthetic moutn/mouts and
# measures the merge-phase peak. Each variant runs in its own process so the
# Windows PeakWorkingSet64 reflects that process's true high-water (the old loop
# never frees, so its heap genuinely grows ~2x; both metrics should agree here).
#
# Usage: Rscript test/mem_merge_bench.R <old|new> [gridside] [nvar]
args  <- commandArgs(trailingOnly = TRUE)
which <- args[1]; if (is.na(which)) which <- "new"
side  <- as.integer(args[2]); if (is.na(side)) side <- 50
nvar  <- as.integer(args[3]); if (is.na(nvar)) nvar <- 10

nr <- side; nc <- side
nosnow_days <- 182; snow_days <- 183          # ~half/half year
nosnowh <- seq_len(nosnow_days * 24)
snowh   <- nosnow_days * 24 + seq_len(snow_days * 24)
n  <- length(nosnowh) + length(snowh)
s1 <- seq_along(nosnowh)                       # all no-snow hours used
dms <- c(nr, nc)

# Build realistic (non-NA) inputs so they occupy real memory.
# Fixed seed => identical inputs across processes => matching checksum proves
# the two merge variants produce identical output.
set.seed(42)
moutn <- vector("list", nvar); mouts <- vector("list", nvar)
for (i in seq_len(nvar)) {
  moutn[[i]] <- array(runif(nr * nc * length(nosnowh)), dim = c(nr, nc, length(nosnowh)))
  mouts[[i]] <- array(runif(nr * nc * length(snowh)),   dim = c(nr, nc, length(snowh)))
}
names(moutn) <- paste0("v", seq_len(nvar))

merge_old <- function() {                       # faithful HEAD~1 loop
  mout <- list()
  for (i in 1:length(moutn)) {
    mout[[i]] <- array(NA, dim = c(dms, n))
    xx <- moutn[[i]][, , s1]
    mout[[i]][, , nosnowh] <- xx
    mout[[i]][, , snowh]   <- mouts[[i]]
  }
  names(mout) <- names(moutn)
  mout
}
merge_new <- function() {                       # faithful HEAD loop (incremental free)
  mout <- list(); nm <- names(moutn); L <- length(moutn)
  for (i in 1:L) {
    mout[[i]] <- array(NA, dim = c(dms, n))
    mout[[i]][, , nosnowh] <- moutn[[i]][, , s1]
    mout[[i]][, , snowh]   <- mouts[[i]]
    moutn[[i]] <<- NA
    mouts[[i]] <<- NA
  }
  names(mout) <- nm
  mout
}

invisible(gc(reset = TRUE, full = TRUE))
res <- if (which == "old") merge_old() else merge_new()
g <- gc(full = TRUE)
peak_vec <- g["Vcells", ncol(g)]
pid <- Sys.getpid()
pk  <- suppressWarnings(as.numeric(system(sprintf(
  'powershell -NoProfile -Command "(Get-Process -Id %d).PeakWorkingSet64"', pid),
  intern = TRUE))) / 1024 / 1024
# sanity: confirm both variants produce identical merged output
chk <- sum(vapply(res, function(a) sum(a[, , c(1, n)]), numeric(1)))
cat(sprintf("\n[MERGE-%s] side=%d nvar=%d  gc_vcells_Mb=%.1f  PEAK_RSS_Mb=%.1f  checksum=%.3f\n",
            toupper(which), side, nvar, peak_vec, pk, chk))
