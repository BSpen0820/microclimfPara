# Validation_snowdaymov.R
# Synthetic correctness check for snowdaymov(): confirms it produces a
# continuously-varying signal (no day-boundary jump) where snowdayan()
# produces a block-constant signal with a sharp jump at every midnight.
# Run from the package root: Rscript test/Validation_snowdaymov.R

library(microclimfPara)

# 10 days of hourly data, alternating 0/10 degC each day - an extreme
# block signal designed to make snowdayan()'s day-boundary jump obvious.
ndays <- 10
tsteps <- ndays * 24
day_vals <- rep(c(0, 10), length.out = ndays)
Tg <- rep(day_vals, each = 24)
dim(Tg) <- c(1, 1, tsteps)

# snowdaymov() no longer computes a moving-average window (manCpp/nb) - it
# interpolates linearly between per-calendar-day means (each day's mean
# computed with NA hours skipped via na.rm semantics), anchored at hour 12
# of each day-block. meanD just needs to be a valid positive damping-depth
# estimate for that per-day-mean/fill/interpolate branch to run.
reqhgt <- -0.15
meanD <- matrix(0.35, nrow = 1, ncol = 1)

block  <- microclimfPara:::snowdayan(Tg)
moving <- microclimfPara:::snowdaymov(Tg, meanD, reqhgt)

block_v  <- as.vector(block)
moving_v <- as.vector(moving)

# Deltas across each midnight boundary (hour 24->25, 48->49, ...)
midnight_idx <- seq(24, tsteps - 24, by = 24)
block_jump  <- max(abs(diff(block_v)[midnight_idx]))
moving_jump <- max(abs(diff(moving_v)[midnight_idx]))

cat(sprintf("max midnight jump, snowdayan (old):  %.4f\n", block_jump))
cat(sprintf("max midnight jump, snowdaymov (new): %.4f\n", moving_jump))

stopifnot(block_jump > 5)   # old block behavior: near-full 10 degC swing at midnight
stopifnot(moving_jump < 5)  # new behavior: smoothed, well under the full swing

cat("\nPASS: snowdaymov() produces a smoothed midnight transition where snowdayan() jumps.\n")

## ---------------------------------------------------------------------
## Gappy/NA input test
## ---------------------------------------------------------------------
## This is the blind spot that mattered: an earlier moving-average design
## for this same regime-blend passed a clean synthetic test like the one
## above, then produced NaN on real gappy/NA-containing input. Exercise:
##   (a) a fully-NA day block (forward-fill/back-fill logic), and
##   (b) an isolated NA hour within an otherwise-valid day (na.rm day-mean
##       logic, i.e. the mean must exclude the NA hour rather than
##       propagate it).
## Both must produce zero NaN in the output, and (b) is checked against a
## hand-derived value computed independently in R.

Tg_gap <- as.vector(Tg)

# (a) Fully-NA day block: 0-based hours 48:71 (the 3rd 24-hour block) -> NA.
# 1-based array positions 49:72.
full_na_block_1based <- 49:72
Tg_gap[full_na_block_1based] <- NA_real_

# (b) Isolated NA hour: 0-based array position 144 (the first hour of the
# 24-hour block starting there), leaving the block's other 23 hours intact.
block_start_0based  <- (144 %/% 24) * 24   # = 144, already block-aligned
block_idx_1based    <- (block_start_0based + 1):(block_start_0based + 24)
isolated_na_1based   <- block_start_0based + 1  # 1-based index for position 144
stopifnot(isolated_na_1based %in% block_idx_1based)

day6_values_before <- Tg_gap[block_idx_1based]  # snapshot before introducing the NA
Tg_gap[isolated_na_1based] <- NA_real_

dim(Tg_gap) <- c(1, 1, tsteps)

moving_gap   <- microclimfPara:::snowdaymov(Tg_gap, meanD, reqhgt)
moving_gap_v <- as.vector(moving_gap)

stopifnot(!any(is.na(moving_gap_v)))
cat("PASS: snowdaymov() output contains no NA/NaN when given gappy input (fully-NA day + isolated NA hour).\n")

# --- (a) Fully-NA day block should be filled from a neighboring day's mean.
# Day 2 (0-based hours 24:47, 1-based positions 25:48) is untouched and is
# constant at 10 throughout, so its mean is exactly 10. Per snowdaymov()'s
# forward-fill-then-back-fill logic (src/microclimfCpp.cpp), the fully-NA
# day-block is forward-filled from the nearest earlier day with a valid
# mean - here, day 2 - so every hour of the fully-NA block should equal 10.
day2_mean <- mean(Tg_gap[25:48], na.rm = TRUE)
stopifnot(abs(day2_mean - 10) < 1e-9)
filled_block_vals <- moving_gap_v[full_na_block_1based]
stopifnot(all(abs(filled_block_vals - day2_mean) < 1e-9))
cat(sprintf("PASS: fully-NA day block filled at %.4f (previous day's mean), matches hand-derived value.\n", day2_mean))

# --- (b) Isolated NA hour: the day's own mean, computed with the NA hour
# excluded via na.rm, should appear at that day-block's anchor hour (hour
# 12, i.e. array position block_start + 12) since snowdaymov() anchors each
# day's own mean exactly at h == 12 (frac == 0 in both interpolation
# branches at that hour).
day6_mean_handderived <- mean(Tg_gap[block_idx_1based], na.rm = TRUE)
stopifnot(sum(!is.na(Tg_gap[block_idx_1based])) == 23)  # confirm exactly one NA excluded

anchor_1based <- block_start_0based + 12 + 1  # hour 12 of the block, 1-based
anchor_val <- moving_gap_v[anchor_1based]
stopifnot(abs(anchor_val - day6_mean_handderived) < 1e-9)
cat(sprintf("PASS: isolated-NA day's anchor-hour value (%.4f) matches mean(23 valid hours, na.rm=TRUE) (%.4f).\n",
            anchor_val, day6_mean_handderived))

cat("\nPASS: snowdaymov() gappy/NA-containing input test succeeded (fully-NA day block + isolated NA hour, no NaN, values match hand-derived expectations).\n")
