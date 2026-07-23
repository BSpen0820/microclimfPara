# Validation_build_snow_Tgref.R
# Correctness check for .build_snow_Tgref() (R/internal.R): the R-level
# assembly step that builds a full-year, NA-free reference series for
# snowdaymov() to smooth. Splices Tg (snow-surface skin temperature) onto
# snowdays/mixed-day hours and moutn$Tz (the no-snow model's own below-ground
# prediction) onto pure-nosnow-day hours, then fills any remaining NA hour
# by linear interpolation from that pixel's own nearest valid hours.
# Run from the package root: Rscript test/Validation_build_snow_Tgref.R

library(microclimfPara)

ndays <- 10
tsteps <- ndays * 24
rows <- 2; cols <- 2

# Pixel (1,1): fully valid, no gaps at all.
# Pixel (1,2): has an isolated single-hour NA in Tg on a snowdays hour (the
#   real gap shape confirmed present in production tile-101 data).
# Pixel (2,1): has an NA at hour 0 specifically (array position 1) - this
#   must NOT be treated as a fully-masked pixel; it must still be filled.
# Pixel (2,2): fully masked (all-NA throughout) - must stay all-NA.
Tg <- array(rep(seq_len(tsteps), each = rows * cols), dim = c(rows, cols, tsteps))
storage.mode(Tg) <- "double"
Tg[2, 2, ] <- NA_real_  # fully masked pixel

# Days: 1-4 and 8-10 are snowdays; 5-7 are pure nosnow days (the "gap").
snowdays <- c(1, 2, 3, 4, 8, 9, 10)
nosnowdays <- c(5, 6, 7)

# moutn$Tz is indexed compactly over nosnowdays only (3 days x 24h = 72
# hours), with a distinct, easily-recognized value (1000 + day) so the
# splice's calendar-position mapping can be checked by eye.
moutn_Tz <- array(0, dim = c(rows, cols, length(nosnowdays) * 24))
for (p in seq_along(nosnowdays)) {
  hrs <- ((p - 1) * 24 + 1):(p * 24)
  moutn_Tz[, , hrs] <- 1000 + nosnowdays[p]
}
# Pixel (2,2) is also masked in the no-snow model (fully NA)
moutn_Tz[2, 2, ] <- NA_real_

# Isolated single-hour NA in Tg on a snowdays hour: day 2 (a snowdays day),
# hour 5 (0-based), pixel (1,2).
day2_hour5_idx <- (2 - 1) * 24 + 5 + 1  # 1-based array position
Tg[1, 2, day2_hour5_idx] <- NA_real_

# NA at hour 0 specifically (array position 1), pixel (2,1) - otherwise valid.
Tg[2, 1, 1] <- NA_real_

Tgref <- microclimfPara:::.build_snow_Tgref(Tg, moutn_Tz, snowdays, nosnowdays)

# --- (1) Shape unchanged (still full-year, not subsetted).
stopifnot(identical(dim(Tgref), dim(Tg)))
cat("PASS: Tgref keeps the full-year shape (not subsetted to snowdays).\n")

# --- (2) Pure-nosnow days (5,6,7) are replaced by moutn$Tz's value at the
# right calendar position, for a pixel with no other gaps (pixel 1,1).
for (d in nosnowdays) {
  hrs <- ((d - 1) * 24 + 1):(d * 24)
  stopifnot(all(Tgref[1, 1, hrs] == 1000 + d))
}
cat("PASS: pure-nosnow days are spliced from moutn$Tz at the correct calendar position.\n")

# --- (3) Snowdays hours (for a gap-free pixel) keep the original Tg value.
snow_hrs_all <- as.vector(sapply(snowdays, function(d) ((d - 1) * 24 + 1):(d * 24)))
stopifnot(all(Tgref[1, 1, snow_hrs_all] == Tg[1, 1, snow_hrs_all]))
cat("PASS: snowdays/mixed-day hours keep Tg's own value (not overwritten by moutn$Tz).\n")

# --- (4) No NA remains anywhere for the non-masked pixels.
stopifnot(!anyNA(Tgref[1, 1, ]))
stopifnot(!anyNA(Tgref[1, 2, ]))
stopifnot(!anyNA(Tgref[2, 1, ]))
cat("PASS: zero NA remains for valid pixels (isolated Tg gap and hour-0 gap both filled).\n")

# --- (5) The isolated single-hour Tg NA (pixel 1,2, day 2 hour 5) was filled
# by linear interpolation, not left NA and not overwritten by moutn$Tz
# (day 2 is a snowdays day, so moutn$Tz must never appear there).
filled_val <- Tgref[1, 2, day2_hour5_idx]
stopifnot(!is.na(filled_val))
neighbor_before <- Tgref[1, 2, day2_hour5_idx - 1]
neighbor_after  <- Tgref[1, 2, day2_hour5_idx + 1]
stopifnot(filled_val > min(neighbor_before, neighbor_after) - 1e-9)
stopifnot(filled_val < max(neighbor_before, neighbor_after) + 1e-9)
cat(sprintf("PASS: isolated single-hour Tg gap on a snowdays hour filled by interpolation (%.4f, between its neighbours %.4f and %.4f).\n",
            filled_val, neighbor_before, neighbor_after))

# --- (6) The hour-0 NA (pixel 2,1) was filled, NOT treated as a fully-masked
# pixel - this is the critical case: a valid pixel whose first hour happens
# to be NA must still get every other NA hour filled too, since snowdaymov()
# (C++) uses hour 0 as its own per-pixel validity probe.
stopifnot(!is.na(Tgref[2, 1, 1]))
cat(sprintf("PASS: an NA at hour 0 on an otherwise-valid pixel is filled (%.4f), not treated as a masked pixel.\n", Tgref[2, 1, 1]))

# --- (7) The genuinely fully-masked pixel (2,2) stays entirely NA.
stopifnot(all(is.na(Tgref[2, 2, ])))
cat("PASS: a genuinely fully-masked pixel (all-NA throughout) stays all-NA.\n")

# --- (8) No-gap edge case: when nosnowdays is empty (pure snow season, no
# calendar gaps to splice), Tgref must skip the splice entirely but still
# run the NA-fill pass.
snowdays_all <- 1:ndays
Tgref_nogap <- microclimfPara:::.build_snow_Tgref(Tg, array(0, dim = c(rows, cols, 0)), snowdays_all, integer(0))
stopifnot(!anyNA(Tgref_nogap[1, 2, ]))    # isolated gap still filled
stopifnot(all(is.na(Tgref_nogap[2, 2, ])))  # masked pixel still all-NA
cat("PASS: empty-nosnowdays edge case (pure snow season, no gaps to splice) still runs the NA-fill pass correctly.\n")

cat("\nPASS: all .build_snow_Tgref() correctness checks succeeded.\n")
