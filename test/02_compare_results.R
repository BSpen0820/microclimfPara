# 02_compare_results.R
# Run after 01_original_vs_fork.R. Two parts:
#
#   Part 1 — hashes every .rds under Base_Results/, MicroPar_Ser/,
#   MicroPar_Par/ and compares each whole object against the matching
#   Base_Results/ file: exact MD5 match, and (if that fails) an
#   all.equal(tolerance = 1e-10) near-equality check. Writes
#   Validation_Results.csv.
#
#   Part 2 — for mout_nosnow, mout_snow and smod_snow (the grid/snow
#   outputs), breaks the comparison down by individual output variable
#   (Tz, tleaf, Tc, Tg, ...) so it's possible to see which specific
#   variables differ from the original microclimf, and by how much,
#   rather than only a single pass/fail for the whole object. Writes
#   Validation_Results_byvariable.csv. See README.md for why some
#   variables are expected to differ (the snow ground-temperature
#   smoothing fix) and by roughly how much.
#
# Run from the test/ directory, after 01_original_vs_fork.R has populated
# Base_Results/, MicroPar_Ser/, MicroPar_Par/ there.
#
# Run: cd test; Rscript 02_compare_results.R

library(tools)
library(digest)

# =============================================================================
# Part 1: whole-object hash / near-equality comparison
# =============================================================================

files <- list.files(path = "./", pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

files.df <- data.frame(
  file = files,
  name = file_path_sans_ext(basename(files)),
  dir = dirname(files)
)

files.df$hash <- vapply(files.df$file, digest, character(1), algo = "md5", file = TRUE)

files.df$identical_2_Base <- NA
files.df$identical_Toler_Base <- NA

for (unique_name in unique(files.df$name)) {

  files_subset <- files.df[files.df$name == unique_name, ]
  base_row <- files_subset[files_subset$dir == "./Base_Results", ]
  if (nrow(base_row) == 0) next
  base_hash <- base_row$hash[1]

  for (j in seq_len(nrow(files_subset))) {
    if (files_subset$dir[j] == "./Base_Results") next
    compare_hash <- files_subset$hash[j]
    same <- identical(base_hash, compare_hash)
    idx <- files.df$file == files_subset$file[j]
    if (same) {
      files.df$identical_2_Base[idx] <- TRUE
      files.df$identical_Toler_Base[idx] <- TRUE
    } else {
      base_data <- readRDS(base_row$file[1])
      compare_data <- readRDS(files_subset$file[j])
      # learned that cpp with parallel worker will return values a little different due to passing values with "less" precision,
      # but essentially the same results, so I set a tolerance level to check for "near" equality.
      cmp <- all.equal(base_data, compare_data, tolerance = 1e-10)
      files.df$identical_2_Base[idx] <- FALSE
      files.df$identical_Toler_Base[idx] <- if (isTRUE(cmp)) TRUE else paste(cmp, collapse = "; ")
      rm(base_data, compare_data)
      invisible(gc(full = TRUE))
    }
  }
}

write.csv(files.df, "Validation_Results.csv", row.names = FALSE)
cat("\n=== Part 1: whole-object comparison (Validation_Results.csv) ===\n")
print(files.df[, c("name", "dir", "identical_2_Base")], row.names = FALSE)

# =============================================================================
# Part 2: per-variable breakdown for mout_nosnow, mout_snow, smod_snow
# =============================================================================

# Compares one numeric array/vector between the original and the fork.
# NA counts are reported separately from the numeric diff stats because a
# changed NA pattern (e.g. the fork filling in a previously-NA cell) isn't
# comparable via a simple difference. Relative difference excludes
# near-zero base values (unstable/misleading when the base value is close
# to 0, e.g. ground/snowpack temperature straddling 0 degC, or snow depth/SWE
# that's exactly or nearly 0 whenever there's no/little snow). A fixed
# absolute threshold isn't enough for that on its own (e.g. 1e-8 vs 1e-6 SWE
# both "pass" a 1e-6 cutoff and can still blow a relative difference up by
# orders of magnitude for that single cell), so the threshold also scales
# with the variable's own typical magnitude -- and the summary reported is
# the MEDIAN relative difference, not the mean: a mean is dominated by a
# handful of near-zero-denominator outliers exactly like this, while the
# median reflects the typical case.
compare_variable <- function(base_val, fork_val) {
  if (!is.numeric(base_val) || !is.numeric(fork_val)) return(NULL)
  base_na <- is.na(base_val)
  fork_na <- is.na(fork_val)
  both_valid <- !base_na & !fork_na
  n_compared <- sum(both_valid)

  if (n_compared > 0) {
    absdiff <- abs(base_val[both_valid] - fork_val[both_valid])
    scale <- stats::median(abs(base_val[both_valid]))
    if (!is.finite(scale) || scale == 0) scale <- max(abs(base_val[both_valid]), 0)
    thresh <- max(1e-6, 0.01 * scale)
    nonzero <- both_valid & abs(base_val) > thresh
    reldiff <- if (sum(nonzero) >= max(30, 0.01 * n_compared)) {
      abs(base_val[nonzero] - fork_val[nonzero]) / abs(base_val[nonzero])
    } else {
      numeric(0)
    }
  } else {
    absdiff <- numeric(0)
    reldiff <- numeric(0)
  }

  data.frame(
    n_total          = length(base_val),
    n_na_base        = sum(base_na),
    n_na_fork        = sum(fork_na),
    na_mismatch      = sum(base_na != fork_na),
    n_compared       = n_compared,
    exact_match      = sum(base_na != fork_na) == 0 && (n_compared == 0 || max(absdiff) == 0),
    mean_abs_diff      = if (n_compared > 0) mean(absdiff) else NA,
    median_abs_diff    = if (n_compared > 0) stats::median(absdiff) else NA,
    max_abs_diff       = if (n_compared > 0) max(absdiff) else NA,
    median_rel_diff_pct = if (length(reldiff) > 0) 100 * stats::median(reldiff) else NA
  )
}

group_names <- c("mout_nosnow", "mout_snow", "smod_snow")
compare_dirs <- c("MicroPar_Ser", "MicroPar_Par")
byvar_rows <- list()

for (nm in group_names) {
  base_path <- file.path("Base_Results", paste0(nm, ".rds"))
  if (!file.exists(base_path)) next
  base_obj <- readRDS(base_path)

  for (d in compare_dirs) {
    fork_path <- file.path(d, paste0(nm, ".rds"))
    if (!file.exists(fork_path)) next
    fork_obj <- readRDS(fork_path)

    for (var in names(base_obj)) {
      res <- compare_variable(base_obj[[var]], fork_obj[[var]])
      if (is.null(res)) next  # skip non-numeric fields (e.g. tme)
      res$name <- nm
      res$dir <- d
      res$variable <- var
      byvar_rows[[length(byvar_rows) + 1]] <- res
    }

    rm(fork_obj)
    invisible(gc(full = TRUE))
  }

  rm(base_obj)
  invisible(gc(full = TRUE))
}

byvar.df <- do.call(rbind, byvar_rows)
byvar.df <- byvar.df[, c("name", "dir", "variable", "n_total", "n_na_base", "n_na_fork",
                          "na_mismatch", "n_compared", "exact_match",
                          "mean_abs_diff", "median_abs_diff", "max_abs_diff", "median_rel_diff_pct")]

write.csv(byvar.df, "Validation_Results_byvariable.csv", row.names = FALSE)
cat("\n=== Part 2: per-variable comparison, fork vs. original microclimf (Validation_Results_byvariable.csv) ===\n")
print(byvar.df, row.names = FALSE)
