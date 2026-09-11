# 01_quantify_differences.R
#
# PURPOSE: Quantify the magnitude of differences between base and parallel outputs.
# This distinguishes tiny FP rounding (< 1e-8) from large algorithmic errors.
# Run this BEFORE rebuilding the package, using the existing test/ RDS files.
#
# Run from: D:/Code/PhD/R_Packages/microclimfPar/
# setwd("D:/Code/PhD/R_Packages/microclimfPar/")

library(digest)

base_dir <- "./test/Base_Results"
par_dir  <- "./test/MicroPar_Par"

outputs <- c("mout_nosnow", "mout_snow", "smod_snow")

cat("=============================================================\n")
cat("Quantifying differences: Base vs Parallel\n")
cat("=============================================================\n\n")

results <- list()

for (nm in outputs) {
    base <- readRDS(file.path(base_dir, paste0(nm, ".rds")))
    par  <- readRDS(file.path(par_dir,  paste0(nm, ".rds")))

    cat("---", nm, "---\n")
    cat("  identical():", identical(base, par), "\n")

    # all.equal with zero tolerance reports the actual differences
    ae <- all.equal(base, par, tolerance = 0)
    if (isTRUE(ae)) {
        cat("  all.equal(tol=0): TRUE (bit-exact)\n")
    } else {
        cat("  all.equal(tol=0): differences found\n")
        # Print first 10 difference messages
        msgs <- if (length(ae) > 10) c(ae[1:10], "... (truncated)") else ae
        for (m in msgs) cat("    ", m, "\n")
    }

    # Drill into numeric elements to get max absolute and relative diff
    extract_nums <- function(x) {
        if (is.numeric(x)) return(list(x))
        if (is.list(x)) return(unlist(lapply(x, extract_nums), recursive = FALSE))
        return(list())
    }

    base_nums <- unlist(extract_nums(base))
    par_nums  <- unlist(extract_nums(par))

    if (length(base_nums) == length(par_nums) && length(base_nums) > 0) {
        diff <- abs(base_nums - par_nums)
        # Exclude NA pairs
        valid <- !is.na(base_nums) & !is.na(par_nums)
        diff_valid <- diff[valid]
        base_valid <- abs(base_nums[valid])

        max_abs  <- max(diff_valid, na.rm = TRUE)
        max_rel  <- max(diff_valid / pmax(base_valid, 1e-15), na.rm = TRUE)
        n_diff   <- sum(diff_valid > 0)
        n_total  <- length(diff_valid)

        cat("  Max absolute diff:", formatC(max_abs, format = "e", digits = 4), "\n")
        cat("  Max relative diff:", formatC(max_rel, format = "e", digits = 4), "\n")
        cat("  Cells differing:  ", n_diff, "/", n_total, "\n")

        # Classify
        if (max_abs == 0) {
            cat("  VERDICT: Bit-exact (identical)\n")
        } else if (max_abs < 1e-8) {
            cat("  VERDICT: Tiny FP rounding — consistent with compiler optimization differences\n")
            cat("           Fix: -ffp-contract=off in Makevars should resolve\n")
        } else if (max_abs < 1e-3) {
            cat("  VERDICT: Small but non-trivial differences — may indicate a logic bug\n")
            cat("           Investigate further with 02_ncores1_test.R\n")
        } else {
            cat("  VERDICT: LARGE differences — definite logic bug, not FP rounding\n")
            cat("           Run 02_ncores1_test.R immediately\n")
        }
    } else {
        cat("  Could not compare numerics (length mismatch or non-numeric)\n")
    }

    cat("\n")
    results[[nm]] <- list(identical = identical(base, par), ae = ae)
}

cat("=============================================================\n")
cat("Summary: if all diffs < 1e-8, rebuild with -ffp-contract=off\n")
cat("         and then run 03_revalidate_after_fix.R\n")
cat("=============================================================\n")
