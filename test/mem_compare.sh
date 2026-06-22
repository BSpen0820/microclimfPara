#!/usr/bin/env bash
# Compare peak memory before vs after the memory-management changes.
# "after"  = current HEAD (8da01b0)
# "before" = parent commit de2f791, checked out in a throwaway git worktree.
MAIN="D:/Code/PhD/R_Packages/microclimfPar"
BEFORE_DIR="D:/Code/PhD/R_Packages/microclimfPar_before"
BEFORE_COMMIT="de2f791"
RSCRIPT="/c/Program Files/R/R-4.5.2/bin/x64/Rscript.exe"
FIX="$MAIN/test/mem_fixture"
CSV="$MAIN/test/mem_compare_results.csv"
SCRIPT="$MAIN/test/mem_scenario.R"
SCENARIOS=(micro_tz_par micro_full_par runsnow_par micro_snow_par)

cleanup() { git -C "$MAIN" worktree remove --force "$BEFORE_DIR" 2>/dev/null || true; }
trap cleanup EXIT

rm -f "$CSV"

echo "### [1/5] Creating before-worktree at $BEFORE_COMMIT"
git -C "$MAIN" worktree remove --force "$BEFORE_DIR" 2>/dev/null || true
git -C "$MAIN" worktree add --force "$BEFORE_DIR" "$BEFORE_COMMIT"

echo "### [2/5] Building shared fixtures (HEAD)"
"$RSCRIPT" "$SCRIPT" "$MAIN" setup "$FIX" || { echo "fixture build failed"; exit 1; }

echo "### [3/5] Pre-compiling both package versions"
"$RSCRIPT" -e "suppressMessages(pkgload::load_all('$MAIN', quiet=TRUE))" >/dev/null 2>&1
"$RSCRIPT" -e "suppressMessages(pkgload::load_all('$BEFORE_DIR', quiet=TRUE))" >/dev/null 2>&1

echo "### [4/5] Running scenarios"
for s in "${SCENARIOS[@]}"; do "$RSCRIPT" "$SCRIPT" "$MAIN"       "$s" "$FIX" "$CSV" after  2; done
for s in "${SCENARIOS[@]}"; do "$RSCRIPT" "$SCRIPT" "$BEFORE_DIR" "$s" "$FIX" "$CSV" before 2; done

echo "### [5/5] Report"
"$RSCRIPT" "$MAIN/test/mem_report.R" "$CSV"
