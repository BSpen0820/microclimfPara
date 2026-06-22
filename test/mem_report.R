# Summarise before/after peak-memory comparison from mem_scenario.R output.
args <- commandArgs(trailingOnly = TRUE)
csv  <- if (length(args) >= 1) args[1] else "test/mem_compare_results.csv"
d <- read.csv(csv, stringsAsFactors = FALSE)

fmt <- function(s) {
  b <- d[d$scenario == s & d$version == "before", ]
  a <- d[d$scenario == s & d$version == "after",  ]
  if (nrow(b) == 0 || nrow(a) == 0) return(NULL)
  data.frame(
    scenario        = s,
    before_RSS_Mb   = b$peak_rss_mb,
    after_RSS_Mb    = a$peak_rss_mb,
    saved_RSS_Mb    = round(b$peak_rss_mb - a$peak_rss_mb, 1),
    RSS_reduction   = sprintf("%.1f%%", 100 * (b$peak_rss_mb - a$peak_rss_mb) / b$peak_rss_mb),
    before_Rheap_Mb = b$peak_r_mb,
    after_Rheap_Mb  = a$peak_r_mb,
    saved_Rheap_Mb  = round(b$peak_r_mb - a$peak_r_mb, 1)
  )
}
tab <- do.call(rbind, lapply(unique(d$scenario), fmt))
cat("\n================= PEAK MEMORY: before vs after =================\n")
print(tab, row.names = FALSE)
cat("\nRSS  = process peak working set (ps peak_wset) - true OS peak\n")
cat("Rheap = R gc() max-used high-water (Ncells+Vcells)\n")
