extract_ne_data <- function(netable, ci = c("parametric", "jackknife")) {
  ci <- match.arg(ci)
  
  lo <- if (ci == "parametric") "CI low Parametric" else "CI low JackKnife"
  hi <- if (ci == "parametric") "CI high Parametric" else "CI high JackKnife"
  
  tb <- netable[[1]]
  # A couple of NeEstimator tables carry a ragged trailing column with no
  # name (a duplicate of the last frequency column, truncated). tidyselect
  # cannot index a data frame whose names contain NA, so drop it.
  tb <- tb[, !is.na(names(tb)), drop = FALSE]
  
  tb |>
    # One row per frequency threshold, one column per reported statistic
    pivot_longer(-Statistic, names_to = "freq_col", values_to = "value") |>
    pivot_wider(names_from = Statistic, values_from = value) |>
    transmute(
      allele_freq = `Lowest Allele Frequency Used`,
      ne = as.numeric(`Estimated Ne^`),
      ne_min = as.numeric(.data[[lo]]),
      ne_max = as.numeric(.data[[hi]])
    )
}