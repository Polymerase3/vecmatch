# Changelog

## vecmatch (development version)

### Minor changes and bug fixes

- [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)
  now works with any `reference` when `method = "nnm"`. The columns kept
  for matching were selected with a logical mask built on the full
  `csmatrix`, but the mask was applied after the `treatment` column had
  already been dropped, so every position shifted by one. The gps column
  of the reference was therefore only retained when the reference
  happened to be the first gps column of the `csmatrix`; any other
  reference failed with `undefined columns selected`. When it did not
  fail, the shifted mask kept the neighbouring gps column as well, so
  \[Matching::Matchby()\] silently matched on two covariates instead of
  one and returned a different set of matches. Both branches now index
  the gps column of the first group of the current combination by name,
  which also fixes the same error for a user-supplied `combos` whose
  `group1` differs from `reference`.

- [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)
  no longer warns when a method-specific tuning parameter is left at its
  default. Previously, a call such as
  `match_gps(csr, method = "fullopt")` emitted a warning for every unset
  parameter (`order`, `min_controls`, `max_controls`), which made an
  ordinary default call look like it had gone wrong. Defaults are now
  filled in silently. Arguments that do not apply to the selected
  `method` are still reported, but only when they were explicitly
  supplied by the user.

- [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)
  now actually honours the `max_controls` argument for
  `method = "fullopt"`. It was passed to \[optmatch::fullmatch()\] under
  a name that function does not accept and carried the value of
  `caliper`, so it was silently discarded and had no effect on the
  matching. Calls that set `max_controls` to a non-default value will
  now produce different (correct) results.

- [`balqual()`](https://polymerase3.github.io/vecmatch/reference/balqual.md)
  gained two new `type` values, `"desc_full"` and `"desc_reduced"`,
  which report standard descriptive statistics for every balancing
  variable, separately for each treatment level and for the unmatched
  and matched dataset. This makes it possible to compare the covariate
  distributions between both datasets, in addition to the pairwise
  balance metrics. The tables follow the layout of the balance tables,
  with the statistics as rows and the two matching stages as the
  `Before` and `After` columns. Numeric covariates are summarized by
  `N`, mean, standard deviation, minimum, first quartile, median, third
  quartile, maximum, skewness and excess kurtosis for `"desc_full"`, or
  by `N`, minimum, mean, median and maximum for `"desc_reduced"`.
  Categorical covariates are instead cross-tabulated, with the count and
  the percentage of every level within each treatment level, computed
  separately for each matching stage and printed as a single `N (%)`
  cell in a separate table. The two values are mutually exclusive, and
  both are opt-in, so the default output is unchanged.

- The descriptive metrics of
  [`balqual()`](https://polymerase3.github.io/vecmatch/reference/balqual.md)
  describe the covariates as they are given in the `formula`, while the
  balance metrics are computed on the model matrix. A factor is
  therefore described once by its levels rather than once per
  dummy-coded column, and interaction terms are covered by the balance
  metrics only.

- The `statistic` argument of
  [`balqual()`](https://polymerase3.github.io/vecmatch/reference/balqual.md)
  is not used by the descriptive metrics, and a warning is now issued
  when it is supplied together with `type = "desc_full"` or
  `type = "desc_reduced"` alone.

- Fixed a bug in
  [`balqual()`](https://polymerase3.github.io/vecmatch/reference/balqual.md),
  where the values of the `cutoffs` argument were matched to the metrics
  by position across all three balance metrics rather than by name.
  Passing a subset of the metrics, e.g. `type = c("smd", "var_ratio")`,
  raised a recycling warning and could evaluate a metric against the
  cutoff of another one.

- `print.quality()` now returns its argument invisibly, as documented
  for [`print()`](https://rdrr.io/r/base/print.html) methods, instead of
  returning it visibly. The visible return value made the quality tables
  be printed twice whenever the result of
  [`print()`](https://rdrr.io/r/base/print.html) was itself
  auto-printed, e.g. inside
  [`capture.output()`](https://rdrr.io/r/utils/capture.output.html) or a
  knitted document.

- [`csregion()`](https://polymerase3.github.io/vecmatch/reference/csregion.md)
  now stops with an informative message when the common support region
  turns out to be empty, or when a treatment group loses all of its
  observations. Previously such a call failed with
  `non-numeric argument to binary operator`, which came from the
  fallback warning of the internal refitting step rather than from the
  empty region itself, and gave no hint about the actual cause. With
  `refit = FALSE` the same situation returned a zero-row `csr` object
  silently, and the problem only surfaced later in
  [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md).
  Emptying a group is now reported as well, since dropping it would
  silently change the estimand. The new messages name the groups
  involved and suggest a larger sample, fewer treatment groups, a
  different `borders` value or a different model specification in
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md).

- Fixed the fallback warning of
  [`csregion()`](https://polymerase3.github.io/vecmatch/reference/csregion.md),
  which was emitted when refitting the GPS model on the CSR-restricted
  data fails. Its message was split into three fragments that
  [`strwrap()`](https://rdrr.io/r/base/strwrap.html) received as the
  `width` and `indent` arguments instead of as text, so the warning
  aborted instead of being shown.

- [`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md)
  no longer aborts with `$ operator is invalid for atomic vectors` when
  the GPS estimation fails for the first combination of the search
  space. Failed combinations are passed through as the message of the
  underlying error, but the treatment levels were read from the first
  result unconditionally, so a single failure at the head of the
  estimation space brought down the whole sweep. The first result that
  actually holds a GPS matrix is now used instead. If every combination
  fails, the function stops with a message that reports how many were
  tried and the reason given by the first failure.

- The package documentation is now also available as a `pkgdown` website
  at <https://polymerase3.github.io/vecmatch/>. The site collects the
  reference pages of all exported functions, grouped by the five steps
  of the `vecmatch` workflow, together with the vignettes and this
  changelog, and is rebuilt automatically by a GitHub Actions workflow
  on every push to `main`.

## vecmatch 1.3.0

CRAN release: 2025-12-01

## vecmatch 1.3.0

CRAN release: 2025-12-01

### Major changes

- Refactored and unified the S3 class system across the package, and
  added helper methods for inspecting the internal structure of vecmatch
  objects.
- Added
  [`get_select_params()`](https://polymerase3.github.io/vecmatch/reference/get_select_params.md)
  and
  [`run_selected_matching()`](https://polymerase3.github.io/vecmatch/reference/run_selected_matching.md)
  to streamline the re-estimation step after the main optimization
  workflow.
- Reduced and cleaned up package dependencies in `DESCRIPTION`, and
  improved how suggested packages are handled in the code.

### Minor changes and bug fixes

- Fixed a bug in
  [`raincloud()`](https://polymerase3.github.io/vecmatch/reference/raincloud.md),
  where facet labels were reversed when using `facet`.
- Removed backend handling from
  [`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).
  The parallel backend must now be registered outside the function.
- Updated the optimization vignette to use
  [`run_selected_matching()`](https://polymerase3.github.io/vecmatch/reference/run_selected_matching.md).
- Corrected a typo in the `cancer` dataset.
- Added examples for all exported functions and wrapped long-running
  examples in `\donttest{}`.
- Updated badges in the `README`.
- Added automated tests to check reproducibility of results.

## vecmatch 1.2.0

CRAN release: 2025-07-08

### Major changes

- Added
  [`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md),
  [`make_opt_args()`](https://polymerase3.github.io/vecmatch/reference/make_opt_args.md),
  and
  [`select_opt()`](https://polymerase3.github.io/vecmatch/reference/select_opt.md)
  to support a new GPS‐optimization workflow.
- Modified
  [`csregion()`](https://polymerase3.github.io/vecmatch/reference/csregion.md)
  so the GPS can be reestimated after dropping observations.

### Minor changes

- Fixed factor handling in
  [`raincloud()`](https://polymerase3.github.io/vecmatch/reference/raincloud.md)
  and
  [`mosaic()`](https://polymerase3.github.io/vecmatch/reference/mosaic.md),
  now allowing custom facet ordering via releveling.
- Added SMD and p-value labels to
  [`raincloud()`](https://polymerase3.github.io/vecmatch/reference/raincloud.md).
- Updated the
  [`raincloud()`](https://polymerase3.github.io/vecmatch/reference/raincloud.md)
  legend to show group names with their observation counts.

## vecmatch 1.1.0

CRAN release: 2025-04-24

### Major changes

- [`csregion()`](https://polymerase3.github.io/vecmatch/reference/csregion.md)
  now allows specifying how to handle observations at the borders of the
  Common Support Region (CSR) using the new `borders` argument.
- [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)
  has been updated to support datasets with only two unique treatment
  groups.

### Minor changes

- Added a vignette demonstrating usage and functionality.
- Introduced this `NEWS.md` file to document package changes.
