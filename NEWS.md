# vecmatch (development version)

## Minor changes and bug fixes

- `balqual()` gained two new `type` values, `"desc_full"` and `"desc_reduced"`,
  which report standard descriptive statistics for every balancing variable,
  separately for each treatment level and for the unmatched and matched
  dataset. This makes it possible to compare the covariate distributions
  between both datasets, in addition to the pairwise balance metrics. The
  tables follow the layout of the balance tables, with the statistics as rows
  and the two matching stages as the `Before` and `After` columns. Numeric
  covariates are summarized by `N`, mean, standard deviation, minimum, first
  quartile, median, third quartile, maximum, skewness and excess kurtosis for
  `"desc_full"`, or by `N`, minimum, mean, median and maximum for
  `"desc_reduced"`. Categorical covariates are instead cross-tabulated, with
  the count and the percentage of every level within each treatment level,
  computed separately for each matching stage and printed as a single `N (%)`
  cell in a separate table. The two values are mutually exclusive, and both are
  opt-in, so the default output is unchanged.
- The descriptive metrics of `balqual()` describe the covariates as they are
  given in the `formula`, while the balance metrics are computed on the model
  matrix. A factor is therefore described once by its levels rather than once
  per dummy-coded column, and interaction terms are covered by the balance
  metrics only.
- The `statistic` argument of `balqual()` is not used by the descriptive
  metrics, and a warning is now issued when it is supplied together with
  `type = "desc_full"` or `type = "desc_reduced"` alone.
- Fixed a bug in `balqual()`, where the values of the `cutoffs` argument were
  matched to the metrics by position across all three balance metrics rather
  than by name. Passing a subset of the metrics, e.g.
  `type = c("smd", "var_ratio")`, raised a recycling warning and could
  evaluate a metric against the cutoff of another one.

# vecmatch 1.3.0

# vecmatch 1.3.0

## Major changes

- Refactored and unified the S3 class system across the package, and added
  helper methods for inspecting the internal structure of vecmatch objects.
- Added `get_select_params()` and `run_selected_matching()` to streamline
  the re-estimation step after the main optimization workflow.
- Reduced and cleaned up package dependencies in `DESCRIPTION`, and improved
  how suggested packages are handled in the code.

## Minor changes and bug fixes

- Fixed a bug in `raincloud()`, where facet labels were reversed when using
  `facet`.
- Removed backend handling from `optimize_gps()`. The parallel backend must
  now be registered outside the function.
- Updated the optimization vignette to use `run_selected_matching()`.
- Corrected a typo in the `cancer` dataset.
- Added examples for all exported functions and wrapped long-running examples
  in `\donttest{}`.
- Updated badges in the `README`.
- Added automated tests to check reproducibility of results.

# vecmatch 1.2.0

## Major changes

* Added `optimize_gps()`, `make_opt_args()`, and `select_opt()` to support a new
  GPS‐optimization workflow.
* Modified `csregion()` so the GPS can be reestimated after dropping 
  observations.

## Minor changes

* Fixed factor handling in `raincloud()` and `mosaic()`, now allowing custom
  facet ordering via releveling.
* Added SMD and p-value labels to `raincloud()`.
* Updated the `raincloud()` legend to show group names with their observation
  counts.


# vecmatch 1.1.0

## Major changes

* `csregion()` now allows specifying how to handle observations at the borders 
  of the Common Support Region (CSR) using the new `borders` argument.
* `match_gps()` has been updated to support datasets with only two unique 
  treatment groups.

## Minor changes

* Added a vignette demonstrating usage and functionality.
* Introduced this `NEWS.md` file to document package changes.
