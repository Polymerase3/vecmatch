# Optimize the Matching Process via Random Search

The `optimize_gps()` function performs a random search to identify
optimal combinations of parameters for the
[`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)
and
[`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
functions. The goal is to maximize the percentage of matched samples
(`perc_matched`) while minimizing the maximum standardized mean
difference (`smd`), thereby improving the overall balance of covariates
across treatment groups. The function supports parallel execution
through the `foreach` and `future` packages, enabling multithreaded
computation to accelerate the optimization process, particularly when
dealing with large datasets or complex parameter spaces.

## Usage

``` r
optimize_gps(
  data = NULL,
  formula,
  ordinal_treat = NULL,
  n_iter = 1000,
  opt_args = NULL
)
```

## Arguments

- data:

  A `data.frame` containing all variables specified in the `formula`
  argument. If `opt_args` is used, the `data` provided within `opt_args`
  must match this input exactly.

- formula:

  A valid formula object used to estimate the generalized propensity
  scores (GPS). The treatment variable appears on the left-hand side,
  and covariates on the right-hand side. Interactions can be specified
  using `*`. See
  [`stats::formula()`](https://rdrr.io/r/stats/formula.html) and
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
  for more details. If `opt_args` is provided, the formula within it
  must be identical to this argument.

- ordinal_treat:

  An atomic vector defining the ordered levels of the treatment
  variable. This confirms the variable is ordinal and adjusts its levels
  accordingly using
  `factor(treat, levels = ordinal_treat, ordered = TRUE)`. It is passed
  directly to
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md).
  If `NULL`, ordinal GPS estimation methods such as `polr` will be
  excluded from the optimization. See
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
  for details.

- n_iter:

  Integer. Number of unique parameter combinations to evaluate during
  optimization. Higher values generally yield better results but
  increase computation time. For large datasets or high-dimensional
  parameter spaces, increasing `n_iter` is recommended. When using
  parallel processing (`n_cores > 1`), performance gains become more
  apparent with larger `n_iter`. Too many cores with too few iterations
  may introduce overhead and reduce efficiency.

- opt_args:

  An object of class `"opt_args"` containing optimization parameters and
  argument settings. Use
  [`make_opt_args()`](https://polymerase3.github.io/vecmatch/reference/make_opt_args.md)
  to create this object. It specifies the search space for the GPS
  estimation and matching procedure.

## Value

An S3 object of class `best_opt_result`. The core component is a
`data.frame` containing the parameter combinations and results of the
optimization procedure. You can access it using
`attr(result, "opt_results")` or by calling `View(result)`, where
`result` is your `best_opt_result` object.

The object contains the following custom attributes:

- **`opt_results`**: A `data.frame` of optimization results. Each row
  corresponds to a unique parameter combination tested. For a complete
  description of columns, see the ***Details*** section.

- **`optimization_time`**: Time (in seconds) taken by the optimization
  loop (i.e., the core `for`-loop that evaluates combinations). This
  does **not** include the time needed for GPS estimation,
  pre-processing, or merging of results after loop completion. On large
  datasets, these excluded steps can still be substantial.

- **`combinations_tested`**: Total number of unique parameter
  combinations evaluated during optimization.

- **`smd_results`**: A detailed table of standardized mean differences
  (SMDs) for all pairwise treatment group comparisons and for all
  covariates specified in the `formula`. This is used by the
  [`select_opt()`](https://polymerase3.github.io/vecmatch/reference/select_opt.md)
  function to filter optimal models based on covariate-level balance
  across groups.

- **`treat_names`**: A character vector with the names of the unique
  treatment groups.

- **`model_covs`**: A character vector listing the model covariates
  (main effects and interactions) used in the `formula`. These names
  correspond to the variables shown in the `smd_results` table.

## Details

The output is an S3 object of class `best_opt_result`. Its core
component is a `data.frame` containing the parameter settings for the
best-performing models, grouped and ranked based on their balance
quality.

Optimization results are categorized into seven bins based on the
maximum standardized mean difference (SMD):

- 0.00-0.05

- 0.05-0.10

- 0.10-0.15

- 0.15-0.20

- 0.20-0.25

- 0.25-0.30

- Greater than 0.30

Within each SMD group, the parameter combination(s) achieving the
highest `perc_matched` (i.e., percentage of matched samples) is
selected. In cases where multiple combinations yield identical `smd` and
`perc_matched`, all such results are retained. Combinations where
matching failed or GPS estimation did not converge will return `NA` in
the result columns (e.g., `perc_matched`, `smd`).

The returned `data.frame` includes the following columns (depending on
the number of treatment levels):

- `iter_ID`: Unique identifier for each parameter combination

- `method_match`: Matching method used in
  [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md),
  e.g., `"nnm"` or `"fullopt"`

- `caliper`: Caliper value used in
  [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)

- `order`: Ordering of GPS scores prior to matching

- `kmeans_cluster`: Number of k-means clusters used

- `replace`: Whether replacement was used in matching (`nnm` only)

- `ties`: Tie-breaking rule in nearest-neighbor matching (`nnm` only)

- `ratio`: Control-to-treated ratio for `nnm`

- `min_controls`, `max_controls`: Minimum and maximum controls for
  `fullopt`

- `reference`: Reference group used in both
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
  and
  [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)

- `perc_matched`: Percentage of matched samples (from
  [`balqual()`](https://polymerase3.github.io/vecmatch/reference/balqual.md))

- `smd`: Maximum standardized mean difference (from
  [`balqual()`](https://polymerase3.github.io/vecmatch/reference/balqual.md))

- `p_{group_name}`: Percent matched per treatment group (based on group
  s ample size)

- `method_gps`: GPS estimation method used (from
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md))

- `link`: Link function used in GPS model

- `smd_group`: SMD range category for the row

The resulting `best_opt_result` object also includes a custom
[`summary()`](https://rdrr.io/r/base/summary.html) method that
summarizes:

- The number of optimal parameter sets per SMD group

- Their associated SMD and match rates

- Total combinations tested

- Total runtime of the optimization loop

## Examples

``` r
# Define formula for GPS estimation and matching
formula_cancer <- formula(status ~ age * sex)

# Set up the optimization parameter space
opt_args <- make_opt_args(cancer, formula_cancer, gps_method = "m1")

# Run optimization with 2000 random parameter sets and a fixed seed
# \donttest{
withr::with_seed(
  8252,
  {
    optimize_gps(
      data = cancer,
      formula = formula_cancer,
      opt_args = opt_args,
      n_iter = 2000
    )
  }
)
#> ℹ Initiating estimation of the GPS...
#> ℹ Initiating matching optimization...
#> ℹ Optimizing matching...
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: number of columns of result is not a multiple of vector length (arg 2)
#> Warning: number of columns of result is not a multiple of vector length (arg 2)
#> Warning: number of columns of result is not a multiple of vector length (arg 2)
#> Warning: number of columns of result is not a multiple of vector length (arg 2)
#> Warning: number of columns of result is not a multiple of vector length (arg 2)
#> Warning: number of columns of result is not a multiple of vector length (arg 2)
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: Matching failed. (Restrictions impossible to meet?)
#>  Enter ?matchfailed for more info.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
#> Warning: replace==FALSE, but there are more (weighted) treated obs than control obs.  Some treated obs will not be matched.  You may want to estimate ATC instead.
# }
```
