# Define the Optimization Parameter Space for Matching

`make_opt_args()` creates an object of class `"opt_args"` that defines
the parameter search space for
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).

The function accepts **vectors of values** for each customizable
argument involved in GPS estimation and matching. It computes the
**Cartesian product** of all parameter combinations, which serves as the
input search space for the random search algorithm used by
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).

To ensure valid optimization, the `data` and `formula` arguments must
exactly match those passed to
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).

## Usage

``` r
make_opt_args(
  data = NULL,
  formula,
  reference = NULL,
  gps_method = paste0("m", 1:10),
  matching_method = c("fullopt", "nnm"),
  caliper = seq(0.01, 10, 0.01),
  order = c("desc", "asc", "original", "random"),
  cluster = 2,
  replace = c(TRUE, FALSE),
  ties = c(TRUE, FALSE),
  ratio = 1,
  min_controls = 1,
  max_controls = 1
)
```

## Arguments

- data:

  A `data.frame` containing all variables referenced in `formula`. Must
  match the dataset used in
  [`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).

- formula:

  A valid formula specifying the treatment variable (left-hand side) and
  covariates (right-hand side). Interaction terms can be included using
  `*`. Must match the formula used in
  [`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).

- reference:

  A single string or vector of treatment group levels to be used as the
  reference (baseline) group in both GPS estimation and matching.

- gps_method:

  A string or vector of strings specifying GPS estimation methods.
  Allowed values are `"m1"` to `"m10"`. See *Details* below.

- matching_method:

  A string or vector of strings specifying the matching method(s) to
  evaluate. Currently supported options are `"nnm"` and `"fullopt"`. See
  [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md).

- caliper:

  A numeric value or vector of values specifying caliper widths (i.e.,
  maximum allowed GPS distance for matching). Same as in
  [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md),
  but allows multiple values.

- order:

  A string or vector of strings indicating the sorting order of
  logit-transformed GPS values before matching. Options are:

  - `"desc"`: sort from highest to lowest (default),

  - `"asc"`: sort from lowest to highest,

  - `"original"`: keep original order,

  - `"random"`: randomize order (use
    [`set.seed()`](https://rdrr.io/r/base/Random.html) for
    reproducibility).

- cluster:

  An integer or vector of integers specifying the number of clusters for
  k-means clustering (if applicable).

- replace:

  Logical value or vector of logicals indicating whether to allow
  matching with replacement. Same meaning as in
  [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md),
  but supports multiple settings.

- ties:

  Logical value or vector of logicals defining how ties should be
  handled during nearest-neighbor matching.

- ratio:

  A numeric value or vector specifying the ratio of control to treated
  units for matching (used in `"nnm"`).

- min_controls:

  A scalar or vector specifying the **minimum** number of controls to be
  matched to each treated unit (used in `"fullopt"`).

- max_controls:

  A scalar or vector specifying the **maximum** number of controls to be
  matched to each treated unit (used in `"fullopt"`).

## Value

An object of class `"opt_args"`, containing all valid parameter
combinations to be sampled by
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).
Use [`print()`](https://rdrr.io/r/base/print.html) to explore the
defined search space.

## Details

The returned object is of class `"opt_args"` and is intended to be
passed directly to
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md).
Internally, the function calculates the full Cartesian product of all
supplied parameter values and validates the structure of each.

The `gps_method` argument must contain one or more of the following
codes:

    | gps_method |      Method      |       Link Function         |
    |------------|------------------|-----------------------------|
    |    "m1"    |    multinom      |   generalized_logit         |
    |    "m2"    |     polr         |   logistic                  |
    |    "m3"    |     polr         |   probit                    |
    |    "m4"    |     polr         |   loglog                    |
    |    "m5"    |     polr         |   cloglog                   |
    |    "m6"    |     polr         |   cauchit                   |
    |    "m7"    |     vglm         |   multinomial_logit         |
    |    "m8"    |     vglm         |   reduced_rank_ml           |
    |    "m9"    |    brglm2        |   baseline_category_logit   |
    |   "m10"    |    mblogit       |   baseline_category_logit   |

The object includes a custom S3
[`print()`](https://rdrr.io/r/base/print.html) method that displays:

- A summary table of all allowed values for each optimization parameter,

- The total number of unique parameter combinations (i.e., the size of
  the search space).

## See also

[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md),
[`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md),
[`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)

## Examples

``` r

# Create search space with multiple values for GPS and matching
make_opt_args(
  data = cancer,
  formula = formula(status ~ age * sex),
  gps_method = c("m1", "m2", "m9"),
  matching_method = c("nnm", "fullopt"),
  caliper = c(0.1, 0.2),
  order = c("desc", "random"),
  reference = "control"
)
#> Optimization Argument Set (class: opt_args)
#> ---------------------------------------- 
#> gps_method      : m1, m2, m9
#> reference       : control
#> matching_method : nnm, fullopt
#> caliper         : 0.1, 0.2
#> order           : desc, random
#> cluster         : 2
#> replace         : TRUE, FALSE
#> ties            : TRUE, FALSE
#> ratio           : 1
#> min_controls    : 1
#> max_controls    : 1
#> ---------------------------------------- 
#> Total combinations: 60 
```
