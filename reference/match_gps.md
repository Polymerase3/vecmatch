# Match the Data Based on Generalized Propensity Scores

The `match_gps()` function performs sample matching based on generalized
propensity scores (GPS). It applies a k-means clustering step to the GPS
to partition the data into clusters and then matches all treatment
groups within each cluster. This yields balanced comparisons across
treatment levels while preserving the GPS structure. To our knowledge,
`vecmatch` provides the first implementation in R of the vector-matching
algorithm described by Lopez and Gutman (2017).

## Usage

``` r
match_gps(
  csmatrix = NULL,
  method = "nnm",
  caliper = 0.2,
  reference = NULL,
  ratio = NULL,
  replace = NULL,
  order = NULL,
  ties = NULL,
  min_controls = NULL,
  max_controls = NULL,
  kmeans_args = NULL,
  kmeans_cluster = 5,
  verbose_output = FALSE,
  ...
)
```

## Arguments

- csmatrix:

  An object of class `gps` and/or `csr` representing a data frame of
  generalized propensity scores. The first column must be the treatment
  variable, with additional attributes describing the calculation of the
  common support region and the estimation of generalized propensity
  scores. It is crucial that the common support region was calculated
  using the
  [`csregion()`](https://polymerase3.github.io/vecmatch/reference/csregion.md)
  function to ensure compatibility.

- method:

  A single string specifying the matching method to use. The default is
  `"nnm"`, which applies the k-nearest neighbors matching algorithm. See
  the Details section for a full list of available methods.

- caliper:

  A numeric value specifying the caliper width, which defines the
  allowable range within which observations can be matched. It is
  expressed as a percentage of the standard deviation of the
  logit-transformed generalized propensity scores. To perform matching
  without a caliper, set this parameter to a very large value. For exact
  matching, set `caliper = 0` and enable the `exact` option by setting
  it to `TRUE`.

- reference:

  A single string specifying the exact level of the treatment variable
  to be used as the reference in the matching process. All other
  treatment levels will be matched to this reference level. Ideally,
  this should be the control level. If no natural control is present,
  avoid selecting a level with extremely low or high covariate or
  propensity score values. Instead, choose a level with covariate or
  propensity score distributions that are centrally positioned among all
  treatment groups to maximize the number of matches.

- ratio:

  A scalar for the number of matches which should be found for each
  control observation. The default is one-to-one matching. Only
  available for the methods `"nnm"` and `"pairopt"`.

- replace:

  Logical value indicating whether matching should be done with
  replacement. If `FALSE`, the order of matches generally matters.
  Matches are found in the same order as the data is sorted.
  Specifically, the matches for the first observation will be found
  first, followed by those for the second observation, and so on.
  Matching without replacement is generally not recommended as it tends
  to increase bias. However, in cases where the dataset is large and
  there are many potential matches, setting `replace = FALSE` often
  results in a substantial speedup with negligible or no bias. Only
  available for the method `"nnm"`

- order:

  A string specifying the order in which logit-transformed GPS values
  are sorted before matching. The available options are:

  - `"desc"` – sorts GPS values from highest to lowest (default).

  - `"asc"` – sorts GPS values from lowest to highest.

  - `"original"` – preserves the original order of GPS values.

  - `"random"` – randomly shuffles GPS values. To generate different
    random orders, set a seed using
    [`set.seed()`](https://rdrr.io/r/base/Random.html).

- ties:

  A logical flag indicating how tied matches should be handled.
  Available only for the `"nnm"` method, with a default value of `FALSE`
  (all tied matches are included in the final dataset, but only unique
  observations are retained). For more details, see the `ties` argument
  in
  [`Matching::Matchby()`](https://rdrr.io/pkg/Matching/man/Matchby.html).

- min_controls:

  The minimum number of treatment observations that should be matched to
  each control observation. Available only for the `"fullopt"` method.
  For more details, see the `min.controls` argument in
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).

- max_controls:

  The maximum number of treatment observations that can be matched to
  each control observation. Available only for the `"fullopt"` method.
  For more details, see the `max.controls` argument in
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).

- kmeans_args:

  A list of arguments to pass to
  [stats::kmeans](https://rdrr.io/r/stats/kmeans.html). These arguments
  must be provided inside a [`list()`](https://rdrr.io/r/base/list.html)
  in the paired `name = value` format.

- kmeans_cluster:

  An integer specifying the number of clusters to pass to
  [stats::kmeans](https://rdrr.io/r/stats/kmeans.html).

- verbose_output:

  a logical flag. If `TRUE` a more verbose version of the function is
  run and the output is printed out to the console.

- ...:

  Additional arguments to be passed to the matching function.

## Value

A `data.frame` similar to the one provided as the `data` argument in the
[`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
function, containing the same columns but only the observations for
which a match was found. The returned object includes two attributes,
accessible with the [`attr()`](https://rdrr.io/r/base/attr.html)
function:

- `original_data`: A `data.frame` with the original data returned by the
  [`csregion()`](https://polymerase3.github.io/vecmatch/reference/csregion.md)
  or
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
  function, after the estimation of the csr and filtering out
  observations not within the csr.

- `matching_filter`: A logical vector indicating which rows from
  `original_data` were included in the final matched dataset.

## Details

Propensity score matching can be performed using various matching
algorithms. Lopez and Gutman (2017) do not explicitly specify the
matching algorithm used, but it is assumed they applied the commonly
used k-nearest neighbors matching algorithm, implemented as
`method = "nnm"`. However, this algorithm can sometimes be challenging
to use, especially when treatment and control groups have unequal sizes.
When `replace = FALSE`, the number of matches is strictly limited by the
smaller group, and even with `replace = TRUE`, the results may not
always be satisfactory. To address these limitations, we have
implemented an additional matching algorithm to maximize the number of
matched observations within a dataset.

The available matching methods are:

- `"nnm"` – classic k-nearest neighbors matching, implemented using
  [`Matching::Matchby()`](https://rdrr.io/pkg/Matching/man/Matchby.html).
  The tunable parameters in `match_gps()` are `caliper`, `ratio`,
  `replace`, `order`, and `ties`. Additional arguments can be passed to
  [`Matching::Matchby()`](https://rdrr.io/pkg/Matching/man/Matchby.html)
  via the `...` argument.

- `"fullopt"` – optimal full matching algorithm, implemented with
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).
  This method calculates a discrepancy matrix to identify all possible
  matches, often optimizing the percentage of matched observations. The
  available tuning parameters are `caliper`, `min_controls`, and
  `max_controls`.

- `"pairmatch"` – optimal 1:1 and 1:k matching algorithm, implemented
  using
  [`optmatch::pairmatch()`](https://rdrr.io/pkg/optmatch/man/pairmatch.html),
  which is actually a wrapper around
  [`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html).
  Like `"fullopt"`, this method calculates a discrepancy matrix and
  finds matches that minimize its sum. The available tuning parameters
  are `caliper` and `ratio`.

Each method only accepts the tuning parameters listed above. Arguments
that do not apply to the selected `method` are dropped, and a warning
naming them is emitted only when they were explicitly supplied. Tuning
parameters that are left at their default are filled in silently, so a
call relying on the defaults runs without any warnings.

## References

Michael J. Lopez, Roee Gutman "Estimation of Causal Effects with
Multiple Treatments: A Review and New Ideas," Statistical Science,
Statist. Sci. 32(3), 432-454, (August 2017)

## See also

[`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
for the calculation of generalized propensity scores;
`MatchIt::matchit()`,
[`optmatch::fullmatch()`](https://rdrr.io/pkg/optmatch/man/fullmatch.html)
and
[`optmatch::pairmatch()`](https://rdrr.io/pkg/optmatch/man/pairmatch.html)
for the documentation of the matching functions;
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) for the
documentation of the k-Means algorithm.

## Examples

``` r

# Step 1.) Estimation of the generalized propensity scores
gp_scores <- estimate_gps(
  formula(status ~ age * sex),
  data = cancer,
  method = "multinom",
  reference = "control",
  verbose_output = TRUE
)
#> # weights:  20 (12 variable)
#> initial  value 1696.824298 
#> iter  10 value 1675.784597
#> final  value 1675.690057 
#> converged

# Step 2.) Defining the common support region
gps_csr <- csregion(gp_scores)
#> # weights:  20 (12 variable)
#> initial  value 1669.098411 
#> iter  10 value 1648.708826
#> final  value 1646.780796 
#> converged

# Example 1: optimal full matching ("fullopt") with explicit tuning params
matched_cancer_full <- match_gps(
  gps_csr,
  caliper = 0.25,
  reference = "control",
  method = "fullopt",
  kmeans_cluster = 2,
  kmeans_args = list(
    iter.max = 200,
    algorithm = "Forgy"
  ),
  min_controls = 0,
  max_controls = Inf,
  order = "desc",
  verbose_output = TRUE
)

# Example 2: nearest neighbour matching ("nnm") demonstrating ratio /
# replace / order / ties
matched_cancer_nnm <- match_gps(
  gps_csr,
  caliper = 0.25,
  reference = "control",
  method = "nnm",
  ratio = 2,
  replace = TRUE,
  order = "asc",
  ties = TRUE,
  verbose_output = FALSE
)
```
