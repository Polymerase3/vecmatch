# Filter the Data Based on Common Support Region

The `csregion()` function estimates the boundaries of the rectangular
common support region, as defined by Lopez and Gutman (2017), and
filters the matrix of generalized propensity scores based on these
boundaries. The function returns a matrix of observations whose
generalized propensity scores lie within the treatment group-specific
boundaries.

## Usage

``` r
csregion(gps_matrix, borders = "include", refit = TRUE)
```

## Arguments

- gps_matrix:

  An object of classes `gps` and `data.frame` (e.g., created by the
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
  function). The first column corresponds to the treatment or grouping
  variable, while the other columns represent the treatment assignment
  probabilities calculated separately for each hypotetical treatment
  group. The number of columns should therefore be equal to the number
  of unique levels of the treatment variable plus one (for the treatment
  variable itself). The number of rows should correspond to the number
  of subjects for which generalized propensity scores were estimated.

- borders:

  A character string specifying how to handle observations at the edges
  of the Common Support Region (CSR). Acceptable values are `"include"`
  and `"exclude"`. If `"include"` is selected (default), observations
  with Generalized Propensity Scores (GPS) exactly equal to the CSR
  boundaries are retained for further analysis. This corresponds to a
  non-strict inequality: `lower_bound <= GPS <= upper_bound`. If
  `"exclude"` is selected, observations lying exactly on the CSR
  boundaries are removed. This corresponds to a strict inequality:
  `lower_bound < GPS < upper_bound`. Using `"exclude"` will typically
  result in a slightly smaller matched sample size compared to
  `"include"`, but may be preferred for more conservative matching.

- refit:

  Logical. If `TRUE` (default), the model used to estimate the GPS is
  refitted after excluding samples outside the common support region,
  using the same formula and method as in the original
  [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
  call. If `FALSE`, the model is not refitted, but still only samples
  within the CSR are retained. Refitting is recommended, as suggested by
  Lopez and Gutman (2017).

## Value

A numeric matrix similar to the one returned by
[`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md),
but with the number of rows reduced to exclude those observations that
do not fit within the common support region (CSR) boundaries. The
returned object also possesses additional attributes that summarize the
calculation process of the CSR boundaries:

- `filter_matrix` - A logical matrix with the same dimensions as the
  gps-part of `gps_matrix`, indicating which treatment assignment
  probabilities fall within the CSR boundaries,

- `filter_vector` - A vector indicating whether each observation was
  kept (`TRUE`) or removed (`FALSE`), essentially a row-wise sum of
  `filter_matrix`,

- `csr_summary` - A summary of the CSR calculation process, including
  details of the boundaries and the number of observations filtered.

- `csr_data` - The original dataset used for the estimation of
  generalized propensity scores (`original_data` attribute of the `gps`
  object) filtered by the `filter_vector`

## Examples

``` r
# We could estimate simples generalized propensity scores for the `iris`
# dataset
gps <- estimate_gps(Species ~ Sepal.Length, data = iris)

# And then define the common support region boundaries using `csregion()`
gps_csr <- csregion(gps)
#> Warning: Some groups have fewer than 20 observations, which may impact the performance of the matching process. Consider using `replace = TRUE`in `match_gps()` to address this.

# The additional information of the CSR-calculation process are
# accessible through the attributes described in the `*Value*` section
attr(gps_csr, "filter_matrix")
#>        setosa versicolor virginica
#>   [1,]   TRUE       TRUE      TRUE
#>   [2,]   TRUE       TRUE      TRUE
#>   [3,]  FALSE      FALSE     FALSE
#>   [4,]  FALSE      FALSE     FALSE
#>   [5,]   TRUE       TRUE      TRUE
#>   [6,]   TRUE       TRUE      TRUE
#>   [7,]  FALSE      FALSE     FALSE
#>   [8,]   TRUE       TRUE      TRUE
#>   [9,]  FALSE      FALSE     FALSE
#>  [10,]   TRUE       TRUE      TRUE
#>  [11,]   TRUE       TRUE      TRUE
#>  [12,]  FALSE      FALSE     FALSE
#>  [13,]  FALSE      FALSE     FALSE
#>  [14,]  FALSE      FALSE     FALSE
#>  [15,]   TRUE       TRUE      TRUE
#>  [16,]   TRUE       TRUE      TRUE
#>  [17,]   TRUE       TRUE      TRUE
#>  [18,]   TRUE       TRUE      TRUE
#>  [19,]   TRUE       TRUE      TRUE
#>  [20,]   TRUE       TRUE      TRUE
#>  [21,]   TRUE       TRUE      TRUE
#>  [22,]   TRUE       TRUE      TRUE
#>  [23,]  FALSE      FALSE     FALSE
#>  [24,]   TRUE       TRUE      TRUE
#>  [25,]  FALSE      FALSE     FALSE
#>  [26,]   TRUE       TRUE      TRUE
#>  [27,]   TRUE       TRUE      TRUE
#>  [28,]   TRUE       TRUE      TRUE
#>  [29,]   TRUE       TRUE      TRUE
#>  [30,]  FALSE      FALSE     FALSE
#>  [31,]  FALSE      FALSE     FALSE
#>  [32,]   TRUE       TRUE      TRUE
#>  [33,]   TRUE       TRUE      TRUE
#>  [34,]   TRUE       TRUE      TRUE
#>  [35,]   TRUE       TRUE      TRUE
#>  [36,]   TRUE       TRUE      TRUE
#>  [37,]   TRUE       TRUE      TRUE
#>  [38,]   TRUE       TRUE      TRUE
#>  [39,]  FALSE      FALSE     FALSE
#>  [40,]   TRUE       TRUE      TRUE
#>  [41,]   TRUE       TRUE      TRUE
#>  [42,]  FALSE      FALSE     FALSE
#>  [43,]  FALSE      FALSE     FALSE
#>  [44,]   TRUE       TRUE      TRUE
#>  [45,]   TRUE       TRUE      TRUE
#>  [46,]  FALSE      FALSE     FALSE
#>  [47,]   TRUE       TRUE      TRUE
#>  [48,]  FALSE      FALSE     FALSE
#>  [49,]   TRUE       TRUE      TRUE
#>  [50,]   TRUE       TRUE      TRUE
#>  [51,]  FALSE       TRUE     FALSE
#>  [52,]  FALSE       TRUE     FALSE
#>  [53,]  FALSE       TRUE     FALSE
#>  [54,]   TRUE       TRUE      TRUE
#>  [55,]  FALSE       TRUE     FALSE
#>  [56,]   TRUE       TRUE      TRUE
#>  [57,]  FALSE       TRUE     FALSE
#>  [58,]   TRUE       TRUE      TRUE
#>  [59,]  FALSE       TRUE     FALSE
#>  [60,]   TRUE       TRUE      TRUE
#>  [61,]   TRUE       TRUE      TRUE
#>  [62,]  FALSE       TRUE     FALSE
#>  [63,]  FALSE       TRUE     FALSE
#>  [64,]  FALSE       TRUE     FALSE
#>  [65,]   TRUE       TRUE      TRUE
#>  [66,]  FALSE       TRUE     FALSE
#>  [67,]   TRUE       TRUE      TRUE
#>  [68,]   TRUE       TRUE      TRUE
#>  [69,]  FALSE       TRUE     FALSE
#>  [70,]   TRUE       TRUE      TRUE
#>  [71,]  FALSE       TRUE     FALSE
#>  [72,]  FALSE       TRUE     FALSE
#>  [73,]  FALSE       TRUE     FALSE
#>  [74,]  FALSE       TRUE     FALSE
#>  [75,]  FALSE       TRUE     FALSE
#>  [76,]  FALSE       TRUE     FALSE
#>  [77,]  FALSE       TRUE     FALSE
#>  [78,]  FALSE       TRUE     FALSE
#>  [79,]  FALSE       TRUE     FALSE
#>  [80,]   TRUE       TRUE      TRUE
#>  [81,]   TRUE       TRUE      TRUE
#>  [82,]   TRUE       TRUE      TRUE
#>  [83,]   TRUE       TRUE      TRUE
#>  [84,]  FALSE       TRUE     FALSE
#>  [85,]   TRUE       TRUE      TRUE
#>  [86,]  FALSE       TRUE     FALSE
#>  [87,]  FALSE       TRUE     FALSE
#>  [88,]  FALSE       TRUE     FALSE
#>  [89,]   TRUE       TRUE      TRUE
#>  [90,]   TRUE       TRUE      TRUE
#>  [91,]   TRUE       TRUE      TRUE
#>  [92,]  FALSE       TRUE     FALSE
#>  [93,]   TRUE       TRUE      TRUE
#>  [94,]   TRUE       TRUE      TRUE
#>  [95,]   TRUE       TRUE      TRUE
#>  [96,]   TRUE       TRUE      TRUE
#>  [97,]   TRUE       TRUE      TRUE
#>  [98,]  FALSE       TRUE     FALSE
#>  [99,]   TRUE       TRUE      TRUE
#> [100,]   TRUE       TRUE      TRUE
#> [101,]  FALSE       TRUE     FALSE
#> [102,]   TRUE       TRUE      TRUE
#> [103,]  FALSE       TRUE     FALSE
#> [104,]  FALSE       TRUE     FALSE
#> [105,]  FALSE       TRUE     FALSE
#> [106,]  FALSE      FALSE     FALSE
#> [107,]   TRUE       TRUE      TRUE
#> [108,]  FALSE       TRUE     FALSE
#> [109,]  FALSE       TRUE     FALSE
#> [110,]  FALSE       TRUE     FALSE
#> [111,]  FALSE       TRUE     FALSE
#> [112,]  FALSE       TRUE     FALSE
#> [113,]  FALSE       TRUE     FALSE
#> [114,]   TRUE       TRUE      TRUE
#> [115,]   TRUE       TRUE      TRUE
#> [116,]  FALSE       TRUE     FALSE
#> [117,]  FALSE       TRUE     FALSE
#> [118,]  FALSE      FALSE     FALSE
#> [119,]  FALSE      FALSE     FALSE
#> [120,]  FALSE       TRUE     FALSE
#> [121,]  FALSE       TRUE     FALSE
#> [122,]   TRUE       TRUE      TRUE
#> [123,]  FALSE      FALSE     FALSE
#> [124,]  FALSE       TRUE     FALSE
#> [125,]  FALSE       TRUE     FALSE
#> [126,]  FALSE       TRUE     FALSE
#> [127,]  FALSE       TRUE     FALSE
#> [128,]  FALSE       TRUE     FALSE
#> [129,]  FALSE       TRUE     FALSE
#> [130,]  FALSE       TRUE     FALSE
#> [131,]  FALSE       TRUE     FALSE
#> [132,]  FALSE      FALSE     FALSE
#> [133,]  FALSE       TRUE     FALSE
#> [134,]  FALSE       TRUE     FALSE
#> [135,]  FALSE       TRUE     FALSE
#> [136,]  FALSE      FALSE     FALSE
#> [137,]  FALSE       TRUE     FALSE
#> [138,]  FALSE       TRUE     FALSE
#> [139,]  FALSE       TRUE     FALSE
#> [140,]  FALSE       TRUE     FALSE
#> [141,]  FALSE       TRUE     FALSE
#> [142,]  FALSE       TRUE     FALSE
#> [143,]   TRUE       TRUE      TRUE
#> [144,]  FALSE       TRUE     FALSE
#> [145,]  FALSE       TRUE     FALSE
#> [146,]  FALSE       TRUE     FALSE
#> [147,]  FALSE       TRUE     FALSE
#> [148,]  FALSE       TRUE     FALSE
#> [149,]  FALSE       TRUE     FALSE
#> [150,]  FALSE       TRUE     FALSE
attr(gps_csr, "csr_summary")
#>             treatment     csr_low  csr_high n_negative_matrix
#> setosa         setosa 0.100597981 0.9184742                86
#> versicolor versicolor 0.076518141 0.6392385                22
#> virginica   virginica 0.005007632 0.2601635                86
attr(gps_csr, "csr_data")
#>     Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
#> 1            5.1         3.5          1.4         0.2     setosa
#> 2            4.9         3.0          1.4         0.2     setosa
#> 5            5.0         3.6          1.4         0.2     setosa
#> 6            5.4         3.9          1.7         0.4     setosa
#> 8            5.0         3.4          1.5         0.2     setosa
#> 10           4.9         3.1          1.5         0.1     setosa
#> 11           5.4         3.7          1.5         0.2     setosa
#> 15           5.8         4.0          1.2         0.2     setosa
#> 16           5.7         4.4          1.5         0.4     setosa
#> 17           5.4         3.9          1.3         0.4     setosa
#> 18           5.1         3.5          1.4         0.3     setosa
#> 19           5.7         3.8          1.7         0.3     setosa
#> 20           5.1         3.8          1.5         0.3     setosa
#> 21           5.4         3.4          1.7         0.2     setosa
#> 22           5.1         3.7          1.5         0.4     setosa
#> 24           5.1         3.3          1.7         0.5     setosa
#> 26           5.0         3.0          1.6         0.2     setosa
#> 27           5.0         3.4          1.6         0.4     setosa
#> 28           5.2         3.5          1.5         0.2     setosa
#> 29           5.2         3.4          1.4         0.2     setosa
#> 32           5.4         3.4          1.5         0.4     setosa
#> 33           5.2         4.1          1.5         0.1     setosa
#> 34           5.5         4.2          1.4         0.2     setosa
#> 35           4.9         3.1          1.5         0.2     setosa
#> 36           5.0         3.2          1.2         0.2     setosa
#> 37           5.5         3.5          1.3         0.2     setosa
#> 38           4.9         3.6          1.4         0.1     setosa
#> 40           5.1         3.4          1.5         0.2     setosa
#> 41           5.0         3.5          1.3         0.3     setosa
#> 44           5.0         3.5          1.6         0.6     setosa
#> 45           5.1         3.8          1.9         0.4     setosa
#> 47           5.1         3.8          1.6         0.2     setosa
#> 49           5.3         3.7          1.5         0.2     setosa
#> 50           5.0         3.3          1.4         0.2     setosa
#> 54           5.5         2.3          4.0         1.3 versicolor
#> 56           5.7         2.8          4.5         1.3 versicolor
#> 58           4.9         2.4          3.3         1.0 versicolor
#> 60           5.2         2.7          3.9         1.4 versicolor
#> 61           5.0         2.0          3.5         1.0 versicolor
#> 65           5.6         2.9          3.6         1.3 versicolor
#> 67           5.6         3.0          4.5         1.5 versicolor
#> 68           5.8         2.7          4.1         1.0 versicolor
#> 70           5.6         2.5          3.9         1.1 versicolor
#> 80           5.7         2.6          3.5         1.0 versicolor
#> 81           5.5         2.4          3.8         1.1 versicolor
#> 82           5.5         2.4          3.7         1.0 versicolor
#> 83           5.8         2.7          3.9         1.2 versicolor
#> 85           5.4         3.0          4.5         1.5 versicolor
#> 89           5.6         3.0          4.1         1.3 versicolor
#> 90           5.5         2.5          4.0         1.3 versicolor
#> 91           5.5         2.6          4.4         1.2 versicolor
#> 93           5.8         2.6          4.0         1.2 versicolor
#> 94           5.0         2.3          3.3         1.0 versicolor
#> 95           5.6         2.7          4.2         1.3 versicolor
#> 96           5.7         3.0          4.2         1.2 versicolor
#> 97           5.7         2.9          4.2         1.3 versicolor
#> 99           5.1         2.5          3.0         1.1 versicolor
#> 100          5.7         2.8          4.1         1.3 versicolor
#> 102          5.8         2.7          5.1         1.9  virginica
#> 107          4.9         2.5          4.5         1.7  virginica
#> 114          5.7         2.5          5.0         2.0  virginica
#> 115          5.8         2.8          5.1         2.4  virginica
#> 122          5.6         2.8          4.9         2.0  virginica
#> 143          5.8         2.7          5.1         1.9  virginica
```
