# Optimizing the Matching Process with a Random Search Algorithm

## Practical Example: Optimizing the Matching Process

Matching observations based on generalized propensity scores involves
tuning multiple hyperparameters. Running separate workflows with
different parameter combinations can be tedious, and the effects of some
parameters are not always predictable. To streamline this process,
vecmatch provides an automated optimization workflow using a random
search algorithm. The function
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md)
is implemented with multithreading to leverage computational resources
efficiently.

Step 1: Define the Formula, Data, and Optimization Space

In this example, we use the built-in `cancer` dataset and focus on two
predictors: the categorical variable `sex` and the continuous variable
`age`. We first specify the model formula. Note that `data` and
`formula` are fixed inputs; if you want to compare different formulas,
you must run the workflow separately for each.

Next, we define the search space for the hyperparameters. The helper
function
[`make_opt_args()`](https://polymerase3.github.io/vecmatch/reference/make_opt_args.md)
validates your inputs and computes the Cartesian product of all
specified values, reporting the total number of parameter combinations.

``` r

opt_args <- make_opt_args(
  data            = cancer,
  formula         = formula_cancer,
  reference       = c("control", "adenoma", "crc_benign", "crc_malignant"),
  gps_method      = c("m1", "m7", "m8"),
  matching_method = c("fullopt", "nnm"),
  caliper         = seq(0.01, 5, 0.01),
  cluster         = 1:3,
  ratio           = 1:3,
  min_controls    = 1:3,
  max_controls    = 1:3
)

opt_args
#> Optimization Argument Set (class: opt_args)
#> ---------------------------------------- 
#> gps_method      : m1, m7, m8
#> reference       : control, adenoma, crc_benign, crc_malignant
#> matching_method : fullopt, nnm
#> caliper         : [500 values]
#> order           : desc, asc, original, random
#> cluster         : 1, 2, 3
#> replace         : TRUE, FALSE
#> ties            : TRUE, FALSE
#> ratio           : 1, 2, 3
#> min_controls    : 1, 2, 3
#> max_controls    : 1, 2, 3
#> ---------------------------------------- 
#> Total combinations: 1512000
```

The [`print()`](https://rdrr.io/r/base/print.html) method for
[`make_opt_args()`](https://polymerase3.github.io/vecmatch/reference/make_opt_args.md)
provides a clear summary of the search space, including each
hyperparameter’s values and the total number of combinations.

## Step 2: Run the Optimizer

With the search space defined, we can launch the optimization. The
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md)
function performs a random search across the parameter grid and returns
a summary table containing key quality metrics for each tested
combination. You control the number of iterations via the `n_iter`
argument.

The function uses multithreading (via the `future` package) to
parallelize work. As a guideline, aim for at least 1000–1500 iterations
per core for reliable search coverage. Monitor your system’s memory
usage, since the parallel backend can consume substantial amounts of
RAM.

The function uses already registered backends.

By default,
[`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md)
preserves the global random seed. For reproducibility, set a seed before
calling the optimizer.

``` r

library(future)
library(doFuture)

## 1. Register future as the foreach backend
doFuture::registerDoFuture()

## 2. Choose parallel strategy
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(
  future::multisession,
  workers = 4
)

## 3. Seeding before calling optimize_gps()
set.seed(167894)
seed_before <- .Random.seed

opt_results <- optimize_gps(
  data     = cancer,
  formula  = formula_cancer,
  opt_args = opt_args,
  n_iter   = 6000
)

summary(opt_results)
#> Best Optimization Results by SMD Group
#> ======================================
#> 
#> -------------------------------------------------------
#> | smd_group  | unique_configs |  smd   | perc_matched |
#> -------------------------------------------------------
#> |0-0.05      |               2|   0.025|         79.73|
#> |0.05-0.10   |               1|   0.098|         80.73|
#> |0.10-0.15   |               7|   0.145|          91.2|
#> |0.15-0.20   |               3|   0.191|          96.1|
#> |0.20-0.25   |               1|   0.203|         93.52|
#> |0.25-0.30   |               2|   0.266|         79.01|
#> |0.25-0.30   |               2|   0.286|         79.01|
#> |>0.30       |               1|   0.343|         80.73|
#> -------------------------------------------------------
#> 
#> Optimization Summary
#> --------------------
#> Total combinations tested  : 6000
#> Total optimization time [s]: 362.29
```

We ran the optimization on a single core with `n_iter = 1500`; on our
test machine this required 362.29 seconds. Given the size of the
parameter grid, increasing `n_iter` would improve the search’s coverage,
but here we limited iterations to keep the vignette’s build time
reasonable.

When you print `opt_results`, it summarizes the entire search by
grouping parameter sets into bins defined by their maximum standardized
mean difference (SMD), and within each bin it highlights the combination
that achieves the highest proportion of matched observations.

## Step 3: Select Optimal Configurations

After optimization,
[`select_opt()`](https://polymerase3.github.io/vecmatch/reference/select_opt.md)
helps you choose parameter combinations that meet your specific
objectives. For example, you may wish to maximize the number of matched
samples for certain treatment groups while minimizing imbalance on key
covariates.

In this example, we aim to: \* Retain the largest possible number of
observations in the `adenoma` and `crc_malignant` groups. \* Minimize
the standardized mean difference (SMD) for the `age` variable.

We can achieve this by specifying arguments in
[`select_opt()`](https://polymerase3.github.io/vecmatch/reference/select_opt.md):

``` r

select_results <- select_opt(opt_results,
  perc_matched = c("adenoma", "crc_malignant"),
  smd_variables = "age",
  smd_type = "max"
)

summary(select_results)
#> 
#> Optimization Selection Summary
#> ====================
#> 
#> ------------------------------------------------------- 
#> | smd_group  | unique_configs |  smd   | perc_matched |
#> ------------------------------------------------------- 
#> |0-0.05      |               4|   0.047|        83.155|
#> |0.05-0.10   |               2|   0.094|        90.304|
#> |0.10-0.15   |               7|   0.145|        94.737|
#> |0.15-0.20   |               3|   0.191|        97.976|
#> |0.20-0.25   |               1|   0.203|        95.763|
#> |0.25-0.30   |               1|   0.266|        80.137|
#> |0.30-0.35   |               1|   0.343|        81.764|
#> |0.35-0.40   |               1|   0.371|         71.67|
#> |0.40-0.45   |               1|   0.445|        76.449|
#> |0.45-0.50   |               1|   0.465|         77.39|
#> |>0.50       |               1|   0.541|        80.082|
#> -------------------------------------------------------
```

The output shows the SMD bins and highlights the combination within each
bin that best meets our criteria. Suppose the configuration in the
`0.10–0.15` SMD bin offers a desirable balance of matched samples; we
can extract its parameters for manual refitting.

## Step 4: Refit the Optimized Model

To inspect the matched dataset and detailed balance summaries, we rerun
the standard `vecmatch` workflow using the selected parameters.

1.Select the subset with highest percentage of matched samples for
`smd_group` of interest and rematch the original data using the standard
`vecmatch` workflow. At the end you get an output object of class
`Matched`, just like in the standard workflow:

``` r

# estimating gps
final_match <- run_selected_matching(select_results,
  cancer,
  formula_cancer,
  smd_group = "0.05-0.10"
)

summary(final_match)
#> 
#> Treatment          | N before           | N after            | Retained (%)       
#> -------------------------------------------------------------------------------- 
#> adenoma            | 374                | 306                | 81.8%              
#> control            | 315                | 261                | 82.9%              
#> crc_benign         | 273                | 241                | 88.3%              
#> crc_malignant      | 248                | 245                | 98.8%
```

3.  You can call `balqual` on the output object to evaluate balance and
    verify reproducibility at then end:

``` r

balqual(final_match,
  formula = formula_cancer,
  type = "smd",
  statistic = "max",
  round = 3,
  cutoffs = 0.2
)
#> 
#> Matching Quality Evaluation
#> ================================================================================ 
#> 
#> Count table for the treatment variable:
#> -------------------------------------------------- 
#> Treatment                 | Before     | After      
#> -------------------------------------------------- 
#> adenoma                   | 374        | 306        
#> control                   | 315        | 261        
#> crc_benign                | 273        | 241        
#> crc_malignant             | 248        | 245        
#> -------------------------------------------------- 
#> 
#> 
#> Matching summary statistics:
#> ---------------------------------------- 
#> Total n before matching:  1210 
#> Total n after matching:       1053 
#> % of matched observations:    87.02 %
#> Total  maximal   SMD value:   0.13 
#> 
#> 
#> Maximal values :
#> -------------------------------------------------------------------------------- 
#> Variable                  | Coef  | Before       | After        | Quality      
#> -------------------------------------------------------------------------------- 
#> age                       | SMD   | 0.240        | 0.094        | Balanced     
#> sexF                      | SMD   | 0.138        | 0.116        | Balanced     
#> sexM                      | SMD   | 0.138        | 0.116        | Balanced     
#> age:sexF                  | SMD   | 0.143        | 0.130        | Balanced     
#> age:sexM                  | SMD   | 0.156        | 0.126        | Balanced     
#> --------------------------------------------------------------------------------

all.equal(seed_before, .Random.seed)
#> [1] TRUE
```
