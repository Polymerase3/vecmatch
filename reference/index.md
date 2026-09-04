# Package index

## Assessing initial imbalance

Visualize the covariate imbalances in the raw dataset before any
propensity score estimation.

- [`raincloud()`](https://polymerase3.github.io/vecmatch/reference/raincloud.md)
  : Examine the Imbalance of Continuous Covariates
- [`mosaic()`](https://polymerase3.github.io/vecmatch/reference/mosaic.md)
  : Plot the Distribution of Categorical Covariates

## Estimating generalized propensity scores

Estimate the treatment allocation probabilities and define the common
support region.

- [`estimate_gps()`](https://polymerase3.github.io/vecmatch/reference/estimate_gps.md)
  : Calculate Treatment Allocation Probabilities
- [`csregion()`](https://polymerase3.github.io/vecmatch/reference/csregion.md)
  : Filter the Data Based on Common Support Region

## Matching

Match the observations across multiple treatment groups based on the
estimated generalized propensity scores.

- [`match_gps()`](https://polymerase3.github.io/vecmatch/reference/match_gps.md)
  : Match the Data Based on Generalized Propensity Scores

## Evaluating matching quality

Assess the balance of the covariates and the descriptive statistics
before and after matching.

- [`balqual()`](https://polymerase3.github.io/vecmatch/reference/balqual.md)
  : Evaluate Matching Quality

## Optimizing the matching process

Search the parameter space of the estimation and matching functions and
rerun the pipeline for the selected configurations.

- [`optimize_gps()`](https://polymerase3.github.io/vecmatch/reference/optimize_gps.md)
  : Optimize the Matching Process via Random Search
- [`make_opt_args()`](https://polymerase3.github.io/vecmatch/reference/make_opt_args.md)
  : Define the Optimization Parameter Space for Matching
- [`select_opt()`](https://polymerase3.github.io/vecmatch/reference/select_opt.md)
  : Select Optimal Parameter Combinations from Optimization Results
- [`get_select_params()`](https://polymerase3.github.io/vecmatch/reference/get_select_params.md)
  : Extract Parameter Grid for Selected Configurations
- [`run_selected_matching()`](https://polymerase3.github.io/vecmatch/reference/run_selected_matching.md)
  : Rerun GPS Estimation and Matching for a Selected Configuration

## Datasets

- [`cancer`](https://polymerase3.github.io/vecmatch/reference/cancer.md)
  : Patients with Colorectal Cancer and Adenoma

## Package overview

- [`vecmatch`](https://polymerase3.github.io/vecmatch/reference/vecmatch-package.md)
  [`vecmatch-package`](https://polymerase3.github.io/vecmatch/reference/vecmatch-package.md)
  : vecmatch: Vector Matching for Generalized Propensity Scores
