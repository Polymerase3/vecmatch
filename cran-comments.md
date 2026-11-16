## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no reverse dependencies.


## vecmatch 1.4.1

This is a patch release fixing a bug in `summary()` for `csr` objects.

* `summary()` on a `csr` object reported wrong group sizes before filtering and wrong per-group exclusion counts, because both were derived from per-column counts of observations outside the GPS bounds. Both are now taken from the treatment variable itself, so they add up to the sample size and to the number of removed observations.
* The per-treatment bounds table now labels its overlapping count `Outside bounds` instead of `Number excluded`, and the corresponding field returned by `summary()` was renamed from `n_excluded` to `n_outside_bounds`.
* The `csr_summary` attribute of a `csr` object gained the columns `n_before` and `n_excluded_group`.
