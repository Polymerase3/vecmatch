## R CMD check results

0 errors | 0 warnings | 1 note


## vecmatch 1.3.0

* Refactored/unified S3 classes and added helpers for inspecting `vecmatch` objects.
* Added `get_select_params()` and `run_selected_matching()` for streamlined re-estimation.
* Reduced/cleaned dependencies and improved handling of suggested packages.
* Fixed `raincloud()` facet label bug.
* Removed backend handling from `optimize_gps()` (backend must now be registered externally).
* Updated optimization vignette, `cancer` dataset typo, README badges, examples (`\donttest{}`), and added reproducibility tests.
