# vecmatch (development version)

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
