# vecmatch 1.2.0

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
