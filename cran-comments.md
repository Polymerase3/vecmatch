## R CMD check results

0 errors | 0 warnings | 1 note


## Minor release
In this version, I have:

* Added `optimize_gps()`, `make_opt_args()`, and `select_opt()` to support a new GPS-optimization workflow.
* Modified `csregion()` to allow reestimating the GPS after excluding observations.
* Fixed factor handling in `raincloud()` and `mosaic()`, now allowing custom facet ordering via releveling.
* Added SMD and p-value labels to `raincloud()`.
* Updated the `raincloud()` legend to explicitly include group names and the number of observations.
* Added a vignette.
