## R CMD check results

0 errors | 0 warnings | 0 notes


## vecmatch 1.4.0

* `balqual()` gained the `type` values `"desc_full"` and `"desc_reduced"`, reporting descriptive statistics per treatment level for the unmatched and matched dataset.
* Fixed column selection in `match_gps()` for `method = "nnm"`, which failed with `undefined columns selected` for any `reference` other than the first gps column.
* `match_gps()` now honours `max_controls` for `method = "fullopt"`; it was passed under a name `optmatch::fullmatch()` does not accept and carried the value of `caliper`.
* `csregion()` now stops with an informative message when the common support region is empty or a treatment group loses all observations, and its refitting fallback warning is no longer malformed.
* `match_gps()` no longer warns for method-specific tuning parameters left at their defaults.
* Fixed `cutoffs` in `balqual()` being matched to metrics by position rather than by name.
* `print.quality()` now returns its argument invisibly, as documented for `print()` methods.
* `optimize_gps()` no longer aborts when GPS estimation fails for the first combination of the search space.
* Added a `pkgdown` website at <https://polymerase3.github.io/vecmatch/>.
