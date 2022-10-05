## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Changes after CRAN email for resubmission

* Changes to DESCRIPTION to add cited papers, Authors@R field

* Changes to CODE, DOCUMENTATION, and VIGNETTE, `T` into `TRUE` and `F` into `FALSE`, and `warning()` instead of `cat()`

* Changes to DOCUMENTATION, added VALUE field

* After running `R CMD check` the only NOTE is

  New Submission
  
  Possibly misspelled words in DESCRIPTION:
     Bierens (21:49)
     Escanciano (22:10)
     Geng (28:16)
     Gozalo (23:10)
     Lavergne (24:8, 25:8)
     Patilea (24:24, 25:24)
     Zheng (29:14)
  
* The package `AER` is unavailable on Linux and the application of the vignette of `SpeTestNP` uses `AER` thus `SpeTestNP` cannot be built on Linux. Maybe this can be bypassed.