## Test environments

* Local macOS, R release
* GitHub Actions: macOS, Windows, Ubuntu; R release, oldrel, and devel
* R-hub and sanitizer checks: planned before any CRAN submission

## R CMD check results

Development branch checks are tracked in GitHub Actions.

## Notes

This development branch adds experimental deterministic QMC backends while
retaining the existing Fortran/C implementations as the reference path.
