## Submission

This is a new submission.

## Test environments

- Linux (NSCC ASPIRE2A), R 4.5.3: `R CMD check --as-cran` (2026-11-12: 0 errors, 0 warnings; the notes were from
  the build environment: no internet for the time check, and a folder now in .Rbuildignore)
- win-builder (R-devel) and the macOS builder: to be run before submission

## R CMD check results

0 errors | 0 warnings | 1 note

- New submission.

## Notes for the reviewers

- SummarizedExperiment and VariantAnnotation (Bioconductor) are in Suggests. They are used only by `from_proactiv()` and
  by reading VCF files, after `requireNamespace()`; the tests skip without them.
- The package is the implementation of the method of a manuscript under review; the reference will be added to
  DESCRIPTION once it is published.
