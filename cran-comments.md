## Resubmission

This is a resubmission of CLRtools.

Changes in this version:

- Fixed the NOTE regarding `id` being used as a global variable.
- Reduced vignette build time by precomputing Stan model results.
- Stored calibrated models in inst/extdata and loaded them using system.file().
