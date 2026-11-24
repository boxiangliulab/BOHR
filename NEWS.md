# BOHR 0.1.3

* `paper/` now holds the analysis code of the final version of the paper: the run without a junction-read floor,
  the skip rule SD < 0.005 with calibrated admission of the added tests, the boundary (80%) and read-depth (five reads)
  rules, tissue means on GTEx analysis-freeze samples, and an ancestry-matched LD reference for GWAS clumping. See
  "Version of the analysis" in `paper/README.md`. The whole-chain test of the earlier version (`paper/run_paper.sh`,
  `checks/`, `with_package/`, `test_subset/`) is in releases 0.1.1 and 0.1.2 only.
* The package code is unchanged from 0.1.0. Its defaults are the earlier skip rules; `paper/README.md` gives the
  `bohr_control()` settings of the final version.

# BOHR 0.1.2

* Same content as 0.1.1: the release was tagged on the commit of 0.1.1 before the changes of 0.1.3 were pushed.

# BOHR 0.1.1

* Tests: the bit-for-bit comparison with the stored engine output runs only with `BOHR_EXACT_PARITY=true`; otherwise
  values are compared with a relative tolerance of 1e-6. The package code is unchanged; its version number stayed 0.1.0.

# BOHR 0.1.0

* First release.
