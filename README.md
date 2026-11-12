# BOHR

**BOHR** (beta one-zero hybrid regression) is a regression for proportions in [0, 1] that contain exact 0s and 1s. The response is modelled with a zero-and-one-inflated beta distribution, and a predictor's effect is tested separately on three components:
- **p0**: the probability of 0;
- **p1**: the probability of 1;
- **μ**: the mean of the values in between.

Its first application is alternative promoter usage, where p0 is the probability that a promoter is silent, p1 the probability that it is its gene's only active promoter, and μ its mean usage otherwise. The genomics layer maps genetic effects on these components (cis-QTL mapping).

The repository holds the R package at its root, and:
- `workflow/parity/`: scripts that compare the package with the stored results of the paper's original pipeline on the GTEx data;
- `paper/` (to come): the analysis code of the accompanying manuscript, organised by Methods section.

Neither folder is part of the installed package.

## Installation

```r
# install.packages("remotes")
remotes::install_github("boxiangliulab/BOHR")
```

While the repository is private, `install_github()` needs a GitHub personal access token with read access to it, set as the environment variable `GITHUB_PAT`.

## License

MIT
