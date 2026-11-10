# BOHR

**BOHR** (beta one-zero hybrid regression) is a regression for proportions in [0, 1] that contain exact 0s and 1s. The response is modelled with a zero-and-one-inflated beta distribution, and a predictor's effect is tested separately on three components:
- **p0**: the probability of 0;
- **p1**: the probability of 1;
- **μ**: the mean of the values in between.

Its first application is alternative promoter usage, where p0 is the probability that a promoter is silent, p1 the probability that it is its gene's only active promoter, and μ its mean usage otherwise. The genomics layer maps genetic effects on these components (cis-QTL mapping).

The repository is under construction. It will hold the R package, at the repository root, and the analysis workflow of the accompanying manuscript, in `workflow/`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YuhaoDoong/BOHR")
```

## License

MIT
