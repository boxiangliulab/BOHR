# BOHR

**BOHR** (beta one-zero hybrid regression) maps genetic and trait effects on alternative promoter usage. Relative promoter activity is modelled with a zero-and-one-inflated beta distribution, and three components are tested separately:
- **p0**: the probability that a promoter is silent;
- **p1**: the probability that it is the only active promoter of its gene;
- **μ**: the mean usage when it is neither.

The repository is under construction. It will hold the R package, at the repository root, and the analysis workflow of the accompanying manuscript, in `workflow/`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YuhaoDoong/BOHR")
```

## License

MIT
