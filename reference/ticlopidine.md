# Ticlopidine hydrochloride data from a crossover design

Original data were collected in Marzo et. al. (2002) to assess the
bioequivalence of a new formulation of ticlopidine hydrochloride with
the formulation that was marketed at that time. The original study
involved 24 healthy male volunteers who received both formulations in
the form of a tablet containing 250 mg of active ingredient in a 2 by 2
by 2 crossover design. The dataset contains differences between the two
formulations for each pair of pharmacokinetic outcomes after applying
the logarithmic transformation, where evident outliers were removed
bringing the available sample size to 20.

## Usage

``` r
data(ticlopidine)
```

## Format

A \`data.frame\` with 20 rows and 4 columns:

- t_half:

  Differences for the log-transformed elimination half-life.

- AUC:

  Differences for the log-transformed area under the concentration-time
  curve from time zero to the last measurable concentration.

- AUC_inf:

  Differences for the log-transformed area under the concentration-time
  curve from time zero to infinity.

- C_max:

  Differences for the log-transformed maximum plasma concentration.

## References

Marzo A. et al. "Bioequivalence of ticlopidine hydrochloride
administered in single dose to healthy volunteers", Pharmacological
Research, 2002.

Pallmann P. and Jaki T. "Simultaneous confidence regions for
multivariate bioequivalence", Statistics in Medicine, 2017.

Boulaguiem Y. et al. "Multivariate adjustments for average equivalence
testing", Statistics in Medicine, 2025.

## Examples

``` r
data(ticlopidine)
n <- nrow(ticlopidine)
nu <- n - 1
theta <- apply(ticlopidine,2,mean)
Sigma <- cov(ticlopidine)/n
```
