# Log-transformed cutaneous delivery of "Molecule X" on human skin as measured by two operators

Original data were collected using the cutaneous biodistribution method
described in Quartier et. al. (2019), and represents cutaneous delivery
of "Molecule X" (in ng/cm^2) on human abdominal skin. The study in Wu et
al. (2025) compared the reproducibility of an identical experimental
protocol performed by two different operators. The dataset contains
measurements from 6 independent skin samples per operator, after a
log-transformation.

## Usage

``` r
data(biodistribution)
```

## Format

A \`data.frame\` with 6 rows and 2 columns:

- X:

  Log-transformed measurements from the "reference" operator.

- Y:

  Log-transformed measurements from the "target" operator.

## References

Quartier J. et al. "Cutaneous Biodistribution: A High-Resolution
Methodology to Assess Bioequivalence in Topical Skin Delivery",
Pharmaceutics, 2019.

Wu J. et al. "Bridging the gap between experimental burden and
statistical power for quantiles equivalence testing", arXiv, 2025.

## Examples

``` r
data(biodistribution)
x_bar <- mean(biodistribution$X)
x_sd <- sd(biodistribution$X)
y_bar <- mean(biodistribution$Y)
y_sd <- sd(biodistribution$Y)
```
