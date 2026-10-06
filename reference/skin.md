# Log-transformed cutaneous delivery of econazole (ECZ) from bioequivalent products on porcine skin

Original data were collected using the cutaneous biodistribution method
described in Quartier et. al. (2019), and represents cutaneous delivery
of econazole nitrate (ECZ in ng/cm^2) on porcine skin from a reference
medicinal product and an approved bioequivalent product. The dataset
contains 17 pairs of comparable porcine skin samples on which
measurement of ECZ deposition was gathered, and log transformed, using
both creams.

## Usage

``` r
data(skin)
```

## Format

A \`data.frame\` with 17 rows and 2 columns:

- Reference:

  Log-transformed econazole nitrate delivery for the reference product.

- Generic:

  Log-transformed econazole nitrate delivery for the generic
  bioequivalent product.

## References

Quartier J. et al. "Cutaneous Biodistribution: A High-Resolution
Methodology to Assess Bioequivalence in Topical Skin Delivery",
Pharmaceutics, 2019.

Boulaguiem Y. et al. "Finite Sample Adjustments for Average Equivalence
Testing", Statistics in Medicine, 2024.

## Examples

``` r
data(skin)
theta <- diff(apply(skin,2,mean))
n <- nrow(skin)
nu <- n - 1
sigma <- var(apply(skin,1,diff)) / n
```
