# Multivariate measurements for log-transformed cutaneous delivery of econazole (ECZ) from bioequivalent products on porcine skin

Original data were collected using the cutaneous biodistribution method
described in Quartier et. al. (2019), and represents cutaneous delivery
of econazole nitrate (ECZ in ng/cm^2) on porcine skin from a reference
medicinal product and an approved bioequivalent product across multiple
skin layers (from the surface to a depth of approximately 800 \\\mu\\m).
The dataset contains differences of log-transformed measurements of ECZ
deposition for the two creams, across 12 comparable porcine skin samples
and four anatomical regions.

## Usage

``` r
data(skin_mvt)
```

## Format

A \`data.frame\` with 12 rows and 4 columns:

- \`stratum corneum\`:

  Differences of log-transformed econazole nitrate delivery for the two
  creams across the stratum corneum (0-20 \\\mu\\m).

- \`viable epidermis\`:

  Differences of log-transformed econazole nitrate delivery for the two
  creams across the viable epidermis (20-160 \\\mu\\m).

- \`upper dermis\`:

  Differences of log-transformed econazole nitrate delivery for the two
  creams across the upper dermis (160-400 \\\mu\\m).

- \`lower dermis\`:

  Differences of log-transformed econazole nitrate delivery for the two
  creams across the lower dermis (400-800 \\\mu\\m).

## References

Quartier J. et al. "Cutaneous Biodistribution: A High-Resolution
Methodology to Assess Bioequivalence in Topical Skin Delivery",
Pharmaceutics, 2019.

Insolia L. et al. "Bioequivalence Assessment for Locally Acting Drugs: A
Framework for Feasible and Efficient Evaluation", arXiv, 2025.

## Examples

``` r
data(skin_mvt)
n <- nrow(skin_mvt)
nu <- n - 1
theta <- apply(skin_mvt,2,mean)
Sigma <- cov(skin_mvt) / n
```
