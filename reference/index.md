# Package index

## Average Equivalence Testing

Functions for testing equivalence of means in univariate and
multivariate settings using finite sample corrections (alpha-TOST,
delta-TOST, optimal cTOST).

- [`ctost()`](https://stephaneguerrier.github.io/cTOST/reference/ctost.md)
  : Finite Sample Adjustment for Average (Bio)Equivalence Assessment
- [`compare_to_tost()`](https://stephaneguerrier.github.io/cTOST/reference/compare_to_tost.md)
  : Comparison of a Corrective Procedure to the results of the Two
  One-Sided Tests (TOST) in Univariate Settings

## Quantile Equivalence Testing

Functions for testing equivalence at specific quantiles or across
multiple quantiles.

- [`qtost()`](https://stephaneguerrier.github.io/cTOST/reference/qtost.md)
  : Quantile equivalence testing procedures
- [`compare_to_qtost()`](https://stephaneguerrier.github.io/cTOST/reference/compare_to_qtost.md)
  : Comparison of a Corrective Procedure to the Results of the Quantile
  Two One-Sided Tests (qTOST) in Single Quantile Setting

## Differential Privacy Equivalence Testing

Functions for equivalence testing under differential privacy
constraints.

- [`prop_test_equiv_dp()`](https://stephaneguerrier.github.io/cTOST/reference/prop_test_equiv_dp.md)
  : Differentially Private TOST for Proportion Equivalence Testing
- [`tost_equiv_dp()`](https://stephaneguerrier.github.io/cTOST/reference/tost_equiv_dp.md)
  : Differentially Private TOST for Mean Equivalence Testing

## Visualization and Output

Print and plot methods for displaying test results.

- [`plot(`*`<tost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/plot.tost.md)
  [`plot(`*`<mtost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/plot.tost.md)
  [`plot(`*`<qtost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/plot.tost.md)
  [`plot(`*`<mqtost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/plot.tost.md)
  : Plot confidence intervals for one or more \`tost\`, \`mtost\`,
  \`qtost\` or \`mqtost\` objects
- [`print(`*`<tost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/print.tost.md)
  : Print Results of (Bio)Equivalence Assessment in Univariate Settings
- [`print(`*`<mtost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/print.mtost.md)
  : Print Results of (Bio)Equivalence Assessment in Multivariate
  Settings
- [`print(`*`<qtost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/print.qtost.md)
  : Print Results of (Bio)Equivalence testing in Single Quantile Setting
- [`print(`*`<mqtost>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/print.mqtost.md)
  : Print Results of (Bio)Equivalence testing in Two Quantiles Setting
- [`print(`*`<prop_test_dp>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/print.prop_test_dp.md)
  : Print Results of DP-TOST for Proportions
- [`print(`*`<tost_dp>`*`)`](https://stephaneguerrier.github.io/cTOST/reference/print.tost_dp.md)
  : Print method for DP-TOST mean test

## Datasets

Real datasets for bioequivalence studies and examples.

- [`skin`](https://stephaneguerrier.github.io/cTOST/reference/skin.md) :
  Log-transformed cutaneous delivery of econazole (ECZ) from
  bioequivalent products on porcine skin
- [`ticlopidine`](https://stephaneguerrier.github.io/cTOST/reference/ticlopidine.md)
  : Ticlopidine hydrochloride data from a crossover design
- [`skin_mvt`](https://stephaneguerrier.github.io/cTOST/reference/skin_mvt.md)
  : Multivariate measurements for log-transformed cutaneous delivery of
  econazole (ECZ) from bioequivalent products on porcine skin
- [`biodistribution`](https://stephaneguerrier.github.io/cTOST/reference/biodistribution.md)
  : Log-transformed cutaneous delivery of "Molecule X" on human skin as
  measured by two operators
