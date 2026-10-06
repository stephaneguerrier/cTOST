# Core quantile equivalence testing procedure

This is the core the qTOST procedure to obtain CI, etc to test both a
single or multiple quantiles. It takes the test statistic (\`theta\`),
its standard error (\`sigma\`), and other parameters to perform the two
one-sided tests for equivalence. This is a low-level function; users
should typically use the main \`qtost\` wrapper function.

## Usage

``` r
qtost_core(
  theta,
  sigma,
  pi_x,
  delta_l,
  delta_u,
  alpha = 0.05,
  corrected_alpha = NULL
)
```

## Arguments

- theta:

  The calculated test statistic, \\\theta\\.

- sigma:

  The standard error of the test statistic, \\\sigma\_\theta\\.

- pi_x:

  A numeric scalar or vector representing the quantile(s) of interest in
  the reference population (\\X\\).

- delta_l:

  A numeric scalar or vector for the lower equivalence margin(s).

- delta_u:

  A numeric scalar or vector for the upper equivalence margin(s).

- alpha:

  The nominal significance level for the test (e.g., 0.05).

- corrected_alpha:

  An optional corrected significance level for \`alpha-qTOST\`.

## Value

An object of class \`qtost\` or \`mqtost\` containing the results of the
test. This list includes:

- \`decision\`: A boolean value (\`TRUE\` for equivalence, \`FALSE\`
  otherwise).

- \`ci\`: The confidence interval for the estimated quantile in
  population Y, \\\hat{\pi}\_y\\.

- \`pi_y_hat\`: The point estimate of the quantile in population Y,
  calculated as \\\Phi(\theta)\\.

- \`theta\`: The value of the test statistic \\\theta\\.

- \`sigma\`: The standard error of \\\theta\\.

- \`ci_theta\`: The confidence interval for \\\theta\\.

- \`alpha\`: The significance level used for the test.

- \`pi_x\`: The reference quantile(s).

- \`delta\`: The original equivalence margin.

- \`eq_region\`: The equivalence region defined by \`delta_l\` and
  \`delta_u\`.

- \`method\`: The name of the method used ("qTOST" or "alpha-qTOST").

- \`setting\`: The setting used ("single" or "multiple" quantiles).
