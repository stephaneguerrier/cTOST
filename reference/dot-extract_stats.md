# Extract Summary Statistics

An internal helper function that processes input data to ensure it is in
the correct summary statistics format (mean, standard deviation, and
sample size) for equivalence testing on multiple quantiles (see the
\`qtost\` function). If a numeric vector is provided, it calculates
these statistics. If a list is provided, it validates that the required
statistics are present.

## Usage

``` r
.extract_stats(data, group_name)
```

## Arguments

- data:

  A numeric vector, or a list containing the named elements \`mean\`,
  \`sd\`, and \`n\`.

- group_name:

  A character string representing the name of the group (e.g., "x" or
  "y"), used for constructing informative error messages.

## Value

A list containing three named elements:

- \`mean\`: The mean of the data.

- \`sd\`: The standard deviation of the data.

- \`n\`: The number of non-missing observations.
