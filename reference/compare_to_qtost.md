# Comparison of a Corrective Procedure to the Results of the Quantile Two One-Sided Tests (qTOST) in Single Quantile Setting

This function renders a comparison of the qTOST or the alpha-qTOST
outputs obtained with the function \`qtost\`.

## Usage

``` r
compare_to_qtost(x, ticks = 30, rn = 5)
```

## Arguments

- x:

  A `qtost` object, which is the output of one of the function:
  \`qtost\`.

- ticks:

  an integer indicating the number of segments that will be printed to
  represent the confidence intervals.

- rn:

  integer indicating the number of decimals places to be used (see
  function \`round\`) for the printed results.

## Value

Prints a comparison between the qTOST results (i.e., output of
\`qtost\`) and the alpha-qTOST results; returns `x` invisibly.

## Examples

``` r
# Using summary statistics from FDA label
x_bar_orig = 35.6
x_sd_orig = 16.7
n_x = 106
y_bar_orig = 41.6
y_sd_orig = 24.3
n_y = 14
x_bar = log(x_bar_orig^2 / sqrt(x_bar_orig^2 + x_sd_orig^2))
x_sd = sqrt(log(1 + (x_sd_orig^2 / x_bar_orig^2)))
y_bar = log(y_bar_orig^2 / sqrt(y_bar_orig^2 + y_sd_orig^2))
y_sd = sqrt(log(1 + (y_sd_orig^2 / y_bar_orig^2)))
x = list(mean=x_bar, sd=x_sd, n=n_x)
y = list(mean=y_bar, sd=y_sd, n=n_y)

# alpha-qTOST
aqtost <- qtost(x, y, pi_x = 0.8, delta = 0.15, method = "alpha")
compare_to_qtost(aqtost)
#> qTOST: ✖ Can't accept quantile (bio)equivalence
#> alpha-qTOST: ✖ Can't accept quantile (bio)equivalence
#> 
#> Equiv. Region:             |---------------------|
#> qTOST:          (----------x--------------)        
#> alpha-qTOST:     (-----------x------------)        
#> 
#>                CI - low      CI - high
#> qTOST:         0.50105       0.83712
#> alpha-qTOST:   0.51340       0.82937
#> 
#> Equiv. lim. = dw/up 0.65000 0.95000
```
