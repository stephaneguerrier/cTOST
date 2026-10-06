# Comparison of a Corrective Procedure to the results of the Two One-Sided Tests (TOST) in Univariate Settings

This function renders a comparison of the alpha-TOST or the delta-TOST
outputs obtained with the function \`ctost\` to the TOST output obtained
with \`tost\`.

## Usage

``` r
compare_to_tost(x, ticks = 30, rn = 5)
```

## Arguments

- x:

  A `tost` object, which is the output of one of the function:
  \`ctost\`.

- ticks:

  an integer indicating the number of segments that will be printed to
  represent the confidence intervals.

- rn:

  integer indicating the number of decimals places to be used (see
  function \`round\`) for the printed results.

## Value

Prints a comparison between the TOST results (i.e., output of \`tost\`)
and either the alpha-TOST or the delta-TOST results.

## Examples

``` r
data(skin)
theta_hat = diff(apply(skin, 2, mean))
nu = nrow(skin) - 1
sig_hat = var(apply(skin, 1, diff)) / nu

# alpha-TOST
atost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
              alpha = 0.05, delta = log(1.25), method = "alpha")
compare_to_tost(atost)
#> TOST:       
#> ✖ Can't accept (bio)equivalence
#> alpha-TOST: 
#> ✔ Accept (bio)equivalence
#> 
#> Equiv. Region:  |---------------0---------------|  
#> TOST:            (---------------x----------------)
#> alpha-TOST:        (-------------x--------------)  
#> 
#>                  CI - low      CI - high
#> TOST:            -0.21174       0.25715
#> alpha-TOST:      -0.17654       0.22194
#> 
#> Equiv. lim. = +/- 0.22314
```
