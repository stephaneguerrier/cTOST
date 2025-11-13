#' Print Results of (Bio)Equivalence Assessment in Univariate Settings
#'
#' @param x      A \code{tost} object, which is the output of one of the following functions `tost` or `ctost`.
#' @param ticks  Number of ticks to print the confidence interval in the console.
#' @param rn     Number of digits to consider when printing the results.
#' @param ...    Further arguments to be passed to or from methods.
#' @return       Prints object.
#' @importFrom   cli cli_text col_green col_red
#'
#' @rdname print.tost
#'
#' @export
print.tost = function(x, ticks = 30, rn = 5, ...){

  if (x$decision){
    cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
  }else{
    cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
  }

  if (x$method == "delta-TOST"){
    lower_be = x$ci[1] > -x$corrected_delta
    upper_be = x$ci[2] < x$corrected_delta
    rg = range(c(x$ci, x$corrected_delta, -x$corrected_delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$corrected_delta, x$corrected_delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1
    std_fit_interval = round(ticks*(x$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center = round(ticks*(sum(x$ci)/2 - rg[1])/rg_delta) + 1
  }else{
    if (!(x$method %in% c("qTOST", "alpha-qTOST"))) {
      shift_ci = 0
    } else {
      shift_ci = x$pi_x
    }
    lower_be = x$ci[1] - shift_ci > -x$delta
    upper_be = x$ci[2] - shift_ci < x$delta
    rg = range(c(x$ci - shift_ci, x$delta, -x$delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$delta, x$delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1
    std_fit_interval = round(ticks*(x$ci - shift_ci - rg[1])/rg_delta) + 1
    std_fit_interval_center = round(ticks*(sum(x$ci - shift_ci)/2 - rg[1])/rg_delta) + 1
    std_fit_interval = round(ticks*(x$ci - shift_ci - rg[1])/rg_delta) + 1
    std_fit_interval_center = round(ticks*(sum(x$ci - shift_ci)/2 - rg[1])/rg_delta) + 1
  }

  if (x$method == "delta-TOST"){
    cat("Corr. Equiv. Region:  ")
  }else{
    cat("Equiv. Region:  ")
  }

  for (i in 1:(ticks+1)){
    if (i >= std_be_interval[1] && i <= std_be_interval[2]){
      if (i == std_be_interval[1]){
        cat(("|-"))
      }else{
        if (i ==  std_be_interval[2]){
          cat(("-|"))
        }else{
          if (i == std_zero){
            cat(("-0-"))
          }else{
            cat(("-"))
          }
        }
      }
    }else{
      cat(" ")
    }
  }

  cat("\n")

  if (x$method == "delta-TOST"){
    cat("      Estim. Inter.:  ")
  }else{
    cat("Estim. Inter.:  ")
  }

  for (i in 1:(ticks+1)){
    if (i >= std_fit_interval[1] && i <= std_fit_interval[2]){
      if (i == std_fit_interval[1]){
        if (i > std_be_interval[1] && i < std_be_interval[2]){
          cat(col_green("(-"))
        }else{
          if (lower_be){
            cat(col_green("(-"))
          }else{
            cat(col_red("(-"))
          }
        }

      }else{
        if (i ==  std_fit_interval[2]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("-)"))
          }else{
            if (upper_be){
              cat(col_green("-)"))
            }else{
              cat(col_red("-)"))
            }
          }

        }else{
          if (i == std_fit_interval_center){
            if (i >= std_be_interval[1] && i <= std_be_interval[2]){
              cat(col_green("-x-"))
            }else{
              cat(col_red("-x-"))
            }
          }else{
            if (i >= std_be_interval[1] && i <= std_be_interval[2]){
              cat(col_green("-"))
            }else{
              cat(col_red("-"))
            }
          }
        }
      }
    }else{
      cat(" ")
    }
  }
  cat("\n")
  cat("CI =  (")
  cat(format(round(x$ci[1], rn), nsmall = rn))
  cat(" ; ")
  cat(format(round(x$ci[2], rn), nsmall = rn))
  cat(")\n\n")

  cat("Method: ")
  cat(x$method)
  cat("\n")
  cat("alpha = ")
  cat(x$alpha)
  cat("; ")
  if (x$method %in% c("qTOST", "alpha-qTOST")) {
    cat("pi_x = ")
    cat(format(round(shift_ci, 2), nsmall = 2))
    cat("; ")
    cat("Equiv. lim. = +/- ")
    cat(format(round(x$delta, 2), nsmall = 2))
  } else {
    cat("Equiv. lim. = +/- ")
    cat(format(round(x$delta, rn), nsmall = rn))
  }
  cat("\n")
  if (x$method == "alpha-TOST" || x$method == "alpha-qTOST"){
    cat("Corrected alpha = ")
    cat(format(round(x$corrected_alpha, rn), nsmall = rn))
    cat("\n")
  }

  if (x$method == "delta-TOST"){
    cat("Corrected Equiv. lim. = +/- ")
    cat(format(round(x$corrected_delta, rn), nsmall = rn))
    cat("\n")
  }

  if (x$method == "cTOST"){
    cat("Estimated c(0) = ")
    cat(format(round(x$corrected_c, rn), nsmall = rn))
    cat("\n")

    cat("Finite sample correction: ")
    cat(x$correction)
    cat("\n")
    if (x$correction %in% c("bootstrap", "offline")){
      cat("Corrected alpha = ")
      cat(format(round(x$corrected_alpha, rn), nsmall = rn))
      cat("\n")
    }
  }
  if (x$method %in% c("qTOST", "alpha-qTOST")) {
    cat("Estimates: ")
    cat("pi_y = ")
    cat(format(round(x$pi_y_hat, rn), nsmall = rn))
    cat("; ")
    cat("Mean = ")
    cat(format(round(x$theta, rn), nsmall = rn))
    cat("; ")
    cat("Stand. dev. = ")
    cat(format(round(x$sigma, rn), nsmall = rn))
    cat("\n")
  } else {
    cat("Mean = ")
    cat(format(round(x$theta, rn), nsmall = rn))
    cat("; ")
    cat("Stand. dev. = ")
    cat(format(round(sqrt(x$sigma), rn), nsmall = rn))
    cat("; ")
    cat("df = ")
    cat(x$nu)
    cat("\n")
  }
}


#' Print Results of (Bio)Equivalence Assessment in Multivariate Settings
#'
#' @param x      A \code{mtost} object, which is the output of one of the following functions `tost` or `ctost`.
#' @param ticks  Number of ticks to print the confidence interval in the console.
#' @param rn     Number of digits to consider when printing the results.
#' @param ...    Further arguments to be passed to or from methods.
#' @return       Prints object.
#' @importFrom   cli cli_text col_green col_red
#'
#' @rdname print.mtost
#'
#' @export
print.mtost = function(x, ticks = 30, rn = 5, ...){
  p = length(x$decision)

  if (all(x$decision)){
    cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
  }else{
    cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
  }

  if (!(x$method %in% c("qTOST", "alpha-qTOST"))) {
    shift_ci = rep(0, p)
  } else {
    shift_ci = x$pi_x
  }
  rg = range(c(x$ci - shift_ci, x$delta, -x$delta))
  rg_delta = rg[2] - rg[1]
  std_be_interval = round(ticks*(c(-x$delta, x$delta) - rg[1])/rg_delta) + 1
  std_zero = round(-ticks*rg[1]/rg_delta) + 1

  lower_be = upper_be = std_fit_interval_center = rep(NA, p)
  std_fit_interval = matrix(NA, p, 2)
  for (i in 1:p){
    lower_be[i] = x$ci[i,1] - shift_ci[i] > -x$delta
    upper_be[i] = x$ci[i,2] - shift_ci[i] < x$delta
    std_fit_interval[i,] = round(ticks*(x$ci[i,] - shift_ci[i] - rg[1])/rg_delta) + 1
    std_fit_interval_center[i] = round(ticks*(sum(x$ci[i,] - shift_ci[i])/2 - rg[1])/rg_delta) + 1
  }

  if (x$method %in% c("qTOST", "alpha-qTOST")) {
    names_var = paste0("pi_x = ", x$pi_x)
  } else {
    names_var = colnames(x$sigma)
  }
  names_len = nchar(names_var)


  cat("Equiv. Region:  ")

  if (nchar("Equiv. Region:   ") < max(names_len)){
    cat(paste(rep(" ", max(names_len) - nchar("Equiv. Region:   ")), collapse = ""))
  }
  cat(" ")

  for (i in 1:(ticks+1)){
    if (i >= std_be_interval[1] && i <= std_be_interval[2]){
      if (i == std_be_interval[1]){
        cat(("|-"))
      }else{
        if (i ==  std_be_interval[2]){
          cat(("-|"))
        }else{
          if (i == std_zero){
            cat(("-0-"))
          }else{
            cat(("-"))
          }
        }
      }
    }else{
      cat(" ")
    }
  }

  cat("\n")
  for (h in 1:p){
    if (max(names_len) > names_len[h]){
      cat(names_var[h])
      cat(paste(rep(" ", max(names_len) - names_len[h]), collapse = ""))
      cat(" ")
    }else{
      cat(names_var[h])
      cat(" ")
    }

    if (nchar("Equiv. Region:   ") > max(names_len)){
      cat(paste(rep(" ", nchar("Equiv. Region:  ") - max(names_len)), collapse = ""))
    }


    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval[h,1] && i <= std_fit_interval[h,2]){
        if (i == std_fit_interval[h,1]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("(-"))
          }else{
            if (lower_be[h]){
              cat(col_green("(-"))
            }else{
              cat(col_red("(-"))
            }
          }

        }else{
          if (i ==  std_fit_interval[h,2]){
            if (i > std_be_interval[1] && i < std_be_interval[2]){
              cat(col_green("-)"))
            }else{
              if (upper_be[h]){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }

          }else{
            if (i == std_fit_interval_center[h]){
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
  }

  cat("\n")
  cat("CIs:")
  cat("\n")

  for (i in 1:p){
    if (max(names_len) > names_len[i]){
      cat(names_var[i])
      cat(paste(rep(" ", max(names_len) - names_len[i]), collapse = ""))
      cat("  ")
    }else{
      cat(names_var[i])
      cat("  ")
    }

    cat("(")
    cat(format(round(x$ci[i,1], rn), nsmall = rn))
    cat(" ; ")
    cat(format(round(x$ci[i,2], rn), nsmall = rn))
    cat(")     ")
    if (x$decision[i]){
      cli_text(col_green("{symbol$tick}"))
    }else{
      cli_text(col_red("{symbol$cross}"))
    }
  }
  cat("\n")

  cat("Method: ")
  cat(x$method)
  cat("\n")
  cat("alpha = ")
  cat(x$alpha)
  cat("; ")
  if (x$method %in% c("qTOST", "alpha-qTOST")) {
    cat("pi_x = (")
    cat(format(round(shift_ci[1], 2), nsmall = 2))
    cat(", ")
    cat(format(round(shift_ci[2], 2), nsmall = 2))
    cat(")")
    cat("; ")
    cat("Equiv. lim. = +/- ")
    cat(format(round(x$delta, 2), nsmall = 2))
  } else {
    cat("Equiv. lim. = +/- ")
    cat(format(round(x$delta, rn), nsmall = rn))
  }
  cat("\n")
  if (x$method == "alpha-TOST" || x$method == "alpha-qTOST"){
    cat("Corrected alpha = ")
    cat(format(round(x$corrected_alpha, rn), nsmall = rn))
    cat("\n")
  }
}

#' @title Comparison of a Corrective Procedure to the results of the Two One-Sided Tests (TOST) in Univariate Settings
#'
#' @description This function renders a comparison of the alpha-TOST or the delta-TOST outputs obtained with the function `ctost` to the TOST output obtained with `tost`.
#'
#' @param x A \code{tost} object, which is the output of one of the function: `ctost`.
#' @param ticks an integer indicating the number of segments that will be printed to represent the confidence intervals.
#' @param rn integer indicating the number of decimals places to be used (see function `round`) for the printed results.
#' @return Pints a comparison between the TOST results (i.e., output of `tost`) and either the alpha-TOST or the delta-TOST results.
#'
#' @importFrom cli cli_text col_green col_red
#'
#' @export
compare_to_tost = function(x, ticks = 30, rn = 5){
  result_tost = tost(theta = x$theta, sigma = sqrt(x$sigma), nu = x$nu,
                     alpha = x$alpha, delta = x$delta)

  if (!(x$method %in% c("alpha-TOST", "delta-TOST", "cTOST"))){
    stop("This method is not compatible")
  }

  if (x$method == "delta-TOST"){
    cat("TOST:       ")
    if (result_tost$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }

    cat("delta-TOST: ")
    if (x$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }
    cat("\n")

    lower_be_tost = x$ci[1] > -x$delta
    upper_be_tost = x$ci[2] < x$delta

    lower_be_dtost = x$ci[1] > -x$corrected_delta
    upper_be_dtost = x$ci[2] < x$corrected_delta

    rg = range(c(x$ci, result_tost$ci, x$delta, -x$delta, x$corrected_delta, -x$corrected_delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$delta, x$delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1

    std_be_interval_cor = round(ticks*(c(-x$corrected_delta, x$corrected_delta) - rg[1])/rg_delta) + 1
    std_zero_cor = round(-ticks*rg[1]/rg_delta) + 1

    std_fit_interval_tost = round(ticks*(result_tost$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center_tost = round(ticks*(sum(result_tost$ci)/2 - rg[1])/rg_delta) + 1

    cat("Stand. Equiv. Region:  ")

    for (i in 1:(ticks+1)){
      if (i >= std_be_interval[1] && i <= std_be_interval[2]){
        if (i == std_be_interval[1]){
          cat(("|-"))
        }else{
          if (i ==  std_be_interval[2]){
            cat(("-|"))
          }else{
            if (i == std_zero){
              cat(("-0-"))
            }else{
              cat(("-"))
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat(" Corr. Equiv. Region:  ")

    for (i in 1:(ticks+1)){
      if (i >= std_be_interval_cor[1] && i <= std_be_interval_cor[2]){
        if (i == std_be_interval_cor[1]){
          cat(("|-"))
        }else{
          if (i ==  std_be_interval_cor[2]){
            cat(("-|"))
          }else{
            if (i == std_zero_cor){
              cat(("-0-"))
            }else{
              cat(("-"))
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat("                TOST:  ")
    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_tost[1] && i <= std_fit_interval_tost[2]){
        if (i == std_fit_interval_tost[1]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("(-"))
          }else{
            if (lower_be_tost){
              cat(col_green("(-"))
            }else{
              cat(col_red("(-"))
            }
          }

        }else{
          if (i ==  std_fit_interval_tost[2]){
            if (i > std_be_interval[1] && i < std_be_interval[2]){
              cat(col_green("-)"))
            }else{
              if(upper_be_tost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }

          }else{
            if (i == std_fit_interval_center_tost){
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat("          delta-TOST:  ")
    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_tost[1] && i <= std_fit_interval_tost[2]){
        if (i == std_fit_interval_tost[1]){
          if (i > std_be_interval_cor[1] && i < std_be_interval_cor[2]){
            cat(col_green("(-"))
          }else{
            if(lower_be_dtost){
              cat(col_green("-)"))
            }else{
              cat(col_red("-)"))
            }
          }

        }else{
          if (i ==  std_fit_interval_tost[2]){
            if (i > std_be_interval_cor[1] && i < std_be_interval_cor[2]){
              cat(col_green("-)"))
            }else{
              if(upper_be_dtost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }

          }else{
            if (i == std_fit_interval_center_tost){
              if (i >= std_be_interval_cor[1] && i <= std_be_interval_cor[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval_cor[1] && i <= std_be_interval_cor[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }

    cat("\n")
    cat("\n")
    cat(" Standard Equiv. lim. = +/- ")
    cat(format(round(x$delta, rn), nsmall = rn))
    cat("\n")

    cat("Corrected Equiv. lim. = +/- ")
    cat(format(round(x$corrected_delta, rn), nsmall = rn))
    cat("\n")
  }

  if (x$method %in% c("alpha-TOST", "cTOST")){

    lower_be_atost = x$ci[1] > -x$delta
    upper_be_atost = x$ci[2] < x$delta
    lower_be_tost = result_tost$ci[1] > -x$delta
    upper_be_tost = result_tost$ci[2] < x$delta

    if (x$method == "alpha-TOST"){
      cat("TOST:       ")
    }else{
      cat("TOST:   ")
    }
    if (result_tost$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }


    if (x$method == "alpha-TOST"){
      cat("alpha-TOST: ")
    }else{
      cat("cTOST: ")
    }

    if (x$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }
    cat("\n")

    rg = range(c(x$ci, result_tost$ci, x$delta, -x$delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$delta, x$delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1
    std_fit_interval_atost = round(ticks*(x$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center_atost = round(ticks*(sum(x$ci)/2 - rg[1])/rg_delta) + 1
    std_fit_interval_tost = round(ticks*(result_tost$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center_tost = round(ticks*(sum(result_tost$ci)/2 - rg[1])/rg_delta) + 1

    cat("Equiv. Region:  ")

    for (i in 1:(ticks+1)){
      if (i >= std_be_interval[1] && i <= std_be_interval[2]){
        if (i == std_be_interval[1]){
          cat(("|-"))
        }else{
          if (i ==  std_be_interval[2]){
            cat(("-|"))
          }else{
            if (i == std_zero){
              cat(("-0-"))
            }else{
              cat(("-"))
            }
          }
        }
      }else{
        cat(" ")
      }
    }

    cat("\n")
    cat("TOST:           ")

    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_tost[1] && i <= std_fit_interval_tost[2]){
        if (i == std_fit_interval_tost[1]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("(-"))
          }else{
            if (lower_be_tost){
              cat(col_green("(-"))
            }else{
              cat(col_red("(-"))
            }
          }

        }else{
          if (i ==  std_fit_interval_tost[2]){
            if (i > std_be_interval[1] && i < std_be_interval[2]){
              cat(col_green("-)"))
            }else{
              if (upper_be_tost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }

          }else{
            if (i == std_fit_interval_center_tost){
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }

    cat("\n")

    if (x$method == "alpha-TOST"){
      cat("alpha-TOST:     ")
    }else{
      cat("cTOST:         ")
    }

    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_atost[1] && i <= std_fit_interval_atost[2]){
        if (i == std_fit_interval_atost[1]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("(-"))
          }else{
            if (lower_be_atost){
              cat(col_green("(-"))
            }else{
              cat(col_red("(-"))
            }
          }
        }else{
          if (i ==  std_fit_interval_atost[2]){
            if (i > std_be_interval[1] && i < std_be_interval[2]){
              cat(col_green("-)"))
            }else{
              if (upper_be_atost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }
          }else{
            if (i == std_fit_interval_center_atost){
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat("\n")
    cat("                 CI - low      ")
    cat("CI - high")
    cat("\n")

    cat("TOST:            ")
    cat(format(round(result_tost$ci[1], rn), nsmall = rn))
    cat("       ")
    cat(format(round(result_tost$ci[2], rn), nsmall = rn))
    cat("\n")

    if (x$method == "alpha-TOST"){
      cat("alpha-TOST:      ")
    }else{
      cat("cTOST:           ")
    }

    cat(format(round(x$ci[1], rn), nsmall = rn))
    cat("       ")
    cat(format(round(x$ci[2], rn), nsmall = rn))
    cat("\n")

    cat("\n")
    cat("Equiv. lim. = +/- ")
    cat(format(round(x$delta, rn), nsmall = rn))
    cat("\n")
  }
}



