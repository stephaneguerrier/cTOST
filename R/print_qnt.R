#' @title Print Results of (Bio)Equivalence testing in Single Quantile Setting
#'
#' @param x     A \code{qtost} object, which is the output of the function 'qtost'.
#' @param ticks Number of ticks to print the confidence interval in the console.
#' @param rn    Number of digits to consider when printing the results.
#' @param ...   Further arguments to be passed to or from methods.
#'
#' @return      The object \code{x}, invisibly.
#' @importFrom  cli cli_text col_green col_red
#'
#' @rdname print.qtost
#'
#' @export
print.qtost = function(x, ticks = 30, rn = 5, ...){
  if(!is.list(x)){
    cat('"',x,'"',"\n")
  }else{
    if (!(x$method %in% c("qTOST", "alpha-qTOST"))){
      stop("This method is not compatible")
    }
    # qtost_core() stores the interval as a 1 x 2 matrix and the limits in eq_region
    ci = as.numeric(x$ci)
    be = as.numeric(x$eq_region[1, c(1, 3)])
    if (x$decision){
      cli_text(col_green("{symbol$tick} Accept quantile (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept quantile (bio)equivalence"))
    }
    if (x$method == "qTOST"){
      lower_be = ci[1] > be[1]
      upper_be = ci[2] < be[2]
      rg = range(c(ci, be))
      rg_delta = rg[2] - rg[1]
      std_be_interval = round(ticks*(c(be[1], be[2]) - rg[1])/rg_delta) + 1
      std_zero = round(-ticks*rg[1]/rg_delta) + 1
      std_fit_interval = round(ticks*(ci - rg[1])/rg_delta) + 1
      std_fit_interval_center = round(ticks*(sum(ci)/2 - rg[1])/rg_delta) + 1
    }else{
      lower_be = ci[1] > be[1]
      upper_be = ci[2] < be[2]
      rg = range(c(ci, be))
      rg_delta = rg[2] - rg[1]
      std_be_interval = round(ticks*(c(be[1], be[2]) - rg[1])/rg_delta) + 1
      std_zero = round(-ticks*rg[1]/rg_delta) + 1
      std_fit_interval = round(ticks*(ci - rg[1])/rg_delta) + 1
      std_fit_interval_center = round(ticks*(sum(ci)/2 - rg[1])/rg_delta) + 1
    }
    if (x$method == "qTOST"){
      cat("Equiv. Region:  ")
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
    if (x$method == "alpha-qTOST"){
      cat("Estim. Inter.: ")
    }else{
      cat("Estim. Inter.: ")
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
    cat(format(round(ci[1], rn), nsmall = rn))
    cat(" ; ")
    cat(format(round(ci[2], rn), nsmall = rn))
    cat(")\n\n")
    cat("Method: ")
    cat(x$method)
    cat("\n")
    cat("alpha = ")
    cat(x$alpha)
    cat("; ")
    cat("Equiv. lim. = (")
    cat(format(round(be[1], rn), nsmall = rn))
    cat(" ; ")
    cat(format(round(be[2], rn), nsmall = rn))
    cat(")")
    cat("\n")
    if (x$method == "alpha-qTOST"){
      cat("Corrected alpha = ")
      cat(format(round(x$corrected_alpha, rn), nsmall = rn))
      cat("\n")
    }
    cat("theta_hat = ")
    cat(format(round(x$theta, rn), nsmall = rn))
    cat("; ")
    cat("Stand. dev. = ")
    cat(format(round(x$sigma, rn), nsmall = rn))
    cat("\n")
  }
  invisible(x)
}

#' @title Print Results of (Bio)Equivalence testing in Two Quantiles Setting
#'
#' @param x     A \code{qtost} object, which is the output of the function 'qtost'.
#' @param ticks Number of ticks to print the confidence interval in the console.
#' @param rn    Number of digits to consider when printing the results.
#' @param ...   Further arguments to be passed to or from methods.
#'
#' @return      The object \code{x}, invisibly.
#' @importFrom  cli cli_text col_green col_red
#'
#' @rdname print.mqtost
#'
#' @export
print.mqtost = function(x, ticks = 60, rn = 5, ...){
  p = length(x$decision)
  # qtost_core() stores ci as p x 2 and eq_region as p x 3 (lower, pi_x, upper);
  # the display below works column-wise (2 x p)
  ci = t(x$ci)
  be = t(x$eq_region[, c(1, 3), drop = FALSE])
  if (all(x$decision)){
    cli_text(col_green("{symbol$tick} Accept quantile (bio)equivalence"))
  }else{
    cli_text(col_red("{symbol$cross} Can't accept quantile (bio)equivalence"))
  }
  rg = range(c(ci, be))
  rg_delta = rg[2] - rg[1]
  std_zero= round(-ticks*rg[1] / rg_delta) + 1
  std_be_interval = std_fit_interval = matrix(NA,2,p)
  lower_be = upper_be = std_fit_interval_center = std_be_interval_center = rep(NA, p)
  for (i in 1:p) {
    std_be_interval[,i] = round(ticks*(be[,i]-rg[1])/rg_delta) + 1
    std_fit_interval[,i] = round(ticks*(ci[,i] - rg[1])/rg_delta) + 1
    std_fit_interval_center[i] = round(ticks*(sum(ci[,i])/2 - rg[1])/rg_delta) + 1
    std_be_interval_center[i] = round(ticks*(sum(be[,i])/2 - rg[1])/rg_delta) + 1
    lower_be[i] = ci[1,i] > be[1,i]
    upper_be[i] = ci[2,i] < be[2,i]
  }
  adj_center = round(mean(std_fit_interval_center))
  adj_std_be_interval = adj_std_fit_interval = matrix(NA,2,p)
  for (i in 1:p){
    distance_be = adj_center-std_be_interval_center[i]
    adj_std_be_interval[,i] = std_be_interval[,i] + distance_be
    distance_fit = adj_center - std_fit_interval_center[i]
    adj_std_fit_interval[,i] = std_fit_interval[,i] + distance_fit
  }
  std_be_interval = adj_std_be_interval
  std_fit_interval = adj_std_fit_interval
  std_fit_interval_center = rep(adj_center,p)

  names_q = paste0("Equiv. Region for q", 1:p, ":")
  names_len_q = nchar(names_q)
  if (x$method == "qTOST"){
    names_s = paste0("qTOST for q", 1:p, ":")
    names_len_s = nchar(names_s)
    max_names_len = max(names_len_q,names_len_s)
    for (i in 1:p){
      name_q = names_q[i]
      cat(name_q)
      nchar_q = nchar(name_q)
      if (nchar_q <= max_names_len){
        cat(paste(rep(" ", max_names_len - nchar_q+1), collapse = "")) #used for cs
      }
      cat(" ")
      for (j in 1:(ticks + 1)) {
        if (j >= std_be_interval[1,i] && j <= std_be_interval[2,i]) {
          if (j == std_be_interval[1,i]){
            cat(("|-"))
          } else{
            if (j == std_be_interval[2,i]) {
              cat(("-|"))
            } else{
              if (j == std_zero) {
                cat(("-0-"))
              } else {
                cat(("-"))
              }
            }
          }
        } else {
          cat(" ")
        }
      }
      cat("\n")
      if(names_len_s[i] < max_names_len){
        cat(names_s[i])
        cat(paste(rep(" ", max_names_len-names_len_s[i]), collapse = ""))
        cat(" ")
      }else{
        cat(names_s[i])
        cat(" ")
      }
      for (j in 1:(ticks+1)){
        if (j >= std_fit_interval[1,i] && j <= std_fit_interval[2,i]){
          if (j == std_fit_interval[1,i]){
            if (j > std_be_interval[1,i] && j < std_be_interval[2,i]){
              cat(col_green("(-"))
            }else{
              if (lower_be[i]){
                cat(col_green("(-"))
              }else{
                cat(col_red("(-"))
              }
            }
          }else{
            if (j ==  std_fit_interval[2,i]){
              if (j > std_be_interval[1,i] && j < std_be_interval[2,i]){
                cat(col_green("-)"))
              }else{
                if (upper_be[i]){
                  cat(col_green("-)"))
                }else{
                  cat(col_red("-)"))
                }
              }
            }else{
              if (j == std_fit_interval_center[i]){
                if (j >= std_be_interval[1,i] && j <= std_be_interval[2,i]){
                  cat(col_green("-x-"))
                }else{
                  cat(col_red("-x-"))
                }
              }else{
                if (j >= std_be_interval[1,i] && j <= std_be_interval[2,i]){
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
  } else if (x$method == "alpha-qTOST"){
    names_c = paste0("alpha-qTOST for q", 1:p, ":")
    names_len_c = nchar(names_c)
    max_names_len = max(names_len_q,names_len_c)
    for (i in 1:p){
      name_q = names_q[i]
      cat(name_q)
      nchar_q = nchar(name_q)
      if (nchar_q < max_names_len){
        cat(paste(rep(" ", max_names_len - nchar_q), collapse = ""))
      }
      cat(" ")
      for (j in 1:(ticks + 1)) {
        if (j >= std_be_interval[1,i] && j <= std_be_interval[2,i]) {
          if (j == std_be_interval[1,i]){
            cat(("|-"))
          } else{
            if (j == std_be_interval[2,i]) {
              cat(("-|"))
            } else{
              if (j == std_zero) {
                cat(("-0-"))
              } else {
                cat(("-"))
              }
            }
          }
        } else {
          cat(" ")
        }
      }
      cat("\n")
      if(names_len_c[i] < max_names_len){
        cat(names_c[i])
        cat(paste(rep(" ", max_names_len-names_len_c[i]), collapse = ""))
        cat(" ")
      }else{
        cat(names_c[i])
        cat(" ")
      }
      for (j in 1:(ticks+1)){
        if (j >= std_fit_interval[1,i] && j <= std_fit_interval[2,i]){
          if (j == std_fit_interval[1,i]){
            if (j > std_be_interval[1,i] && j < std_be_interval[2,i]){
              cat(col_green("(-"))
            }else{
              if (lower_be[i]){
                cat(col_green("(-"))
              }else{
                cat(col_red("(-"))
              }
            }
          }else{
            if (j ==  std_fit_interval[2,i]){
              if (j > std_be_interval[1,i] && j < std_be_interval[2,i]){
                cat(col_green("-)"))
              }else{
                if (upper_be[i]){
                  cat(col_green("-)"))
                }else{
                  cat(col_red("-)"))
                }
              }
            }else{
              if (j == std_fit_interval_center[i]){
                if (j >= std_be_interval[1,i] && j <= std_be_interval[2,i]){
                  cat(col_green("-x-"))
                }else{
                  cat(col_red("-x-"))
                }
              }else{
                if (j >= std_be_interval[1,i] && j <= std_be_interval[2,i]){
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
  }
  cat("\n")
  cat("CIs:")
  cat("\n")
  for (i in 1:p){
    cat(paste0("q",i," = ("))
    cat(format(round(ci[1,i], rn), nsmall = rn))
    cat("; ")
    cat(format(round(ci[2,i], rn), nsmall = rn))
    cat(") ")
    if (x$decision[i]) {
      cat(col_green(cli::symbol$tick))
    } else {
      cat(col_red(cli::symbol$cross))
    }
    cat("\n")
  }
  cat("\n")
  cat("Equivalence limits: ")
  cat("\n")
  for (i in 1:p){
    cat(paste0("q",i," = ("))
    cat(format(round(be[1,i], rn), nsmall = rn))
    cat("; ")
    cat(format(round(be[2,i], rn), nsmall = rn))
    cat(")\n")
  }
  cat("\n")
  cat("Method: ")
  cat(x$method)
  cat("\n")
  cat("alpha = ")
  cat(x$alpha)
  if (x$method == "alpha-qTOST"){
    cat("; ")
    cat("Corrected alpha = ")
    cat(format(round(x$corrected_alpha, rn), nsmall = rn))
    cat("\n")
  }
  invisible(x)
}

#' @title Comparison of a Corrective Procedure to the Results of the Quantile Two One-Sided Tests (qTOST) in Single Quantile Setting
#'
#' @description This function renders a comparison of the qTOST or the alpha-qTOST outputs obtained with the function `qtost`.
#'
#' @param x A \code{qtost} object, which is the output of one of the function: `qtost`.
#' @param ticks an integer indicating the number of segments that will be printed to represent the confidence intervals.
#' @param rn integer indicating the number of decimals places to be used (see function `round`) for the printed results.
#' @return Prints a comparison between the qTOST results (i.e., output of `qtost`) and the alpha-qTOST results; returns \code{x} invisibly.
#'
#' @examples
#' # Using summary statistics from FDA label
#' x_bar_orig = 35.6
#' x_sd_orig = 16.7
#' n_x = 106
#' y_bar_orig = 41.6
#' y_sd_orig = 24.3
#' n_y = 14
#' x_bar = log(x_bar_orig^2 / sqrt(x_bar_orig^2 + x_sd_orig^2))
#' x_sd = sqrt(log(1 + (x_sd_orig^2 / x_bar_orig^2)))
#' y_bar = log(y_bar_orig^2 / sqrt(y_bar_orig^2 + y_sd_orig^2))
#' y_sd = sqrt(log(1 + (y_sd_orig^2 / y_bar_orig^2)))
#' x = list(mean=x_bar, sd=x_sd, n=n_x)
#' y = list(mean=y_bar, sd=y_sd, n=n_y)
#'
#' # alpha-qTOST
#' aqtost <- qtost(x, y, pi_x = 0.8, delta = 0.15, method = "alpha")
#' compare_to_qtost(aqtost)
#'
#' @importFrom cli cli_text col_green col_red
#'
#' @export
compare_to_qtost = function(x, ticks = 30, rn = 5) {
  if (!inherits(x, "qtost")) {
    stop("'x' must be a 'qtost' object, i.e. the output of 'qtost' for a single quantile.")
  }
  if (x$method != "alpha-qTOST") {
    stop("This method is not compatible")
  }
  # unadjusted qTOST on the same estimates, for comparison
  result_qtost = qtost_core(theta = x$theta, sigma = x$sigma, pi_x = x$pi_x,
                            delta_l = x$eq_region[, 1], delta_u = x$eq_region[, 3],
                            alpha = x$alpha)
  x_ci = as.numeric(x$ci)
  q_ci = as.numeric(result_qtost$ci)
  be = as.numeric(x$eq_region[1, c(1, 3)])
  name_q = "Equiv. Region: "
  name_len_q = nchar(name_q)
  name_s = "qTOST: "
  name_len_s = nchar(name_s)
  name_c = "alpha-qTOST: "
  name_len_c = nchar(name_c)
  max_names_len = max(name_len_q, name_len_s, name_len_c)

  if (!(x$method %in% c("qTOST", "alpha-qTOST"))) {
    stop("This method is not compatible")
  } else {
    lower_be_pitost = x_ci[1] > be[1]
    upper_be_pitost = x_ci[2] < be[2]
    lower_be_qtost = q_ci[1] > be[1]
    upper_be_qtost = q_ci[2] < be[2]
    if (x$method == "qTOST") {
      if (name_len_s < name_len_c) {
        # cat(name_s)
        # cat(paste(rep("", max_names_len - name_len_s), collapse = ""))
        # cat(" ")
      } else {
        # cat(name_s)
        # cat(" ")
      }
    } else {
      if (name_len_s < max_names_len) {
        # cat(name_s)
        # cat(paste(rep("", max_names_len - name_len_s), collapse = ""))
        # cat("")
      } else {
        # cat(name_s)
        # cat("")
      }
    }
  }
  if (result_qtost$decision) {
    if (name_len_s < name_len_c) {
      cli_text(paste0(
        name_s,
        paste(rep(" ", name_len_c - name_len_s), collapse = ""),
        col_green("{symbol$tick} Accept quantile (bio)equivalence")
      ))
    } else {
      cli_text(paste0(
        name_s,
        col_green("{symbol$tick} Accept quantile (bio)equivalence")
      ))
    }
  } else {
    if (name_len_s < name_len_c) {
      cli_text(paste0(
        name_s,
        paste(rep(" ", name_len_c - name_len_s), collapse = ""),
        col_red("{symbol$cross} Can't accept quantile (bio)equivalence")
      ))
    } else {
      cli_text(paste0(
        name_s,
        col_red("{symbol$cross} Can't accept quantile (bio)equivalence")
      ))
    }
  }

  if (x$method == "alpha-qTOST") {
    if (name_len_c < max_names_len) {
      # cat(name_c)
      # cat(paste(rep("", max_names_len - name_len_c), collapse = ""))
      if (x$decision) {
        # cat(" ")
      } else {
        # cat(name_c)
        # cat(" ")
      }
    } else {
      cat("only one adjustment method available:          ")
    }
    if (x$decision) {
      cli_text(paste0(name_c, col_green("{symbol$tick} Accept quantile (bio)equivalence")))
    } else {
      cli_text(paste0(name_c, col_red("{symbol$cross} Can't accept quantile (bio)equivalence")))
    }

    cat("\n")
    rg = range(c(x_ci, q_ci, be))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks * (be - rg[1]) / rg_delta) + 1
    std_zero = round(-ticks * rg[1] / rg_delta) + 1
    std_fit_interval_pitost = round(ticks * (x_ci - rg[1]) / rg_delta) + 1
    std_fit_interval_center_pitost = round(ticks * (sum(x_ci) / 2 - rg[1]) / rg_delta) + 1
    std_fit_interval_qtost = round(ticks * (q_ci - rg[1]) / rg_delta) + 1
    # std_fit_interval_center_qtost = round(ticks * (sum(q_ci) / 2 - rg[1]) / rg_delta) + 1
    std_fit_interval_center_qtost = round(ticks * (sum(q_ci) / 2 - rg[1]) / rg_delta)-1

    if (name_len_q < max_names_len) {
      cat(name_q)
      cat(paste(rep(" ", max_names_len - name_len_q), collapse = ""))
      cat("  ")
    } else {
      cat(name_q)
      cat("  ")
    }
  }
  for (i in 1:(ticks + 1)) {
    if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
      if (i == std_be_interval[1]) {
        cat(("|-"))
      } else {
        if (i == std_be_interval[2]) {
          cat(("-|"))
        } else {
          if (i == std_zero) {
            cat(("-0-"))
          } else {
            cat(("-"))
          }
        }
      }
    } else {
      cat(" ")
    }
  }
  cat("\n")
  if (name_len_s < max_names_len) {
    cat(name_s)
    cat(paste(rep(" ", max_names_len - name_len_s), collapse = ""))
    cat(" ")
  } else {
    cat(name_s)
    cat(" ")
  }
  for (i in 1:(ticks + 1)) {
    if (i >= std_fit_interval_qtost[1] && i <= std_fit_interval_qtost[2]) {
      if (i == std_fit_interval_qtost[1]) {
        if (i > std_be_interval[1] && i < std_be_interval[2]) {
          cat(col_green("(-"))
        } else {
          if (lower_be_qtost) {
            cat(col_green("(-"))
          } else {
            cat(col_red("(-"))
          }
        }
      } else {
        if (i == std_fit_interval_qtost[2]) {
          if (i > std_be_interval[1] && i < std_be_interval[2]) {
            cat(col_green("-)"))
          } else {
            if (upper_be_qtost) {
              cat(col_green("-)"))
            } else {
              cat(col_red("-)"))
            }
          }
        } else {
          if (i == std_fit_interval_center_qtost) {
            if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
              cat(col_green("-x-"))
            } else {
              cat(col_red("-x-"))
            }
          } else {
            if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
              cat(col_green("-"))
            } else {
              cat(col_red("-"))
            }
          }
        }
      }
    } else {
      cat(" ")
    }
  }
  cat("\n")
  if (x$method == "alpha-qTOST") {
    if (name_len_c < max_names_len) {
      cat(name_c)
      cat(paste(rep(" ", max_names_len - name_len_c), collapse = ""))
      cat(" ")
    } else {
      cat(name_c)
      cat(" ")
    }
  } else {
    cat("only one adjustment method available:          ")
  }
  for (i in 1:(ticks + 1)) {
    if (i >= std_fit_interval_pitost[1] && i <= std_fit_interval_pitost[2]) {
      if (i == std_fit_interval_pitost[1]) {
        if (i > std_be_interval[1] && i < std_be_interval[2]) {
          cat(col_green("(-"))
        } else {
          if (lower_be_pitost) {
            cat(col_green("(-"))
          } else {
            cat(col_red("(-"))
          }
        }
      } else {
        if (i == std_fit_interval_pitost[2]) {
          if (i > std_be_interval[1] && i < std_be_interval[2]) {
            cat(col_green("-)"))
          } else {
            if (upper_be_pitost) {
              cat(col_green("-)"))
            } else {
              cat(col_red("-)"))
            }
          }
        } else {
          if (i == std_fit_interval_center_pitost) {
            if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
              cat(col_green("-x-"))
            } else {
              cat(col_red("-x-"))
            }
          } else {
            if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
              cat(col_green("-"))
            } else {
              cat(col_red("-"))
            }
          }
        }
      }
    } else {
      cat(" ")
    }
  }
  name_low = "               CI - low      "
  name_high = "CI - high"
  cat("\n")
  cat("\n")
  cat(name_low)
  cat(name_high)
  cat("\n")
  if (name_len_s < nchar(name_low)) {
    cat(name_s)
    cat(paste(rep(" ", 7), collapse = ""))
    cat(" ")
  } else {
    cat(name_s)
    cat(" ")
  }
  cat(format(round(q_ci[1], rn), nsmall = rn))
  cat("       ")
  cat(format(round(q_ci[2], rn), nsmall = rn))
  cat("\n")
  if (x$method == "alpha-qTOST") {
    if (name_len_c < nchar(name_low)) {
      cat(name_c)
      cat(paste(rep(" ", 1), collapse = ""))
      cat(" ")
    } else {
      cat(name_c)
      cat(" ")
    }
  } else {
    cat("only one adjustment method available:          ")
  }
  cat(format(round(x_ci[1], rn), nsmall = rn))
  cat("       ")
  cat(format(round(x_ci[2], rn), nsmall = rn))
  cat("\n")
  cat("\n")
  cat("Equiv. lim. = dw/up ")
  cat(format(round(be, rn), nsmall = rn))
  cat("\n")
  invisible(x)
}

