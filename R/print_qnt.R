library(cli)

#' @title Print Results of (Bio)Equivalence testing in Single Quantile Setting.
#'
#' @param x     A \code{qtost} object, which is the output of the function 'qtost'.
#' @param ticks Number of ticks to print the confidence interval in the console.
#' @param rn    Number of digits to consider when printing the results.
#' @param ...   Further arguments to be passed to or from methods.
#'
#' @return      Prints object.
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
    if (x$decision){
      cli_text(col_green("{symbol$tick} Accept quantile (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept quantile (bio)equivalence"))
    }
    if (x$method == "qTOST"){
      lower_be = x$ci[1] > x$delta[1]
      upper_be = x$ci[2] < x$delta[2]
      rg = range(c(x$ci, x$delta))
      rg_delta = rg[2] - rg[1]
      std_be_interval = round(ticks*(c(x$delta[1], x$delta[2]) - rg[1])/rg_delta) + 1
      std_zero = round(-ticks*rg[1]/rg_delta) + 1
      std_fit_interval = round(ticks*(x$ci - rg[1])/rg_delta) + 1
      std_fit_interval_center = round(ticks*(sum(x$ci)/2 - rg[1])/rg_delta) + 1
    }else{
      lower_be = x$ci[1] > x$delta[1]
      upper_be = x$ci[2] < x$delta[2]
      rg = range(c(x$ci, x$delta))
      rg_delta = rg[2] - rg[1]
      std_be_interval = round(ticks*(c(x$delta[1], x$delta[2]) - rg[1])/rg_delta) + 1
      std_zero = round(-ticks*rg[1]/rg_delta) + 1
      std_fit_interval = round(ticks*(x$ci - rg[1])/rg_delta) + 1
      std_fit_interval_center = round(ticks*(sum(x$ci)/2 - rg[1])/rg_delta) + 1
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
    cat("Equiv. lim. = (")
    cat(format(round(x$delta[1], rn), nsmall = rn))
    cat(" ; ")
    cat(format(round(x$delta[2], rn), nsmall = rn))
    cat(")")
    cat("\n")
    if (x$method == "alpha-qTOST"){
      cat("Corrected alpha = ")
      cat(format(round(x$corrected_alpha, rn), nsmall = rn))
      cat("\n")
    }
    cat("theta_hat = ")
    cat(format(round(x$theta_hat, rn), nsmall = rn))
    cat("; ")
    cat("Stand. dev. = ")
    cat(format(round(x$sigma_hat, rn), nsmall = rn))
    cat("\n")
  }
}

#' @title Print Results of (Bio)Equivalence testing in Two Quantiles Setting.
#'
#' @param x     A \code{qtost} object, which is the output of the function 'qtost'.
#' @param ticks Number of ticks to print the confidence interval in the console.
#' @param rn    Number of digits to consider when printing the results.
#' @param ...   Further arguments to be passed to or from methods.
#'
#' @return      Prints object.
#' @importFrom  cli cli_text col_green col_red
#'
#' @rdname print.m_qtost
#'
#' @export
print.m_qtost = function(x, ticks = 60, rn = 5,...){
  p = length(x$decision)
  if (all(x$decision)){
    cli_text(col_green("{symbol$tick} Accept quantile (bio)equivalence"))
  }else{
    cli_text(col_red("{symbol$cross} Can't accept quantile (bio)equivalence"))
  }
  rg = range(c(x$ci, x$delta))
  rg_delta = rg[2] - rg[1]
  std_zero= round(-ticks*rg[1] / rg_delta) + 1
  std_be_interval = std_fit_interval = matrix(NA,2,p)
  lower_be = upper_be = std_fit_interval_center = std_be_interval_center = rep(NA, p)
  for (i in 1:p) {
    std_be_interval[,i] = round(ticks*(x$delta[,i]-rg[1])/rg_delta) + 1
    std_fit_interval[,i] = round(ticks*(x$ci[,i] - rg[1])/rg_delta) + 1
    std_fit_interval_center[i] = round(ticks*(sum(x$ci[,i])/2 - rg[1])/rg_delta) + 1
    std_be_interval_center[i] = round(ticks*(sum(x$delta[,i])/2 - rg[1])/rg_delta) + 1
    lower_be[i] = x$ci[1,i] > x$delta[1,i]
    upper_be[i] = x$ci[2,i] < x$delta[2,i]
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
    cat(format(round(x$ci[1,i], rn), nsmall = rn))
    cat("; ")
    cat(format(round(x$ci[2,i], rn), nsmall = rn))
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
    cat(format(round(x$delta[1,i], rn), nsmall = rn))
    cat("; ")
    cat(format(round(x$delta[2,i], rn), nsmall = rn))
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
}

#' @title Comparison of a Corrective Procedure to the Results of the Quantile Two One-Sided Tests (qTOST) in Single Quantile Setting
#'
#' @description This function renders a comparison of the qTOST or the alpha-qTOST outputs obtained with the function `qtost`.
#'
#' @param x A \code{qtost} object, which is the output of one of the function: `qtost`.
#' @param ticks an integer indicating the number of segments that will be printed to represent the confidence intervals.
#' @param rn integer indicating the number of decimals places to be used (see function `round`) for the printed results.
#' @return Prints a comparison between the qTOST results (i.e., output of `qtost`) and the alpha-qTOST results.
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
compare_to_qtost = function(...) {
  args = list(...)
  if (!is.null(args)) {
    if (!is.null(args$x_data) && !is.null(args$y_data)) {
      x_data = args$x_data
      y_data = args$y_data
      x_bar = mean(x_data)
      sd_x = sd(x_data)
      n_x = length(x_data)
      y_bar = mean(y_data)
      sd_y = sd(y_data)
      n_y = length(y_data)
    } else {
      if (any(sapply(list(args$x_bar, args$sd_x, args$n_x, args$y_bar, args$sd_y, args$n_y), is.null))) {
        stop("Either provide raw data 'x' (x_data) and 'y' (y_data) or all summary statistics for x (x_bar, sd_x, n_x) and y (y_bar, sd_y, n_y).")
      } else {
        x_bar = args$x_bar
        sd_x = args$sd_x
        n_x = args$n_x
        y_bar = args$y_bar
        sd_y = args$sd_y
        n_y = args$n_y
      }
    }
    if (any(sapply(list(args$pi_x, args$delta_l, args$delta_u), is.null))) {
      stop("Please provide the quantile(s) of interest (pi_x), the lower limit(s) (delta_l) and the upper limit(s) (delta_u).\n")
    } else {
      pi_x = args$pi_x
      delta_l = args$delta_l
      delta_u = args$delta_u
      alpha = args$alpha
      B = args$B
      seed = args$seed
      tol = args$tol
      if (is.null(alpha)) alpha = 0.05
      if (is.null(B)) B = 10^5
      if (is.null(seed)) seed = 12345
      if (is.null(tol)) tol = .Machine$double.eps^0.5
      method = args$method
      if (is.null(method)) method = "alpha-qTOST"
      ticks = 30
      rn = 5
    }
  }
  if (exists("x_data") && exists("y_data")) {
    x = qtost(x_data = x_data, y_data = y_data,
              x_bar = x_bar, sd_x = sd_x, n_x = n_x,
              y_bar = y_bar, sd_y = sd_y, n_y = n_y,
              pi_x = pi_x, delta_l = delta_l, delta_u = delta_u,
              alpha = alpha, method = "alpha-qTOST")
    result_qtost = qtost(x_data = x_data, y_data = y_data,
                         x_bar = x_bar, sd_x = sd_x, n_x = n_x,
                         y_bar = y_bar, sd_y = sd_y, n_y = n_y,
                         pi_x = pi_x, delta_l = delta_l, delta_u = delta_u,
                         alpha = alpha, method = "unadjusted")
  } else {
    x = qtost(x_bar = x_bar, sd_x = sd_x, n_x = n_x,
              y_bar = y_bar, sd_y = sd_y, n_y = n_y,
              pi_x = pi_x, delta_l = delta_l, delta_u = delta_u,
              alpha = alpha, method = "alpha-qTOST")
    result_qtost = qtost(x_bar = x_bar, sd_x = sd_x, n_x = n_x,
                         y_bar = y_bar, sd_y = sd_y, n_y = n_y,
                         pi_x = pi_x, delta_l = delta_l, delta_u = delta_u,
                         alpha = alpha, method = "unadjusted")
  }
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
    lower_be_pitost = x$ci[1] > x$delta[1]
    upper_be_pitost = x$ci[2] < x$delta[2]
    lower_be_qtost = result_qtost$ci[1] > x$delta[1]
    upper_be_qtost = result_qtost$ci[2] < x$delta[2]
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
    rg = range(c(x$ci, result_qtost$ci, x$delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks * (x$delta - rg[1]) / rg_delta) + 1
    std_zero = round(-ticks * rg[1] / rg_delta) + 1
    std_fit_interval_pitost = round(ticks * (x$ci - rg[1]) / rg_delta) + 1
    std_fit_interval_center_pitost = round(ticks * (sum(x$ci) / 2 - rg[1]) / rg_delta) + 1
    std_fit_interval_qtost = round(ticks * (result_qtost$ci - rg[1]) / rg_delta) + 1
    # std_fit_interval_center_qtost = round(ticks * (sum(result_qtost$ci) / 2 - rg[1]) / rg_delta) + 1
    std_fit_interval_center_qtost = round(ticks * (sum(result_qtost$ci) / 2 - rg[1]) / rg_delta)-1

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
  cat(format(round(result_qtost$ci[1], rn), nsmall = rn))
  cat("       ")
  cat(format(round(result_qtost$ci[2], rn), nsmall = rn))
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
  cat(format(round(x$ci[1], rn), nsmall = rn))
  cat("       ")
  cat(format(round(x$ci[2], rn), nsmall = rn))
  cat("\n")
  cat("\n")
  cat("Equiv. lim. = dw/up ")
  cat(format(round(x$delta, rn), nsmall = rn))
  cat("\n")
}

