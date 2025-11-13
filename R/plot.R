#' Plot confidence intervals for one or more `tost` or `mtost` objects
#'
#' @description
#' Displays and compares confidence intervals and rejections regions from one or more
#' equivalence test objects (`tost` or `mtost`). The plot displays point estimates
#' and confidence intervals against a shaded equivalence region, providing a clear
#' visual summary of the test results.
#'
#' @param ... One or more `tost` or `mtost` objects to be plotted, passed as
#'   separate arguments.
#' @param plot_params A named list of parameters to customize the plot's appearance.
#'   See the "Customization" section for details on key options.
#' @param legend_params A named list of parameters to control the legend. The legend
#'   is only displayed for multivariate plots (where K > 1). Set to `NULL` to
#'   disable. See the "Customization" section for details.
#'
#' @details
#' ## Customization via Lists
#' The function's appearance is controlled by passing named lists to the
#' `plot_params` and `legend_params` arguments.
#'
#' ### Key `plot_params` options:
#' - `main`, `xlab`: Title and x-axis label for the plot.
#' - `lwd`, `pch`, `cex`: Line width, point character, and base size for CIs.
#' - `col`: A vector of colors for the different methods/objects.
#' - `var_names`: Character vector to override the y-axis labels.
#' - `c0_lab`: A vector of two expressions for the equivalence boundary labels.
#' - `eq_region_fill`: Color for the shaded equivalence region.
#' - `cex.axis`, `cex.main`, etc.: Size controls for specific text elements.
#' - `manage_par`: A logical value. If `TRUE` (the default), the function manages
#'   its own graphical parameters (`par`) and resets them upon exiting. Set to
#'   `FALSE` when arranging multiple plots in a grid (e.g., with `par(mfrow=...))`.
#' - `add_decision`: A logical value. If `TRUE` (the default), ticks showing the
#'    decisions for equivalence assessment are added next to confidence intervals.
#'
#' ### Key `legend_params` options:
#' - `x`, `title`, `cex`: Standard `legend()` arguments for position, title, and size.
#' - `equal_spacing`: If `TRUE`, ensures constant spacing between legend items.
#' - `spacing_vec`: A numeric vector to provide exact custom spacing after each item,
#'   overriding `equal_spacing`.
#'
#' @return NULL. The function generates a plot.
#'
#' @export
#'
#' @examples
#' # Univariate assessment
#' data(skin)
#' theta_hat = diff(apply(skin,2,mean))
#' n = nrow(skin)
#' nu = n-1
#' sig_hat = var(apply(skin,1,diff))/n
#' alpha0 = 0.05
#' c0 = log(1.25)
#' # Univariate comparison with default inputs
#' stost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = alpha0, delta = c0, method = "unadjusted")
#' atost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = alpha0, delta = c0, method = "alpha")
#' otost = ctost(theta = theta_hat, sigma = sig_hat,
#'               nu = nu, alpha = alpha0, delta = c0, method = "optimal")
#' plot.tost(stost, atost, otost)
#' # Multivariate assessment
#' data(skin_mvt)
#' n = nrow(skin_mvt)
#' nu = n-1
#' theta_hat = apply(skin_mvt, 2, mean)
#' Sigma_hat = cov(skin_mvt)/n
#' alpha0 = 0.05
#' c0 = log(1.25)
#' # Multivariate assessment of a single method with default inputs
#' (mvt_stost = ctost(theta = theta_hat, sigma = Sigma_hat, nu = nu,
#'                    alpha = alpha0, delta = c0, method = "unadjusted"))
#' plot.tost(mvt_stost)
#' # Multivariate comparison with default inputs
#' (mvt_atost = ctost(theta = theta_hat, sigma = Sigma_hat, nu = nu,
#'                    alpha = alpha0, delta = c0, method = "alpha", B=1e3))
#' (mvt_ctost = ctost(theta = theta_hat, sigma = Sigma_hat, nu = nu,
#'                    alpha = alpha0, delta = c0, method = "optimal"))
#' plot.tost(mvt_stost, mvt_atost, mvt_ctost)
#' # Multivariate comparison with custom inputs
#' plot.tost(
#'   mvt_stost, mvt_atost, mvt_ctost,
#'   plot_params = list(
#'     main = "Bioequivalence Assessment",
#'     xlab = "Effect size",
#'     var_names = c(expression(italic("Stratum corneum")),
#'                   expression(italic("Viable epidermis")),
#'                   expression(italic("Upper dermis")),
#'                   expression(italic("Lower dermis"))),
#'     pch = 15,
#'     lwd = 3,
#'     mar_adj = c(0, 8, 4, 0),
#'     cex.axis=1.5,
#'     cex.main=2,
#'     line.main=4,
#'     line.ylab=2,
#'     eq_region_fill = scales::alpha("grey60", 0.15),
#'     eq_region_lines = "grey60",
#'     add_decision = F
#'   ),
#'   legend_params = list(
#'     legend = c("TOST", bquote(alpha*"-TOST"), "cTOST"),
#'     x = -1.15,
#'     y = 10,
#'     inset = -0.7,
#'     title = "Method:",
#'     bty = "o",
#'     horiz = F
#'   )
#' )
plot.tost = function(..., plot_params = list(), legend_params = list()) {
  x = list(...)
  if ((!all(sapply(x, class) == "tost") && !all(sapply(x, class) == "mtost") &&
       !all(sapply(x, class) == "qtost") && !all(sapply(x, class) == "mqtost"))){
    stop("Input 'x' must be a list of the same 'tost', 'qtost', 'mtost'  or 'mqtost' objects.")
  }
  qnt_plot = if (any(sapply(x, "[[", "method") %in% c("qTOST", "alpha-qTOST"))) T else F
  M = length(x)
  sapply(x, class)
  K = length(x[[1]]$decision)
  # Set default plotting parameters
  if (qnt_plot) {
    if (K>1) {
      text_expressions = sprintf("pi[x] == %.2f", x[[1]]$pi_x)
      var_names_def = sapply(text_expressions, function(x) parse(text = x))
    } else {
      var_names_def = sapply(x, "[[", "method")
    }
    # mar_adj_def = if (K>1 && M>1) c(0, 3, 3, 1) else c(0, 3, 0, 1)
    all_cis = do.call(rbind, lapply(x, "[[", "ci")) - x[[1]]$pi_x
  } else {
    var_names_def = if (K>1) names(x[[1]]$decision) else sapply(x, "[[", "method")
    # mar_adj_def = if (K>1 && M>1) c(0, 6, 3, 1) else c(0, 6, 0, 1)
    all_cis = do.call(rbind, lapply(x, "[[", "ci"))
  }
  # l_spac = max(sapply(var_names_def, strwidth, units = "user", cex = par("cex"))) * 50
  l_spac = max(strwidth(var_names_def, "inches")) / par("cin")[1] + 1.5
  p = list(
    main = "",
    xlab = if (qnt_plot) expression(pi[y]-pi[x]) else expression(theta),
    lwd = 5,
    pch = 19,
    col = gg_color_hue(M),
    ylim_offset = 1.5,
    cap_width = 0.05,
    dodge_width = 1,
    mar_adj = if (K>1 && M>1) c(0, l_spac, 3, 1) else c(0, 3, 0, 1),
    var_names = var_names_def,
    c0_lab = c(expression(-c[0]), expression(c[0])),
    eq_region_fill = scales::alpha("#FF9900", 0.1),
    eq_region_lines = "#FF9900",
    add_decision = T,
    manage_par = TRUE,
    cex = 1.5,
    cex.axis = 1.5,
    cex.main = 1.6,
    cex.xlab = 1.7,
    cex.ylab = 1.6,
    line.main = 2.5,
    line.xlab = 3.5,
    line.ylab = 1,
    line.c0 = 0.1
  )
  p$var_names[p$var_names=="alpha-TOST"] = expression(paste(alpha, "-TOST"))
  p$var_names[p$var_names=="delta-TOST"] = expression(paste(delta, "-TOST"))
  p$var_names[p$var_names=="alpha-qTOST"] = expression(paste(alpha, "-qTOST"))
  p = utils::modifyList(p, plot_params)
  # Calculate plot coordinates and limits
  if (K>1){
    y_coords = if (is.null(p$ylim_offset)) K:1 else seq(K + p$ylim_offset, 1 - p$ylim_offset, len = K)
  } else {
    y_coords = if (is.null(p$ylim_offset)) M:1 else seq(M + p$ylim_offset/3, 1 - p$ylim_offset/3, len = M)
  }
  c0 = x[[1]]$delta
  xlim_range = abs(range(all_cis, c0 * 2.5))
  xlim = c(-max(xlim_range), max(xlim_range))
  # xlim = c(-c0*2.5, c0*2.5)
  ylim = c(min(y_coords) - p$dodge_width, max(y_coords) + p$dodge_width)
  # Plotting
  if (p$manage_par) {
    old_par = par(no.readonly = TRUE)
    on.exit(par(old_par))
  }
  par(mar = c(5.1 + p$mar_adj[1], 4.1 + p$mar_adj[2], 4.1 + p$mar_adj[3], 2.1 + p$mar_adj[4]))
  plot(NA, NA, xlim = xlim, ylim = ylim, axes = FALSE, main = "", xlab = "", ylab = "")
  rect(-c0, par("usr")[3], c0, par("usr")[4], border = NA, col = p$eq_region_fill)
  abline(v = 0, lty = 3, col = "black", lwd = 2)
  abline(v = c(-c0, c0), lty = 3, col = p$eq_region_lines, lwd = 1)
  # abline(h = y_coords, lty = 1, col = "grey80")
  plot_coords = par("usr")
  left_edge = plot_coords[1]
  right_edge = plot_coords[2]
  symbol_x_pos = xlim[2] + (right_edge - xlim[2]) * 0.5
  line_end_x = symbol_x_pos - (right_edge - xlim[2]) * 0.2
  segments(x0 = left_edge, y0 = y_coords, x1 = line_end_x, y1 = y_coords,
           col = "grey80", lty = 1)
  if (K>1) j_all = 1:K else j_all = 1:M
  for (j in j_all) {
    mtext(p$var_names[j], side = 2, at = y_coords[j], las = 1, cex = p$cex.ylab, adj = 1.05, line=p$line.ylab)
  }
  # Add CI's
  offsets = if (M>1 && K>1) seq(p$dodge_width/2, -p$dodge_width/2, length.out = M) else 0
  for (i in 1:M) {
    obj = x[[i]]
    for (j in 1:K) {
      if (K>1) {
        y_pos = y_coords[j] + offsets[i]
      } else {
        y_pos = y_coords[i] + offsets[j]
      }
      if (qnt_plot) {
        mean_point = obj$pi_y_hat[j] - obj$pi_x[j]
        CIs = matrix(obj$ci, nrow=K, ncol=2) - obj$pi_x
      } else {
        mean_point = obj$theta[j]
        CIs = matrix(obj$ci, nrow=K, ncol=2)
      }
      arrows(CIs[j, 1], y_pos, CIs[j, 2],
             y_pos, angle = 90, code = 3,
             length = p$cap_width, col = p$col[i], lwd = p$lwd)
      points(mean_point, y_pos, col = p$col[i], pch = p$pch, cex = p$cex * 1.25)
      if (p$add_decision){
        is_accepted = obj$decision[j]
        symbol_char = if (obj$decision[j]) "✓" else "✗"
        symbol_col = if (obj$decision[j]) "darkgreen" else "red"
        par(xpd = T)
        text(x = symbol_x_pos*1.035, y = y_pos,
             labels = symbol_char,
             col = symbol_col,
             cex = p$cex * 1.4,
             font = 2,
             adj = c(0, 0.2))
        par(xpd = F)
      }
    }
  }
  # Add axes and titles
  axis(1, cex.axis = p$cex.axis)
  if (!is.null(p$main) && p$main != "") mtext(p$main, at = 0, side = 3, line = p$line.main, cex = p$cex.main)
  if (!is.null(p$xlab) && p$xlab != "") mtext(p$xlab, at = 0, side = 1, line = p$line.xlab, cex = p$cex.xlab)
  if (!is.null(p$c0_lab)) mtext(p$c0_lab, at = c(-c0, c0), side = 3, line = p$line.c0, cex = p$cex.axis)
  # Add legend
  if (!is.null(legend_params) && K>1 && M>1) {
    par(xpd = TRUE)
    leg_defaults = list(
      x = mean(par("usr")[1:2]), # "top",
      y = par("usr")[4],
      legend = sapply(x, "[[", "method"),
      title = "",
      col = p$col,
      pch = p$pch,
      lwd = p$lwd,
      bty = "n",
      cex = p$cex,
      inset = -0.35,
      x.intersp = 1.5,
      equal_spacing = T,
      symbol_space = 1.8,
      spacing_vec = NULL,
      bar_width_ratio = 0.7,
      cap_length = 0.04,
      horiz = T,
      xjust = 0.5,
      yjust = -0.7
    )
    leg_defaults$legend[leg_defaults$legend=="alpha-TOST"] = expression(paste(alpha, "-TOST"))
    leg_defaults$legend[leg_defaults$legend=="alpha-qTOST"] = expression(paste(alpha, "-qTOST"))
    leg_args = utils::modifyList(leg_defaults, legend_params)
    if (!is.null(leg_args$spacing_vec)) {
      # If the user provided a vector, use it directly
      if (length(leg_args$spacing_vec) != M) {
        warning(paste("`spacing_vec` length must be", M,
                      "- using automatic spacing instead."))
      } else {
        # Use the user's vector and disable the automatic calculation
        leg_args$x.intersp = leg_args$spacing_vec
        leg_args$equal_spacing = FALSE
      }
    }
    add_ci_legend(leg_args)
    par(xpd = FALSE)
  }
}

#' Generate ggplot2-like Colors
#'
#' @description Creates a vector of colors matching the default ggplot2 palette
#' by spacing hues evenly around the HCL color wheel.
#'
#' @param n The number of colors to generate.
#'
#' @return A character vector of hex color codes.
#'
#' @keywords internal
#'
gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Null-Coalescing Operator
#'
#' @description Returns the left-hand side if it is not NULL, otherwise returns
#' the right-hand side. This provides a convenient shorthand for setting default values.
#'
#' @param a The value to check for NULL.
#' @param b The default value to return if 'a' is NULL.
#'
#' @return The first non-NULL value from left to right.
#'
#' @keywords internal
#'
`%||%` = function(a, b) {
  if (is.null(a)) b else a
}

#' Add a Custom Legend with Confidence Interval Symbols
#'
#' @description
#' This is a helper function designed to be called by a main plotting function. It
#' draws a plot legend where the default keys (points, lines) are replaced with
#' custom symbols that mimic confidence intervals (a point with capped error bars).
#' It also dynamically calculate horizontal spacing for legend items to ensure
#' a clean, constant layout, even when labels have variable lengths.
#'
#' @param args A single list containing arguments. This list should include
#'   standard parameters to be passed to the base `legend()` function (e.g., `x`,
#'   `legend`, `title`, `col`, `pch`, `lwd`, `bty`, `cex`) as well as the custom
#'   parameters detailed below.
#'
#' @details
#' ## How It Works
#' The function operates in three main steps:
#' 1.  It calls `legend(..., plot = FALSE)` to calculate the coordinates for the
#'     legend box and text without drawing anything.
#' 2.  It calls `legend()` a second time to draw *only* the text labels.
#' 3.  It uses the coordinates from step 1 to manually draw the custom CI
#'     symbols in the space where the normal legend keys would be.
#'
#' ## Custom Spacing and Symbols
#' The behavior is controlled by several custom list items within the `args` parameter:
#' - **`equal_spacing`**: Logical. If `TRUE` (the default) and the legend is
#'   horizontal (`horiz = TRUE`), the function dynamically calculates a vector of
#'   `x.intersp` values. This ensures each legend item occupies the same
#'   total width, creating constant visual spacing between them.
#' - **`symbol_space`**: Numeric. A factor controlling the total horizontal
#'   space reserved for each custom symbol to the left of the text.
#' - **`bar_width_ratio`**: Numeric. The width of the CI error bar in the symbol,
#'   expressed as a ratio of the total `symbol_space`.
#' - **`cap_length`**: Numeric. The length of the end caps on the CI symbol,
#'   which is passed to the `arrows()` function.
#'
#' @return NULL. The function updates a plot.
#'
#' @keywords internal
#'
add_ci_legend = function(args) {
  # Define custom arguments
  custom_arg_names = c("symbol_space", "bar_width_ratio", "cap_length", "equal_spacing", "spacing_vec")
  # Clean list to pass to the base legend() function
  base_legend_args = args[!names(args) %in% custom_arg_names]
  # Fix spacing
  # if (!identical(args$equal_spacing, FALSE)) {
  #   pt.cex = args$cex %||% par("cex")
  #   max_width = max(strwidth(args$legend, units = "user", cex = pt.cex))
  #   base_legend_args$text.width = max_width
  # }
  pt.cex = args$cex %||% par("cex")
  if (!identical(args$equal_spacing, FALSE) && (args$horiz %||% FALSE)) {
    # Get the width of a standard space character
    space_char_width = strwidth("M", units = "user", cex = pt.cex)
    # Get the width of each legend label
    label_widths = strwidth(args$legend, units = "user", cex = pt.cex)
    max_label_width = max(label_widths)
    # Define a base interspacing (the minimum gap after the longest item)
    base_intersp = args$x.intersp %||% 1
    # Calculate the required extra space for each label to match the max width
    padding_needed = max_label_width - label_widths
    # Convert this padding into an x.intersp value and add the base gap
    # This creates a vector where shorter labels get a huge interspace
    dynamic_intersp = base_intersp + (padding_needed / space_char_width)
    # Update the arguments that will be passed to the legend() function
    base_legend_args$x.intersp = dynamic_intersp
  }
  # Step 1: Get coordinates for the legend elements without drawing
  coord_args = base_legend_args
  coord_args$plot = FALSE
  leg_coords = do.call(legend, coord_args)
  # Step 2: Draw the text part of the legend only
  text_args = base_legend_args
  text_args$lty = 0
  text_args$pch = NA
  do.call(legend, text_args)
  # Step 3: Draw the custom CI symbols
  symbol_space = args$symbol_space %||% 1.8
  bar_ratio = args$bar_width_ratio %||% 0.4
  cap_len = args$cap_length %||% 0.05
  # Calculate positions for drawing the symbols
  # browser()
  offset = symbol_space * strwidth("M", cex = pt.cex, font = par("font"))
  y_coords = leg_coords$text$y # par("usr")[4]+2  #
  x_coords = leg_coords$text$x - offset
  bar_width = offset * bar_ratio
  # Get graphical parameters, recycling them if necessary
  cols = args$col %||% par("col")
  pchs = args$pch %||% par("pch")
  lwds = args$lwd %||% par("lwd")
  # Loop through each legend item and draw its custom symbol
  for (i in seq_along(y_coords)) {
    current_col = cols[((i - 1) %% length(cols)) + 1]
    current_pch = pchs[((i - 1) %% length(pchs)) + 1]
    current_lwd = lwds[((i - 1) %% length(lwds)) + 1]

    arrows(x_coords[i] - bar_width, y_coords[i], x_coords[i] + bar_width, y_coords[i],
           angle = 90, code = 3, length = cap_len, col = current_col, lwd = current_lwd)
    points(x_coords[i], y_coords[i], pch = current_pch, col = current_col, cex = pt.cex * 1.2)
  }
}
