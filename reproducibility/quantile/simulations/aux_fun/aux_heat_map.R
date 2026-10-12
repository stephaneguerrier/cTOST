get_color <- function(value, color_palette, color_breaks, n_colors) {
  if (is.na(value)) return("white")
  color_idx <- findInterval(value, color_breaks, rightmost.closed = TRUE)
  color_idx <- max(1, min(color_idx, n_colors))
  return(color_palette[color_idx])
}

custom_heatmap <- function(matrix,  
                           xlab = "X", ylab = "Y", main = "Heatmap",
                           show_values = FALSE, my_cex = 1, 
                           legend = TRUE, legend_title = "Value",
                           vline = NULL, hline = NULL, lambda=NULL, alpha=0.05,
                           custom_mfrow=NULL, max_val=NULL, min_val=NULL) {
  nrows <- nrow(matrix)
  ncols <- ncol(matrix)
  if (is.null(min_val)) min_val <- min(matrix, na.rm = TRUE)
  if (is.null(max_val)) max_val <- max(matrix, na.rm = TRUE)
  if (alpha>0){
    lb = 10
    blues = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(lb*2)[1:lb]
    la = lb*ceiling(max_val*10)/(alpha*10)
    reds = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(la*2)[-(1:la)]
    reds = reds[-(1:5)]
    la = length(reds)
    seq_below = seq(0, alpha, length.out = lb + 1)
    seq_above = seq(alpha, max_val+0.02, length.out = la + 1)[-1]
    fixed_breaks = c(seq_below, seq_above)
    colors = c(blues, reds)
    start_idx = which(fixed_breaks >= min_val)[1]
    end_idx = which(fixed_breaks <= max_val)
    end_idx = end_idx[length(end_idx)]
    color_palette = colors[start_idx:(end_idx-1)]
  } else {
    lb = 10
    blues = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(lb*2)[1:lb]
    la = 10
    reds = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(la*2)[-(1:la)]
    seq_below = seq(min(0, min(matrix)), alpha, length.out = lb)
    inds_below = !duplicated(seq_below, fromLast=T)
    seq_above = seq(alpha, max(0, max(matrix)), length.out = la)[-1]
    inds_above = !duplicated(seq_above)
    fixed_breaks = c(seq_below[inds_below], seq_above[inds_above])
    colors = c(blues[inds_below], reds[inds_above])
    start_idx = which(fixed_breaks >= min_val)[1]
    end_idx = which(fixed_breaks <= max_val)
    end_idx = end_idx[length(end_idx)]
    color_palette = colors[start_idx:(end_idx-1)]
  }
  n_colors <- length(color_palette)
  color_breaks <- seq(min_val, max_val, length.out = n_colors + 1)
  if (legend & is.null(custom_mfrow)) {
    layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  }
  par(mar = c(5.1, 4.6, 4.1, 0.5))
  ax_off = 1.8
  plot(0, 0, type = "n", 
       xlim = c(ax_off, ncols-ax_off), ylim = c(ax_off, nrows-ax_off),
       xlab = xlab, ylab = ylab, main = main,
       axes = FALSE, bty="n",
       cex.axis=my_cex+0.3, 
       cex.main=my_cex+0.3,
       cex.lab=my_cex+0.3)
  for (i in 1:nrows) {
    for (j in 1:ncols) {
      x_left <- j - 1
      x_right <- j
      y_bottom <- nrows - i
      y_top <- nrows - i + 1
      cell_color <- get_color(matrix[i, j], color_palette, color_breaks, n_colors)
      rect(x_left, y_bottom, x_right, y_top, 
           col = cell_color, border = "white", lwd = 0.5)
      if (show_values && !is.na(matrix[i, j])) {
        text(x_left + 0.5, y_bottom + 0.5, 
             round(matrix[i, j], 3), 
             cex = my_cex, col = "black")
      }
    }
  }
  inds_ax = seq(0, ncols, length.out = 6)
  inds_ax[1] = 1
  axis(1, at = seq(0.5, ncols - 0.5)[inds_ax], labels = colnames(matrix)[inds_ax], cex.axis= my_cex-0.1)
  axis(2, at = seq(0.5, nrows - 0.5)[inds_ax], labels = rev(rownames(matrix))[inds_ax], cex.axis=my_cex-0.1)
  axis(1, at = seq(0.5, ncols - 0.5), labels = NA)
  axis(2, at = seq(0.5, ncols - 0.5), labels = NA)
  line_col = "orange"
  line_lwd = 4
  if (!is.null(vline)) {
    for (v in vline) {
      if (!is.null(colnames(matrix))) {
        col_vals <- as.numeric(colnames(matrix))
        if (v >= min(col_vals) && v <= max(col_vals)) {
          col_idx <- which.min(abs(col_vals - v))
          x_pos <- col_idx - 0.5
          abline(v = x_pos, col = line_col, lwd = line_lwd, lty=1)
        }
      }
    }
  }
  if (!is.null(hline)) {
    for (h in hline) {
      if (!is.null(rownames(matrix))) {
        row_vals <- as.numeric(rownames(matrix))
        if (h >= min(row_vals) && h <= max(row_vals)) {
          row_idx <- which.min(abs(row_vals - h))
          y_pos <- nrows - row_idx + 0.5
          abline(h = y_pos, col = line_col, lwd = line_lwd, lty=1)
        }
      }
    }
  }
  box()
  col_idx <- which.min(abs(col_vals - lambda[1]))
  row_idx <- which.min(abs(row_vals - lambda[2]))
  y_pos <- nrows - row_idx + 0.5
  x_pos <- col_idx - 0.5
  points(x_pos, y_pos, pch = 19, col=line_col, cex = line_lwd)
  if (legend) {
    par(mar = c(5.1, 0.1, 4.1, 2.1))
    plot(0, 0, type = "n", xlim = c(-1, 1.7), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "", main = "")
    n_segments <- length(color_palette)
    col_thrs = 15
    legend_border = if (n_segments > col_thrs) NA else "white"
    segment_height <- 0.8 / n_segments
    legend_start <- 0.1
    for (i in 1:n_segments) {
      y_bottom <- legend_start + (i - 1) * segment_height
      y_top <- y_bottom + segment_height
      rect(0.2, y_bottom, 0.5, y_top,
           col = color_palette[i], border = legend_border, lwd = 0.5)
    }
    nleg = min(n_segments, col_thrs)+1
    legend_values <- seq(min_val, max_val, length.out = nleg)
    legend_positions <- seq(legend_start, legend_start + 0.8, length.out = nleg)
    for (i in 1:length(legend_values)) {
      lines(c(0.5, 0.55), c(legend_positions[i], legend_positions[i]))
      text(0.6, legend_positions[i], round(legend_values[i], 3), 
           cex = my_cex-0.45, adj = 0)
    }
    text(0.35, 0.95, legend_title, cex = my_cex-0.45, adj = 0.5, font = 2)
    if (is.null(custom_mfrow)) {
      layout(1)
    }
  }
}

