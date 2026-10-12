rm(list = ls())

if (!require("viridisLite")) install.packages("viridisLite")
library(viridisLite)

tikzPlot = 0
if (tikzPlot == 1) {
  if (!require("tikzDevice")) install.packages("tikzDevice")
  library(tikzDevice)
}

n_simu = 50000
n_array = 5000
loadDir = "simulations/HPC/sim_univ/sim_HPC_out/sim_setting_2"
name_sim = paste0("/results_sim_setting_2_n_simu_", n_simu, "_n_array_", n_array)
load(file=paste0(loadDir, name_sim, ".rdata"))
sol_qtost_avg = sol_TOST_avg
sol_alpha_qtost_avg = sol_alphaTOST_avg

sol = sol_alpha_qtost_avg - sol_qtost_avg

m = 500
q1 = seq(from = 0, to = 27, length.out = m)

r2 = function(k) viridisLite::inferno(k, direction = -1)

col2 = r2(m)

cex_tit = 1.2
cex_lab = 1.2 - 0.225
cex_lab_n = cex_lab
cell_scale = 1.15

n = dim(sol)[5]
m = prod(dim(sol)[1:4])
m = m / 2

off = 50
lab_y = m + 4
label_offset = 1
col_step = 3.5
pi_gap = 1
pick_l_range = 1:length(dimnames(sol_qtost_avg)[[2]])
pick_gamma_range = 1:length(dimnames(sol_qtost_avg)[[4]])
pick_n_range = 1:length(dimnames(sol_qtost_avg)[[1]])
ny_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[1]]))
l_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[2]]))
gamma_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[4]]))
assignments = strsplit(dimnames(sol_qtost_avg)[[3]], ";")
pi_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 1)))
delta_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 2)))

y_leg = -5 - pi_gap
xlo = -12.3
xhi = 51
ylo = y_leg - 5
yhi = lab_y + 0.7
fig_w = 6 * cell_scale
fig_h = fig_w * (yhi - ylo) / (xhi - xlo)

draw_panel = function(pi_range) {
  gy = function(gb) -pi_gap * gb

  i = 1
  for (pick_pi_x in pi_range) {
    g = gy(pick_pi_x - pi_range[1])
    for (pick_gamma in pick_gamma_range) {
      for (pick_l in pick_l_range) {
        for (pick_n in pick_n_range) {
          for (j in 1:n) {
            x1 = sol_alpha_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, j] - sol_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, j]
            rect(j - 0.5, n + 1 - i - 0.5 + g, j + 0.5, n + 1 - i + 0.5 + g, col = col2[which.min(abs(q1 - x1))], border = "white")
          }
          i = i + 1
        }
      }
    }
  }

  text(0, lab_y, "$n_y$", cex = cex_tit)
  for (i in 0:(m-1)) {
    g = gy(i %/% 24)
    lines(c(0.5, 0.1), -c(i,i) + off + g)
    text(0.7, -i + off + g, ny_vec[(i %% 4) + 1], cex = cex_lab_n, pos = 2)
  }

  arrow1x = -col_step
  text(arrow1x, lab_y, "$l$", cex = cex_tit)
  for (i in 1:(6*2)) {
    j = i - 1
    g = gy((i - 1) %/% 6)
    arrows(arrow1x, 50 - 4*(i-1) + g, arrow1x, 47 - 4*(i-1) + g, code = 3, length = 0.03)
    if (l_vec[(j %% 3) + 1] > 1) {
      tt = paste0("$1/", l_vec[(j %% 3) + 1], "$")
    } else {
      tt = l_vec[(j %% 3) + 1]
    }
    text(arrow1x - label_offset, mean(c(50 - 4*(i-1), 47 - 4*(i-1))) + g, tt, cex = cex_lab, srt = 90)
  }

  arrow1x = -2 * col_step
  text(arrow1x, lab_y, "$\\gamma$", cex = cex_tit)
  for (i in 1:4) {
    j = i - 1
    g = gy((i - 1) %/% 2)
    arrows(arrow1x, 50 - 3*4*(i-1) + g, arrow1x, 39 - 3*4*(i-1) + g, code = 3, length = 0.03)
    text(arrow1x - label_offset, mean(c(50 - 3*4*(i-1), 39 - 3*4*(i-1))) + g, gamma_vec[(j %% 2) + 1], cex = cex_lab, srt = 90)
  }

  arrow1x = -3 * col_step
  text(arrow1x, lab_y, "$\\pi_x$", cex = cex_tit)
  for (i in 1:2) {
    g = gy(i - 1)
    arrows(arrow1x, 50 - 2*3*4*(i-1) + g, arrow1x, 27 - 2*3*4*(i-1) + g, code = 3, length = 0.03)
    tt = paste0(pi_vec[pi_range[i]], " ($c=", delta_vec[pi_range[i]], "$)")
    text(arrow1x - label_offset, mean(c(50 - 2*3*4*(i-1), 27 - 2*3*4*(i-1))) + g, tt, cex = cex_lab, srt = 90)
  }

  x_xax = c(0.5, 50.5)
  seq_vals = seq(from = min(x_xax) + 0.5, to = max(x_xax) - 0.5, l = 50)
  vals = round(seq(from = -1.2*10, to = 1.2*10, l = 50)) / 10
  y_bot = 2.5 + gy(1)
  lines(x_xax, rep(y_bot, 2))
  for (i in 1:50) {
    lines(rep(seq_vals[i], 2), y_bot + c(0, -0.5))
    if (i == 1 || !(i %% 5)) {
      text(seq_vals[i], y_bot - 1.5, vals[i], cex = cex_lab)
      if (vals[i] == 0) text(seq_vals[i], y_bot - 3.5, "$\\kappa$", cex = cex_tit)
      if (vals[i] == -1 || vals[i] == 1) {
        segments(seq_vals[i], 26.5 + gy(0), seq_vals[i], 50.5 + gy(0), lty = 2)
        segments(seq_vals[i],  2.5 + gy(1), seq_vals[i], 26.5 + gy(1), lty = 2)
      }
    }
  }
}

draw_legend = function() {
  l = 20 - 2
  cols = col2[vapply(seq(0, 18, length.out = l),
                    function(value) which.min(abs(q1 - value)), integer(1))]
  k = 100
  mm = length(cols)
  myseq = (n + 1)/2 - (k/2.8)/2 + seq(from = 0, to = k/2.8, length.out = mm)
  dseq = myseq[2] - myseq[1]
  sq_gap = 0.15
  rect(myseq - dseq/2 + sq_gap/2, y_leg - 0.9, myseq + dseq/2 - sq_gap/2, y_leg + 0.9, col = cols, border = NA)
  text((n + 1)/2, y_leg - 5, "Rejection probability gain of $\\alpha$-qTOST over qTOST (in \\%)", cex = cex_lab)
  line_range = range(myseq) + c(-1, 1)
  lines(line_range, c(y_leg - 0.9, y_leg - 0.9))
  jj = 13
  vlab = seq(from = 0, to = 18, length.out = jj)
  xlab = seq(from = line_range[1], to = line_range[2], length.out = jj)
  for (i in 1:jj) {
    if (i %% 2) {
      lines(rep(xlab[i], 2), y_leg + c(-0.9, -1.4))
      text(xlab[i], y_leg - 2.4, vlab[i], cex = cex_lab)
    }
  }
}

figs = list(a = 1:2, b = 3:4)
for (nm in names(figs)) {
  if (tikzPlot == 1) {
    tikz(paste0("tex/", if (nm == "a") "fig1" else "figS1", ".tex"), width = fig_w, height = fig_h, standAlone = TRUE,
         packages = c("\\usepackage{tikz}",
                      "\\usepackage[active,tightpage,psfixbb]{preview}",
                      "\\PreviewEnvironment{pgfpicture}",
                      "\\setlength\\PreviewBorder{0pt}",
                      "\\usepackage{amssymb}",
                      "\\usepackage{bm}","\\usepackage{amsthm}","\\usepackage{amsbsy}",
                      "\\usepackage{amsfonts}"))
  }
  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(xlo, xhi), ylim = c(ylo, yhi), axes = FALSE, xlab = "", ylab = "")
  par(xpd = TRUE)
  draw_panel(figs[[nm]])
  draw_legend()
  if (tikzPlot == 1) dev.off()
}
