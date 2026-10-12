rm(list = ls())

tikzPlot = 0
if (tikzPlot == 1) {
  if (!require("tikzDevice")) install.packages("tikzDevice")
  library(tikzDevice)
}

n_simu = 50000
n_array = 5000
loadDir = "simulations/HPC/sim_univ/sim_HPC_out/sim_setting_1"
name_sim = paste0("/results_sim_setting_size_1_n_simu_", n_simu, "_n_array_", n_array)
load(file = paste0(loadDir, name_sim, ".rdata"))
sol_qtost_avg = sol_TOST_avg
sol_alpha_qtost_avg = sol_alphaTOST_avg

texDir = paste0("tex/fig3.tex")

sol = sol_alpha_qtost_avg - sol_qtost_avg

m = 50
r1 = colorRampPalette(c("#1f618d", "#eaf2f8"))
col1 = r1(m)
r2 = colorRampPalette(c("#f9ebea", "#c0392b"))
col2 = r2(m)
q1 = seq(from = 0, to = 5, length.out = length(col1))
q2 = seq(from = 5, to = 10, length.out = length(col2))

n= 50
m = prod(dim(sol)[1:4])

if (tikzPlot == 1) {
  tikz(texDir,
       width = 10,
       height = 5.2,
       standAlone = TRUE,
       packages = c(
         "\\usepackage{tikz}",
         "\\usepackage[active,tightpage,psfixbb]{preview}",
         "\\PreviewEnvironment{pgfpicture}",
         "\\setlength\\PreviewBorder{0pt}",
         "\\usepackage{amssymb}",
         "\\usepackage{bm}",
         "\\usepackage{amsthm}",
         "\\usepackage{amsbsy}",
         "\\usepackage{amsfonts}"
       ))
}

m = m/2
par(mar = c(0, 0.2, 0, 0.2))
plot(NA, ylim = c(-3.65, 1.5), xlim = c(0, n-1), axes = FALSE, xlab = "", ylab = "")

pick_pi_x_range = 1:2
pick_l_range = 1:length(dimnames(sol_qtost_avg)[[2]])
pick_gamma_range = 1:length(dimnames(sol_qtost_avg)[[4]])
pick_n_range = 1:length(dimnames(sol_qtost_avg)[[1]])

i = 1
for (pick_pi_x in pick_pi_x_range) {
  for (pick_gamma in pick_gamma_range) {
    for (pick_l in pick_l_range) {
      for (pick_n in pick_n_range) {
        assignments = strsplit(dimnames(sol_qtost_avg)[[3]], ";")
        pi_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 1)))
        delta_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 2)))
        size_alpha_qtost = max(sol_alpha_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 1],
                       sol_alpha_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 2])
        size_qtost = max(sol_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 1],
                      sol_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 2])
        alpha = 5
        j = 1.2
        uy = 0.1
        if (size_alpha_qtost > alpha) {
          rect(2 + i - 0.5, j - uy, 2 + i + 0.5, j + uy,
               col = col2[which.min(abs(q2 - size_alpha_qtost))], border = "white")
        } else {
          rect(2 + i - 0.5, j - uy, 2 + i + 0.5, j + uy,
               col = col1[which.min(abs(q1 - size_alpha_qtost))], border = "white")
        }
        j = 1.4
        if (size_qtost > alpha) {
          rect(2 + i - 0.5, j - uy, 2 + i + 0.5, j + uy,
               col = col2[which.min(abs(q2 - size_qtost))], border = "white")
        } else {
          rect(2 + i - 0.5, j - uy, 2 + i + 0.5, j + uy,
               col = col1[which.min(abs(q1 - size_qtost))], border = "white")
        }
        i = i + 1
      }
    }
  }
}

my_cex = 1
lab_y = j - 0.5
lab_x = 1.8

text(lab_x, lab_y + 0.05, "$n_y$", cex = my_cex + 0.2)
text(lab_x + 0.4, 1.15, "$\\alpha$-qTOST", cex = my_cex + 0.2, adj = 1)
text(lab_x + 0.4, 1.35, "qTOST", cex = my_cex + 0.2, adj = 1)

ny_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[1]]))
off = 50

for (i in 0:(m-1)) {
  lines(-c(i,i) + off, c(lab_y + 0.2, lab_y + 0.15))
  if (rev(ny_vec)[((i) %% 4) + 1] == 120) {
    text(-c(i,i) + off + 1, lab_y, rev(ny_vec)[((i) %% 4) + 1],
         cex = my_cex - 0.1, pos = 2)
  } else {
    text(-c(i,i) + off + 0.85, lab_y, rev(ny_vec)[((i) %% 4) + 1],
         cex = my_cex - 0.1, pos = 2)
  }
}

arrow1x = 0.7
label_offset = 0.09
text(lab_x, arrow1x, "$l$", cex = my_cex + 0.2)

l_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[2]]))
for (i in 1:(6*2)) {
  j = i - 1
  arrows(50 - 4*(i-1), arrow1x, 47 - 4*(i-1), arrow1x, code = 3, length = 0.03)
  if (rev(l_vec)[(j %% 3) + 1] > 1){
    tt = paste0("$1/", rev(l_vec)[(j %% 3) + 1], "$")
  } else {
    tt = rev(l_vec)[(j %% 3) + 1]
  }
  text(mean(c(50 - 4*(i-1), 47 - 4*(i-1))), arrow1x - label_offset,
       tt, cex = my_cex, srt = 0)
}

arrow1x = arrow1x - 0.25
text(lab_x, arrow1x, "$\\gamma$", cex = my_cex + 0.2)

gamma_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[4]]))
for (i in 1:4) {
  j = i - 1
  arrows(50 - 3*4*(i-1), arrow1x, 39 - 3*4*(i-1), arrow1x, code = 3, length = 0.03)
  text(mean(c(50 - 3*4*(i-1), 39 - 3*4*(i-1))), arrow1x - label_offset,
       rev(gamma_vec)[(j %% 2) + 1], cex = my_cex, srt = 0)
}

arrow1x = arrow1x - 0.25
text(lab_x, arrow1x, "$\\pi_x$", cex = my_cex + 0.2)

assignments = strsplit(dimnames(sol_qtost_avg)[[3]], ";")
pi_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 1)))
delta_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 2)))

for (i in 1:2) {
  j = i - 1
  arrows(50 - 2*3*4*(i-1), arrow1x, 27 - 2*3*4*(i-1), arrow1x, code = 3, length = 0.03)
  tt = paste0(rev(pi_vec[1:2])[(j %% 2) + 1], " ($c=", rev(delta_vec[1:2])[(j %% 2) + 1], "$)")
  text(mean(c(50 - 2*3*4*(i-1), 27 - 2*3*4*(i-1))), arrow1x - label_offset,
       tt, cex = my_cex, srt = 0)
}

deltay = 2
pick_pi_x_range = 3:4
i = 1

for (pick_pi_x in pick_pi_x_range) {
  for (pick_gamma in pick_gamma_range) {
    for (pick_l in pick_l_range) {
      for (pick_n in pick_n_range) {
        assignments = strsplit(dimnames(sol_qtost_avg)[[3]], ";")
        pi_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 1)))
        delta_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 2)))
        size_alpha_qtost = max(sol_alpha_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 1],
                       sol_alpha_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 2])
        size_qtost = max(sol_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 1],
                      sol_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, 2])
        alpha = 5
        j = 1.2
        uy = 0.1
        if (size_alpha_qtost > alpha) {
          rect(2 + i - 0.5, j - uy - deltay, 2 + i + 0.5, j + uy - deltay,
               col = col2[which.min(abs(q2 - size_alpha_qtost))], border = "white")
        } else {
          rect(2 + i - 0.5, j - uy - deltay, 2 + i + 0.5, j + uy - deltay,
               col = col1[which.min(abs(q1 - size_alpha_qtost))], border = "white")
        }
        j = 1.4
        if (size_qtost > alpha) {
          rect(2 + i - 0.5, j - uy - deltay, 2 + i + 0.5, j + uy - deltay,
               col = col2[which.min(abs(q2 - size_qtost))], border = "white")
        } else {
          rect(2 + i - 0.5, j - uy - deltay, 2 + i + 0.5, j + uy - deltay,
               col = col1[which.min(abs(q1 - size_qtost))], border = "white")
        }
        i = i + 1
      }
    }
  }
}

lab_y = j - 0.5 - deltay
text(lab_x, lab_y + 0.05, "$n_y$", cex = my_cex + 0.2)
text(lab_x + 0.4, 1.15 - deltay, "$\\alpha$-qTOST", cex = my_cex + 0.2, adj = 1)
text(lab_x + 0.4, 1.35 - deltay, "qTOST", cex = my_cex + 0.2, adj = 1)

ny_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[1]]))
off = 50

for (i in 0:(m-1)) {
  lines(-c(i,i) + off, c(lab_y + 0.2, lab_y + 0.15))
  if (rev(ny_vec)[((i) %% 4) + 1] == 120) {
    text(-c(i,i) + off + 1, lab_y, rev(ny_vec)[((i) %% 4) + 1],
         cex = my_cex - 0.1, pos = 2)
  } else {
    text(-c(i,i) + off + 0.85, lab_y, rev(ny_vec)[((i) %% 4) + 1],
         cex = my_cex, pos = 2)
  }
}

arrow1x = 0.7 - deltay
label_offset = 0.09
text(lab_x, arrow1x, "$l$", cex = my_cex + 0.2)

l_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[2]]))
for (i in 1:(6*2)) {
  j = i - 1
  arrows(50 - 4*(i-1), arrow1x, 47 - 4*(i-1), arrow1x, code = 3, length = 0.03)
  if (rev(l_vec)[(j %% 3) + 1] > 1){
    tt = paste0("$1/", rev(l_vec)[(j %% 3) + 1], "$")
  } else {
    tt = rev(l_vec)[(j %% 3) + 1]
  }
  text(mean(c(50 - 4*(i-1), 47 - 4*(i-1))), arrow1x - label_offset,
       tt, cex = my_cex, srt = 0)
}

arrow1x = arrow1x - 0.25
text(lab_x, arrow1x, "$\\gamma$", cex = my_cex + 0.2)

gamma_vec = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[4]]))
for (i in 1:4) {
  j = i - 1
  arrows(50 - 3*4*(i-1), arrow1x, 39 - 3*4*(i-1), arrow1x, code = 3, length = 0.03)
  text(mean(c(50 - 3*4*(i-1), 39 - 3*4*(i-1))), arrow1x - label_offset,
       rev(gamma_vec)[(j %% 2) + 1], cex = my_cex, srt = 0)
}

arrow1x = arrow1x - 0.25
text(lab_x, arrow1x, "$\\pi_x$", cex = my_cex + 0.2)

assignments = strsplit(dimnames(sol_qtost_avg)[[3]], ";")
pi_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 1)))
delta_vec = as.numeric(sub(".*=", "", sapply(assignments, `[[`, 2)))

for (i in 1:2) {
  j = i - 1
  arrows(50 - 2*3*4*(i-1), arrow1x, 27 - 2*3*4*(i-1), arrow1x, code = 3, length = 0.03)
  tt = paste0(rev(pi_vec[3:4])[(j %% 2) + 1], " ($c=", rev(delta_vec[3:4])[(j %% 2) + 1], "$)")
  text(mean(c(50 - 2*3*4*(i-1), 27 - 2*3*4*(i-1))), arrow1x - label_offset,
       tt, cex = my_cex, srt = 0)
}

alpha = 0.05
nrep = 50000

lower = qbinom(0.005, nrep, alpha)/nrep * 100
upper = qnorm(0.995, mean = alpha, sd = sqrt(alpha*(1-alpha))/sqrt(nrep)) * 100

q1 = seq(from = 0, to = 5, length.out = length(col1))
q2 = seq(from = 5, to = 10, length.out = length(col2))

legend_breaks = sort(c(0, 1, 2, 3.5, lower, upper, 5.5))
legend_colors = numeric(length(legend_breaks))
for(i in 1:length(legend_breaks)) {
  if(legend_breaks[i] <= 5) {
    idx = which.min(abs(q1 - legend_breaks[i]))
    legend_colors[i] = col1[idx]
  } else {
    idx = which.min(abs(q2 - legend_breaks[i]))
    legend_colors[i] = col2[idx]
  }
}

cols = legend_colors
mm = length(cols)

y_leg = -3.05
legend_width = 16
x_start = 25
myseq = seq(from = x_start - legend_width/2, to = x_start + legend_width/2, length.out = mm)
line_range = range(myseq) + c(-0.2, 0.2)

rect_width = myseq[which(legend_breaks >= upper)[1]] - myseq[which(legend_breaks <= lower)[length(which(legend_breaks <= lower))]]
rect_center = mean(c(myseq[which(legend_breaks >= upper)[1]], myseq[which(legend_breaks <= lower)[length(which(legend_breaks <= lower))]]))

rect(rect_center - rect_width, y_leg - 0.22,
     rect_center + rect_width, y_leg + 0.28,
     border = "grey80", col = "grey80")

text(x_start, y_leg - 0.55, "Empirical size (in \\%)", cex = my_cex + 0.1)

arrows(rect_center, y_leg + 0.54, rect_center, y_leg + 0.3, length = 0.06, col = "grey60", lwd = 1.2)
text(x_start, y_leg + 0.66,
     paste0("$99\\%$ Simulation error tolerance under H$_0$ ($\\alpha=0.05$)"),
     col = "grey30", cex = 0.95)

points(myseq, rep(y_leg, mm), pch = 15, col = cols, cex = 5)
hlin = 0.165
lines(line_range, c(y_leg, y_leg) - hlin)

for (i in 1:length(legend_breaks)) {
  lines(rep(myseq[i], 2), y_leg + c(-hlin, -hlin-0.04))
  label = round(legend_breaks[i], 2)
  text(myseq[i], y_leg - 0.3, label, cex = my_cex)
}

if (tikzPlot == 1) {
  dev.off()
}
