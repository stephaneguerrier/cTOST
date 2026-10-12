rm(list=ls())

stopifnot(requireNamespace("cTOST",quietly=TRUE),
          as.character(packageVersion("cTOST"))=="1.1.0")
if (!require("scales")) install.packages("scales")
library(cTOST)
library(scales)

tikzPlot = 0
texDir = paste0("tex/")

if (tikzPlot == 1) {
  if (!require("tikzDevice")) install.packages("tikzDevice")
  library(tikzDevice)
}

biodistribution = cTOST::biodistribution
X = biodistribution$X
Y = biodistribution$Y

##############################

delta = 0.15
pi_x = c(0.2, 0.8)
delta_l = pi_x - delta
delta_u = pi_x + delta
seed = 12345

M = 20
alphas = seq(0.01, 0.2, l=M)
B = 10000
sol_qtost = sol_alpha_qtost = array(NA, dim = c(2, 3, M))

for (alpha in alphas){
  
  tmp = matrix(NA, 2, 3)
  colnames(tmp) = c("ci_l", "pi_y_hat", "ci_u")
  rownames(tmp) = c("pi_x1", "pi_x2")
  res = res_star = tmp
  # qTOST
  qtost_sol = cTOST::qtost(x=X, y=Y, pi_x=pi_x, delta=delta,
                        alpha=alpha, method="unadjusted")
  # qTOST CI's
  res[,2] = qtost_sol$pi_y_hat
  res[,1] = qtost_sol$ci[,1]
  res[,3] = qtost_sol$ci[,2]
  # alpha-qTOST
  alpha_qtost_sol = cTOST::qtost(x=X, y=Y, pi_x=pi_x, delta=delta,
                          alpha=alpha, method="alpha", B=B, seed=seed,
                          tol=.Machine$double.eps^0.5, max_iter=5,
                          tolpower=1e-3, MC_sup=TRUE)
  alpha_star = alpha_qtost_sol$corrected_alpha
  cat(paste0("Adjusted alpha: ", round(alpha_star, 4), "\n"))
  # alpha-qTOST CI's
  print(as.integer(all(alpha_qtost_sol$decision)))
  res_star[,2] = alpha_qtost_sol$pi_y_hat
  res_star[,1] = alpha_qtost_sol$ci[,1]
  res_star[,3] = alpha_qtost_sol$ci[,2]

  sol_qtost[,,which(alphas %in% alpha)] = res
  sol_alpha_qtost[,,which(alphas %in% alpha)] = res_star
}

#########
#plots 
#########

cexsize = 1.9
mylwd = 2
main_color = "black"
qtost_color = main_color
eq_margin_color = "#FF9900"
grid_color = "#DDDDDD"
vlin_col = "darkred"

qtost_fill = alpha("darkblue", 0.075)
alpha_qtost_fill = qtost_fill
qtost_fill_equiv = alpha("darkblue", 0.25)
alpha_qtost_fill_equiv = qtost_fill_equiv
eq_region_fill = alpha(eq_margin_color, 0.15)

if (tikzPlot == 1) {
  tikz(paste0(texDir, "fig8.tex"), 
       width = 13, 
       height = 6, 
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

par(mfrow = c(1, 2),
    mai = c(0.4, 0.1, 0.4, 0.2),
    omi = c(1, 0.8, 0, 0.05))

plot(NA, NA, xlim = range(alphas), ylim = c(0, 1), axes = FALSE,
     ylab = "", xlab = "", main="qTOST", cex.main=cexsize,
     panel.first = {
       grid(col = grid_color, lty = 1)
     })
axis(1, at = seq(0, 0.2, by = 0.05), cex.axis = cexsize-0.2)
axis(2, at = seq(0, 1, by = 0.2), cex.axis = cexsize-0.2)
mtext(paste0("Confidence interval for $\\pi_y$"), cex = cexsize+0.3, line = 3.7, side=2)
mtext(paste0("$\\alpha$"), cex = cexsize+0.3, line = 2.95, side=1)
polygon(c(alphas, rev(alphas)),
        c(rep(delta_l[1], M), rev(rep(delta_u[1], M))),
        col = eq_region_fill,
        border = NA)
polygon(c(alphas, rev(alphas)),
        c(rep(delta_l[2], M), rev(rep(delta_u[2], M))),
        col = eq_region_fill,
        border = NA)
lines(alphas, rep(delta_l[1], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, rep(delta_u[1], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, rep(delta_l[2], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, rep(delta_u[2], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, sol_qtost[1, 1, ], lwd=2, col=1, lty=2)
lines(alphas, sol_qtost[1, 3, ], lwd=2, col=1, lty=2)
lines(alphas, sol_qtost[1, 2, ], lwd=2, col=1, lty=1)
polygon(c(alphas, rev(alphas)),
        c(sol_qtost[1, 1, ], rev(sol_qtost[1, 3, ])),
        col = qtost_fill, border = NA)
qnt_sel = 1
inds_eq = sol_qtost[qnt_sel,1,] > (pi_x[qnt_sel]-delta) & sol_qtost[qnt_sel,3,] < (pi_x[qnt_sel]+delta)
polygon(c(alphas[inds_eq], rev(alphas[inds_eq])),
        c(sol_qtost[qnt_sel, 1, inds_eq], rev(sol_qtost[qnt_sel, 3, inds_eq])),
        col = qtost_fill_equiv,
        border = NA)
lines(alphas, sol_qtost[2, 1, ], lwd=2, col=1, lty=2)
lines(alphas, sol_qtost[2, 3, ], lwd=2, col=1, lty=2)
lines(alphas, sol_qtost[2, 2, ], lwd=2, col=1, lty=1)
polygon(c(alphas, rev(alphas)),
        c(sol_qtost[2, 1, ], rev(sol_qtost[2, 3, ])),
        col = qtost_fill, border = NA)
qnt_sel = 2
inds_eq = sol_qtost[qnt_sel,1,] > (pi_x[qnt_sel]-delta) & sol_qtost[qnt_sel,3,] < (pi_x[qnt_sel]+delta)
polygon(c(alphas[inds_eq], rev(alphas[inds_eq])),
        c(sol_qtost[qnt_sel, 1, inds_eq], rev(sol_qtost[qnt_sel, 3, inds_eq])),
        col = qtost_fill_equiv,
        border = NA)
abline(v = 0.05, col = vlin_col, lty=3, lwd=3)
abline(v = 0.1, col = vlin_col, lty=3, lwd=3)
box()

#####
# panel 2
#####

plot(NA, NA, xlim = range(alphas), ylim = c(0, 1), axes = FALSE,
     ylab = "", xlab = "", main="$\\bm{\\alpha}$-qTOST", cex.main=cexsize,
     panel.first = {
       grid(col = grid_color, lty = 1)
     })
axis(1, at = seq(0, 0.2, by = 0.05), cex.axis = cexsize-0.2)
mtext(paste0("$\\alpha$"), cex = cexsize+0.3, line = 2.95, side=1)
polygon(c(alphas, rev(alphas)),
        c(rep(delta_l[1], M), rev(rep(delta_u[1], M))),
        col = eq_region_fill,
        border = NA)
polygon(c(alphas, rev(alphas)),
        c(rep(delta_l[2], M), rev(rep(delta_u[2], M))),
        col = eq_region_fill,
        border = NA)
lines(alphas, rep(delta_l[1], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, rep(delta_u[1], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, rep(delta_l[2], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, rep(delta_u[2], M), col = eq_margin_color, lwd = 2, lty=1)
lines(alphas, sol_alpha_qtost[1, 1, ], lwd=2, col=1, lty=2)
lines(alphas, sol_alpha_qtost[1, 3, ], lwd=2, col=1, lty=2)
lines(alphas, sol_alpha_qtost[1, 2, ], lwd=2, col=1, lty=1)
polygon(c(alphas, rev(alphas)),
        c(sol_alpha_qtost[1, 1, ], rev(sol_alpha_qtost[1, 3, ])),
        col = alpha_qtost_fill, border = NA)
qnt_sel = 1
inds_eq = sol_alpha_qtost[qnt_sel,1,] > (pi_x[qnt_sel]-delta) & sol_alpha_qtost[qnt_sel,3,] < (pi_x[qnt_sel]+delta)
polygon(c(alphas[inds_eq], rev(alphas[inds_eq])),
        c(sol_alpha_qtost[qnt_sel, 1, inds_eq], rev(sol_alpha_qtost[qnt_sel, 3, inds_eq])),
        col = alpha_qtost_fill_equiv,
        border = NA)
lines(alphas, sol_alpha_qtost[2, 1, ], lwd=2, col=1, lty=2)
lines(alphas, sol_alpha_qtost[2, 3, ], lwd=2, col=1, lty=2)
lines(alphas, sol_alpha_qtost[2, 2, ], lwd=2, col=1, lty=1)
polygon(c(alphas, rev(alphas)),
        c(sol_alpha_qtost[2, 1, ], rev(sol_alpha_qtost[2, 3, ])),
        col = alpha_qtost_fill, border = NA)
qnt_sel = 2
inds_eq = sol_alpha_qtost[qnt_sel,1,] > (pi_x[qnt_sel]-delta) & sol_alpha_qtost[qnt_sel,3,] < (pi_x[qnt_sel]+delta)
polygon(c(alphas[inds_eq], rev(alphas[inds_eq])),
        c(sol_alpha_qtost[qnt_sel, 1, inds_eq], rev(sol_alpha_qtost[qnt_sel, 3, inds_eq])),
        col = alpha_qtost_fill_equiv,
        border = NA)
abline(v = 0.05, col = vlin_col, lty=3, lwd=3)
abline(v = 0.1, col = vlin_col, lty=3, lwd=3)
box()

leg_nam = c("$\\widehat\\pi_y$", 
            "Boundary CI for $\\pi_y$", 
            "Equivalence Region", 
            "Cannot reject H$_0$", 
            "Accept H$_1$")
leg_col = c(main_color, 
            qtost_color, 
            eq_region_fill, 
            alpha_qtost_fill, 
            alpha_qtost_fill_equiv)
custom_widths = c(0.02, 0.085, 0.085, 0.076, 0.06)
legend("bottomleft",
       legend = leg_nam,
       col = leg_col,
       y.intersp = 1,
       x.intersp = c(1.4, 1.4, 0.5, 0.5, 0.5)/2,
       xjust=4,
       yjust=0,
       text.width=custom_widths,
       lty = c(1, 2, 0, 0, 0), 
       seg.len = 1.8,
       pch = c(NA, NA, 15, 15, 15),
       pt.cex = 3,
       lwd = c(mylwd+1, mylwd+1, NA, NA, NA),
       pt.bg = c(NA, NA, eq_region_fill, alpha_qtost_fill, alpha_qtost_fill_equiv),
       cex = cexsize -0.1,
       bty = "n",
       xpd = NA,
       inset = c(-1.075, -0.35),
       ncol = length(leg_nam))

if (tikzPlot == 1) {
  dev.off()
}

#########
# plot 2
#########

cexsize = 2
cols = c("#501d8a", "#ee8c7d")
yall = c(2.55, 2.4)
lwds = c(2.5, 2.5)*3
cap_height = 0.02  
dim_nam = c("qTOST", "$\\bm{\\alpha}$-qTOST")

if (tikzPlot == 1) {
  tikz(paste0(texDir, "fig7.tex"), 
       width = 13.5, 
       height = 6, 
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

par(mfrow = c(2, 2),
    mai = c(0.85, 1.7, 0.2, 0.1),
    omi = c(0, 0, 0, 0.05))

for (u in 1:2){
  if (u==1){
    res_star = sol_alpha_qtost[,,which(abs(alphas - 0.05) < 1e-10)]
    res = sol_qtost[,,which(abs(alphas - 0.05) < 1e-10)]
    alp_print = 5
  } else {
    res_star = sol_alpha_qtost[,,which(abs(alphas - 0.1) < 1e-10)]
    res = sol_qtost[,,which(abs(alphas - 0.1) < 1e-10)]
    alp_print = 10
  }  
  
  qnt = 1
  x1 = seq(0, 0.6, length.out = 11)
  plot(NA, NA, xlim = range(x1), ylim = c(2.3, 2.7), axes = FALSE,
       ylab = "", xlab = "")
  polygon(x=c(pi_x[qnt]-delta, pi_x[qnt]+delta, pi_x[qnt]+delta, pi_x[qnt]-delta),
          y=rep(c(2.1, max(yall)+0.1), each=2),
          col = eq_region_fill,
          border = NA)
  abline(h = yall, lty = 1, col = "grey")
  lines(c(pi_x[qnt],pi_x[qnt]), c(2.1, max(yall)+0.1), lwd=4,lty=3,col="gray")
  mtext(paste0("$\\bm{\\alpha=", alp_print, "\\%", "}$, $\\bm{\\pi_x=", pi_x[qnt], "}$"), at=mean(x1), side=3, line=-1.5, cex=cexsize)
  lines(c(pi_x[qnt]-delta, pi_x[qnt]-delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
  lines(c(pi_x[qnt]+delta, pi_x[qnt]+delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
  axis(1,  cex.axis=cexsize+0.3, padj = 0.5)
  mtext("$\\pi_y$", side = 1, line=4.1, cex=cexsize+0.3)
  ind = 1
  lines(res[qnt,-2], c(yall[ind], yall[ind]), 
        lwd = lwds[ind], col = cols[ind])
  lines(rep(res[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  lines(rep(res[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  points(res[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3.5)
  mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)  
  ind = 2
  lines(res_star[qnt,-2], c(yall[ind], yall[ind]), 
        lwd = lwds[ind], col = cols[ind])
  lines(rep(res_star[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  lines(rep(res_star[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  points(res_star[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3.5)
  mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)
  
  #########
  # panel 2
  #########
  qnt = 2
  x1 = seq(0.4, 1, length.out = 7)
  plot(NA, NA, xlim = range(x1), ylim = c(2.3, 2.7), axes = FALSE,
       ylab = "", xlab = "")
  polygon(x=c(pi_x[qnt]-delta, pi_x[qnt]+delta, pi_x[qnt]+delta, pi_x[qnt]-delta),
          y=rep(c(2.1, max(yall)+0.1), each=2),
          col = eq_region_fill,
          border = NA)
  abline(h = yall, lty = 1, col = "grey")
  lines(c(pi_x[qnt],pi_x[qnt]), c(2.1, max(yall)+0.1), lwd=4,lty=3,col="gray")
  mtext(paste0("$\\bm{\\alpha=", alp_print, "\\%", "}$, $\\bm{\\pi_x=", pi_x[qnt], "}$"), at=mean(x1), side=3, line=-1.5, cex=cexsize)
  lines(c(pi_x[qnt]-delta, pi_x[qnt]-delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
  lines(c(pi_x[qnt]+delta, pi_x[qnt]+delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
  axis(1,  cex.axis=cexsize+0.3, padj = 0.5, at = x1)
  mtext("$\\pi_y$", side = 1, line=4.1, cex=cexsize+0.3)
  ind = 1
  lines(res[qnt,-2], c(yall[ind], yall[ind]), 
        lwd = lwds[ind], col = cols[ind])
  lines(rep(res[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  lines(rep(res[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  points(res[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3.5)
  mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)
  ind = 2
  lines(res_star[qnt,-2], c(yall[ind], yall[ind]), 
        lwd = lwds[ind], col = cols[ind])
  lines(rep(res_star[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  lines(rep(res_star[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
        col = cols[ind], lwd = lwds[ind])
  points(res_star[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3.5)
  mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)
}


if (tikzPlot == 1) {
  dev.off()
}
