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

# https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/021814s030lbl.pdf
# C_trough: table 5
# reference male group
mu_x_orig = 35.6
sd_x_orig = 16.7
n_x = 106
# target female group
mu_y_orig = 41.6
sd_y_orig = 24.3
n_y = 14

# data transformation
mu_x = log(mu_x_orig^2 / sqrt(mu_x_orig^2 + sd_x_orig^2))
sd_x = sqrt(log(1 + (sd_x_orig^2 / mu_x_orig^2)))
mu_y = log(mu_y_orig^2 / sqrt(mu_y_orig^2 + sd_y_orig^2))
sd_y = sqrt(log(1 + (sd_y_orig^2 / mu_y_orig^2)))

##############################

alpha = 0.05
delta = 0.1
B = 100000
seed = 12345

pi_x = seq(from = 0.11, to = 0.3, by = 0.001)
M = length(pi_x)

res = res_star = matrix(NA, M, 3)
for (i in 1:M){
  
  # qTOST CI's
  qtost_sol = cTOST::qtost(x = list(mean=mu_x, sd=sd_x, n=n_x),
                        y = list(mean=mu_y, sd=sd_y, n=n_y),
                        pi_x=pi_x[i], delta=delta, alpha=alpha, method="unadjusted")
  res[i,2] = qtost_sol$pi_y_hat
  res[i,1] = qtost_sol$ci[1,1]
  res[i,3] = qtost_sol$ci[1,2]
  
  # alpha-qTOST
  alpha_qtost_sol = cTOST::qtost(x = list(mean=mu_x, sd=sd_x, n=n_x),
                          y = list(mean=mu_y, sd=sd_y, n=n_y),
                          pi_x=pi_x[i], delta=delta, alpha=alpha, method="alpha",
                          B=B, seed=seed, tol=.Machine$double.eps)
  alpha_star = alpha_qtost_sol$corrected_alpha
  cat(paste0("Adjusted alpha: ", round(alpha_star, 3), "\n"))
  # alpha-qTOST CI's
  res_star[i,2] = alpha_qtost_sol$pi_y_hat
  res_star[i,1] = alpha_qtost_sol$ci[1,1]
  res_star[i,3] = alpha_qtost_sol$ci[1,2]
}

#########
#plots 
#########

main_color = "black"
qtost_color = main_color
alpha_qtost_color = main_color
eq_margin_color = "#FF9900"
grid_color = "#DDDDDD"
vlin_col = "darkred"

qtost_fill = alpha("darkblue", 0.075)
alpha_qtost_fill = qtost_fill
qtost_fill_equiv = alpha("darkblue", 0.25)
alpha_qtost_fill_equiv = qtost_fill_equiv
eq_region_fill = alpha(eq_margin_color, 0.15)

if (tikzPlot == 1) {
  tikz(paste0(texDir, "fig5.tex"), 
       width = 13.5, 
       height = 3, 
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
    mai = c(0.9, 1.7, 0, 0.1),
    omi = c(0, 0, 0, 0.05))


cexsize = 2
cols = c("#501d8a", "#ee8c7d")
yall = c(2.55, 2.4)
lwds = c(2.5, 2.5)*3
cap_height = 0.02  
dim_nam = c("qTOST", "$\\bm{\\alpha}$-qTOST")

qnt = which(pi_x == 0.15)
x1 = seq(0, 0.35, length.out = 11)
plot(NA, NA, xlim = range(x1), ylim = c(2.3, 2.7), axes = FALSE, ylab = "", xlab = "")
polygon(x=c(pi_x[qnt]-delta, pi_x[qnt]+delta, pi_x[qnt]+delta, pi_x[qnt]-delta),
        y=rep(c(2.1, max(yall)+0.1), each=2),
        col = eq_region_fill,
        border = NA)
abline(h = yall, lty = 1, col = "grey")
lines(c(pi_x[qnt],pi_x[qnt]), c(2.1, max(yall)+0.1), lwd=4,lty=3,col="gray")
mtext(paste0("$\\bm{\\pi_x=", pi_x[qnt], "}$"), at=c(pi_x[qnt]), side=3, line=-1.5, cex=cexsize)
lines(c(pi_x[qnt]-delta, pi_x[qnt]-delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
lines(c(pi_x[qnt]+delta, pi_x[qnt]+delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
axis(1, padj = 0.5, labels=NA, lwd.tick=0)
axis(1, at = seq(0.05, 0.35, l=4), labels=seq(0.05, 0.35, l=4), cex.axis=cexsize-0.1, padj = 0.5)
mtext("$\\pi_y$", side = 1, line=3.5, cex=cexsize+0.3)
ind = 1
lines(res[qnt,-2], c(yall[ind], yall[ind]), 
      lwd = lwds[ind], col = cols[ind])
lines(rep(res[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
lines(rep(res[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
points(res[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3)
mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)
ind = 2
lines(res_star[qnt,-2], c(yall[ind], yall[ind]), 
      lwd = lwds[ind], col = cols[ind])
lines(rep(res_star[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
lines(rep(res_star[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
points(res_star[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3)
mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)

########################################################

qnt = which(pi_x == 0.2)
x1 = seq(0, 0.4, length.out = 11)
plot(NA, NA, xlim = range(x1), ylim = c(2.3, 2.7), axes = FALSE,
     ylab = "", xlab = "")
polygon(x=c(pi_x[qnt]-delta, pi_x[qnt]+delta, pi_x[qnt]+delta, pi_x[qnt]-delta),
        y=rep(c(2.1, max(yall)+0.1), each=2),
        col = eq_region_fill,
        border = NA)
abline(h = yall, lty = 1, col = "grey")
lines(c(pi_x[qnt],pi_x[qnt]), c(2.1, max(yall)+0.1), lwd=4,lty=3,col="gray")
mtext(paste0("$\\bm{\\pi_x=", pi_x[qnt], "}$"), at=c(pi_x[qnt]), side=3, line=-1.5, cex=cexsize)
lines(c(pi_x[qnt]-delta, pi_x[qnt]-delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
lines(c(pi_x[qnt]+delta, pi_x[qnt]+delta), c(2.1, max(yall)+0.1), lty = 3, col = "black", lwd=4)
axis(1, at = seq(0, 0.4, by = 0.1), labels=seq(0, 0.4, by = 0.1), cex.axis=cexsize-0.1, padj = 0.5)
mtext("$\\pi_y$", side = 1, line=3.5, cex=cexsize+0.3)
ind = 1
lines(res[qnt,-2], c(yall[ind], yall[ind]), 
      lwd = lwds[ind], col = cols[ind])
lines(rep(res[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
lines(rep(res[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
points(res[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3)
mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)
ind = 2
lines(res_star[qnt,-2], c(yall[ind], yall[ind]), 
      lwd = lwds[ind], col = cols[ind])
lines(rep(res_star[qnt,1], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
lines(rep(res_star[qnt,3], 2), c(yall[ind] - cap_height, yall[ind] + cap_height), 
      col = cols[ind], lwd = lwds[ind])
points(res_star[qnt,2], yall[ind], col = cols[ind], pch = 16, cex = 3)
mtext(dim_nam[ind], side = 2, at = yall[ind], adj = 1.1, las = 1, cex=cexsize, font=2)

if (tikzPlot == 1) {
  dev.off()
}

###################################################

if (tikzPlot == 1) {
  tikz(paste0(texDir, "fig6.tex"), 
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

cexsize = 1.9
mylwd = 2

y_max = 0.5
y_min = 0

# Plot 1: CI qTOST
plot(pi_x, res[,2], type = "l", 
     ylim = c(y_min, y_max), 
     xlim = range(pi_x),
     xlab = NA,
     ylab = NA,
     col = main_color, 
     lwd = mylwd,
     cex.main=cexsize,
     cex.axis=cexsize-0.2,
     main = "qTOST",
     panel.first = {
       grid(col = grid_color, lty = 1)
     })
abline(v=0.2, lty=3, lwd=3, col=vlin_col)
abline(v=0.15, lty=3, lwd=3, col=vlin_col)
mtext(paste0("Confidence interval for $\\pi_y$"), cex = cexsize+0.3, line = 3.7, side=2)
mtext(paste0("$\\pi_x$"), cex = cexsize+0.3, line = 2.95, side=1)
# equivalence region
lines(pi_x, pi_x + delta, col = eq_margin_color, lwd = mylwd * 1.2)
lines(pi_x, pi_x - delta, col = eq_margin_color, lwd = mylwd * 1.2)
polygon(c(pi_x, rev(pi_x)),
        c(pi_x + delta, rev(pi_x - delta)),
        col = eq_region_fill,
        border = NA)
# CI qTOST
lines(pi_x, res[,1], lty = 2, lwd = mylwd, col = qtost_color)
lines(pi_x, res[,3], lty = 2, lwd = mylwd, col = qtost_color)
polygon(c(pi_x, rev(pi_x)), c(res[,1], rev(res[,3])),
        col = qtost_fill, border = NA)
polygon(c(pi_x[res[,1] > (pi_x-delta) & res[,3] < (pi_x+delta)],
          rev(pi_x[res[,1] > (pi_x-delta) & res[,3] < (pi_x+delta)])),
        c(res[res[,1] > (pi_x-delta) & res[,3] < (pi_x+delta), 1],
          rev(res[res[,1] > (pi_x-delta) & res[,3] < (pi_x+delta), 3])),
        col = qtost_fill_equiv,
        border = NA)


# Plot 2: CI alpha-qTOST
plot(pi_x, res[,2], type = "l", 
     ylim = c(y_min, y_max), 
     xlim = range(pi_x),
     xlab = NA, 
     ylab = NA,
     col = main_color, 
     lwd = mylwd,
     cex.main=cexsize,
     cex.axis=cexsize-0.2,
     yaxt = "n",
     main = "$\\bm\\alpha$-qTOST",
     panel.first = {
       grid(col = grid_color, lty = 1)
     })

abline(v=0.2, lty=3, lwd=3, col=vlin_col)
abline(v=0.15, lty=3, lwd=3, col=vlin_col)
mtext(paste0("$\\pi_x$"), cex = cexsize+0.3, line = 2.95,side=1)

# equivalence region
lines(pi_x, pi_x + delta, col = eq_margin_color, lwd = mylwd * 1.2)
lines(pi_x, pi_x - delta, col = eq_margin_color, lwd = mylwd * 1.2)
polygon(c(pi_x, rev(pi_x)),
        c(pi_x + delta, rev(pi_x - delta)),
        col = eq_region_fill,
        border = NA)

# CI alpha-qTOST
lines(pi_x, res_star[,1], lty = 2, col = alpha_qtost_color, lwd = mylwd)
lines(pi_x, res_star[,3], lty = 2, col = alpha_qtost_color, lwd = mylwd)
polygon(c(pi_x, rev(pi_x)), c(res_star[,1], rev(res_star[,3])),
        col = qtost_fill, border = NA)
polygon(c(pi_x[res_star[,1] > (pi_x-delta) & res_star[,3] < (pi_x+delta)],
          rev(pi_x[res_star[,1] > (pi_x-delta) & res_star[,3] < (pi_x+delta)])),
        c(res_star[res_star[,1] > (pi_x-delta) & res_star[,3] < (pi_x+delta), 1],
          rev(res_star[res_star[,1] > (pi_x-delta) & res_star[,3] < (pi_x+delta), 3])),
        col = alpha_qtost_fill_equiv, 
        border = NA)

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
