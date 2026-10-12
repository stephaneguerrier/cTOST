rm(list=ls())

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

texDir = paste0("tex/")

alpha = 0.05
require(grDevices)
color = palette.colors(n = 7, palette = "Okabe-Ito")
cols = color[c(6, 7)]

for (pick_gamma in 1:length(dimnames(sol_qtost_avg)[[4]])) {
  for (pick_l in 1:length(dimnames(sol_qtost_avg)[[2]])) {

    if (tikzPlot==1) {
      mypath_output = paste0(texDir, c("figS2", "figS3", "figS4", "figS5", "figS6", "fig2")[(pick_gamma-1)*3+pick_l], ".tex")
      tikz(file=mypath_output,standAlone=TRUE,width=7,height=5)
    }

    par(mfrow = c(4, 4),
        mai = c(0.1, 0.1, 0, 0),
        omi = c(0.65, 0.8, 0.27, 0.05))

    pick_pi_x_range = 1:4
    pick_n_range = 1:4
    cexsize = 1
    for (pick_n in pick_n_range){
      for (pick_pi_x in pick_pi_x_range) {

        assignments = strsplit(dimnames(sol_qtost_avg)[[3]][pick_pi_x], ";")[[1]]
        pi_x = as.numeric(sub(".*=", "", assignments[1]))
        delta = as.numeric(sub(".*=", "", assignments[2]))

        delta_u = pi_x + delta
        delta_l = pi_x - delta

        n = as.numeric(sub(".*=", "", dimnames(sol_qtost_avg)[[1]][pick_n]))
        pi_y_grid = pnorm(seq(qnorm(pi_x-delta*1.2), qnorm(pi_x+delta*1.2), len=50))
        max_y = max(sol_alpha_qtost_avg[pick_n, pick_l, , pick_gamma, ])+4
        plot(NA,t="l",xlim=range(pi_y_grid),ylim=c(0,max_y),
             xlab=NA, ylab=NA, main=NA, axes=F)
        grid()
        box()

        lines(pi_y_grid, sol_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, ],
              col=cols[1], lwd=2)
        lines(pi_y_grid, sol_alpha_qtost_avg[pick_n, pick_l, pick_pi_x, pick_gamma, ],
              col=cols[2], lwd=2)

        abline(h=alpha*100, lty=2, col="black")
        abline(v=delta_u, lty=2, col="black")
        abline(v=delta_l, lty=2, col="black")
        abline(v=pi_x, lty=2, col="black")

        if (pick_pi_x == pick_pi_x_range[1]){
          axis(2, cex.axis=cexsize+0.1)
          mtext(paste0("$ n_y=", n, "$"), cex = cexsize, line = 2.5, side=2)
        }
        if (pick_n == pick_n_range[1]){
          mtext(paste0("$\\pi_x=", pi_x, "$ ($c=", delta, "$)"), cex = cexsize, line = 0.6, side=3)
        }
        if (pick_n == pick_n_range[length(pick_n_range)]){
          axis(1, cex.axis=cexsize+0.1)
        }
      }
    }

    mtext("Probability of rejecting H$_0$ (in \\%)", cex = cexsize+0.1, line = 4.6, side=2, outer=T)
    mtext("$\\pi_y$", cex = cexsize+0.1, side=1, line = 2.5, outer=TRUE)

    ltys = rep(1, 5)
    par(fig = c(0, 1, 0, 1),
        oma = c(0, 0, 0, 0),
        mar = c(0, 0, 0, 0),
        omi = c(0, 0.8, 0.27, 0.05),
        new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    par(xpd=TRUE)
    est_names = c(" qTOST", "$ \\alpha$-qTOST")
    legend("bottom",
           est_names,
           col = cols[1:2],
           lty = ltys, lwd = 2, cex = cexsize+0.2,
           bty = "n",
           ncol=1, xpd=NA,
           y.intersp=0,
           x.intersp=0.5,
           inset=c(0.1, 0),
           border = "grey",
           xjust=2,
           seg.len=2, horiz = T)
    par(xpd=F)

    if (tikzPlot==1) {
      dev.off()
    }
  }
}
