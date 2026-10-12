rm(list=ls())

main_path = "simulations/HPC/sim_mvt/sim_HPC_out"
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)

tikzPlot = 0
texDir = paste0("tex/")

if (tikzPlot == 1) {
  if (!require("tikzDevice")) install.packages("tikzDevice")
  library(tikzDevice)
}

for (sim_setting in 1:3){

  n_simu = 50000
  n_array = 5000
  load(file=paste0(main_path, "/sim_setting_", sim_setting,"/results_mvt_sim_setting_", sim_setting, "_n_simu_", n_simu, "_n_array_", n_array, ".rdata"))

  sol_qtost_avg = sol_TOST_avg
  sol_alpha_qtost_avg = sol_piTOST_avg
  siz_qtost = sol_qtost_avg / 100
  siz_alpha_qtost = sol_alpha_qtost_avg / 100
  heatmap_range = range(siz_qtost, siz_alpha_qtost)

  rnam = rownames(siz_qtost)
  rnam = round(as.numeric(gsub("theta_2=", "", rnam)), 2)
  cnam = colnames(siz_qtost)
  cnam = round(as.numeric(gsub("theta_1=", "", cnam)), 2)
  rownames(siz_qtost) = rownames(siz_alpha_qtost) = rnam
  colnames(siz_qtost) = colnames(siz_alpha_qtost) = cnam

  alpha = 0.05
  B = n_simu
  delta = 0.15
  pi_x = c(0.2, 0.8)
  delta_l = pi_x - delta
  delta_u = pi_x + delta

  my_cex = 1.7
  if (tikzPlot){
    fil_nam = paste0(texDir, c("fig4", "figS7", "figS8")[sim_setting], ".tex")
    tikz(fil_nam, width = 10, height = 8, standAlone = TRUE,
         packages = c("\\usepackage{tikz}",
                      "\\usepackage[active,tightpage,psfixbb]{preview}",
                      "\\PreviewEnvironment{pgfpicture}",
                      "\\setlength\\PreviewBorder{0pt}",
                      "\\usepackage{amssymb}",
                      "\\usepackage{bm}","\\usepackage{amsthm}","\\usepackage{amsbsy}",
                      "\\usepackage{amsbsy}",
                      "\\usepackage{amsbsy}",
                      "\\usepackage{amsfonts}"))
  }

  source("simulations/aux_fun/aux_heat_map.R")
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow=2, byrow=T), widths = c(4, 1, 4, 1, 4, 1, 5))
  custom_heatmap(
    siz_qtost,
    my_cex=my_cex,
    xlab = "$\\theta_1$",
    ylab = "$\\theta_2$",
    main = "qTOST",
    legend = F,
    vline = qnorm(c(delta_l[1], delta_u[1])),
    hline = qnorm(c(delta_l[2], delta_u[2])),
    lambda=NULL,
    alpha=alpha,
    custom_mfrow=T,
    min_val=heatmap_range[1], max_val=heatmap_range[2]
  )
  plot.new()
  custom_heatmap(
    siz_alpha_qtost,
    my_cex=my_cex,
    xlab = "$\\theta_1$",
    ylab = "$\\theta_2$",
    main = "$\\bm{\\alpha}$-qTOST",
    legend = TRUE,
    legend_title = "Pr(reject H$_0$)",
    vline = qnorm(c(delta_l[1], delta_u[1])),
    hline = qnorm(c(delta_l[2], delta_u[2])),
    lambda=NULL,
    alpha=alpha,
    custom_mfrow=T,
    min_val=heatmap_range[1], max_val=heatmap_range[2]
  )

  custom_heatmap(
    siz_alpha_qtost-siz_qtost,
    my_cex=my_cex,
    xlab = "$\\theta_1$",
    ylab = "$\\theta_2$",
    main = "$\\bm{\\alpha}$-qTOST$-$qTOST",
    legend = TRUE,
    legend_title = "Diff. in \n Pr(reject H$_0$)\n",
    vline = qnorm(c(delta_l[1], delta_u[1])),
    hline = qnorm(c(delta_l[2], delta_u[2])),
    lambda=NULL,
    alpha=0,
    custom_mfrow=T
  )

  qtost_bound = siz_qtost[3:48, 3:48]
  qtost_bound = c(qtost_bound[1,], qtost_bound[nrow(qtost_bound),],
                 qtost_bound[,1], qtost_bound[,ncol(qtost_bound)])
  alpha_qtost_bound = siz_alpha_qtost[3:48, 3:48]
  alpha_qtost_bound = c(alpha_qtost_bound[1,], alpha_qtost_bound[nrow(alpha_qtost_bound),],
                   alpha_qtost_bound[,1], alpha_qtost_bound[,ncol(alpha_qtost_bound)])
  source("simulations/aux_fun/aux_hist.R")
  custom_hists(qtost_bound, alpha_qtost_bound, B, alpha, my_cex-0.2)

  if (tikzPlot){
    dev.off()
  }
}
