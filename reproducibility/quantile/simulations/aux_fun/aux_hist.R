gg_color_hue = function(n, alpha = 1) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100, alpha = alpha)[1:n]
}

custom_hists = function(qtost_bound, alpha_qtost_bound, B, alpha=0.05, my_cex=1, ntst=200){
  
  nbrk = 20
  qtost_hist = hist(qtost_bound, breaks = nbrk, plot=F)
  alpha_qtost_hist = hist(alpha_qtost_bound, breaks = nbrk, plot=F)
  range_break = range(alpha_qtost_hist$breaks,max(0, qtost_hist$breaks))
  breaks_=seq(range_break[1],range_break[2],l=nbrk)
  qtost_hist = hist(qtost_bound, breaks = breaks_, plot=F)
  alpha_qtost_hist = hist(alpha_qtost_bound, breaks = breaks_, plot=F)
  
  alphacompl=1-alpha
  qnt = qnorm(1-0.01/ntst)
  upper = alpha+qnt*sqrt(alphacompl*alpha/B)
  lower = alpha-qnt*sqrt(alphacompl*alpha/B)
  
  mygrey = "grey50"
  grey_transp = adjustcolor("grey30", alpha.f = 0.15)
  cols = gg_color_hue(4)
  
  xmax = max(c(qtost_bound, alpha_qtost_bound, 0.06)) * 100
  plot(NA, xlim = c(-1.1, xmax+0.2), ylim = c(5.9, 10.4), axes = FALSE, xlab = "", ylab = "")
  max_height = 0.5*max(c(max(qtost_hist$density), max(alpha_qtost_hist$density)))
  K = length(qtost_hist$density)
  for (i in 1:K){
    rect(100*qtost_hist$breaks[i], 8, 100*qtost_hist$breaks[i+1], 8 + qtost_hist$density[i]/max_height, col = cols[1], border = mygrey)
    rect(100*alpha_qtost_hist$breaks[i], 6, 100*alpha_qtost_hist$breaks[i+1], 6 + alpha_qtost_hist$density[i]/max_height, col = cols[2], border = mygrey)
  }
  lines(c(-0.15,6), c(6,6))
  lines(c(-0.15,6), c(8,8))
  lines(c(5,5), c(5.5, 10), lwd = 1.5, lty = 2)
  rect(100*lower, 5.5, 100*upper, 10, col = grey_transp, border = NA)
  axis(1, at = 0:6, cex.axis = my_cex)
  mtext("Empirical Size (\\%)", side = 1, line = 2.65, adj = 0.6, cex = my_cex-0.3)
  text(-0.15, 6, "$\\alpha$-qTOST", pos = 2, cex = my_cex+0.1)
  text(-0.15, 8, "qTOST", pos = 2, cex = my_cex+0.1)
  text(5, 10.25, "99\\% Simulation Error\ntolerance under H$_0$ ($\\alpha = 5\\%$)", col = mygrey,
       cex= my_cex)
}

