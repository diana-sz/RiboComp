library(here)

fluxes_bremer <- read.csv(here("data", "fluxes_bremer.csv"), row.names = 1)
fluxes1 <- read.csv(here("data", "05_fluxes_RNAPmax.csv"), row.names = 1)
fluxes2 <- read.csv(here("data", "06_fluxes_RNAPmax_noacc.csv"), row.names = 1)
fluxes3 <- read.csv(here("data", "07_fluxes_RNAPmax_noacc_act.csv"), row.names = 1)

fig_size <- c(18,13)
cex_labs <- 1.8
cex_axis <- 1.6

uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_green <- "#94c154"

plot_RNAP <- function(data, data2, xlim, cex_labs, cex_axis, col_fluxes){
  par(mar = c(4.5,5.5,0.5,0.5))
  ylim <- c(0.00000001, 0.0015)
  plot(vRNAP ~ mu, data, log = "y",
       ylim = ylim, xlim = xlim,
       col = col_fluxes, pch = 19, cex = 0.9,
       xlab = NA, ylab = NA, axes = FALSE)
  par(new = TRUE)
  plot(fluxes ~ mu, data2, log = "y", 
       ylim = ylim, xlim = xlim,
       col = uni_green, pch = 18, cex = 2,
       ylab = expression(paste("RNAP flux [mmol g"^"-1"* "h"^"-1"*"]")),
       xlab = expression(paste("Growth rate [h"^"-1"*"]")),
       cex.lab = cex_labs, cex.axis = cex_axis)
  par(new = FALSE)
}


png(filename = here("plots", "fluxes_RNAPmax.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
plot_RNAP(fluxes1, fluxes_bremer, c(0,3.5), cex_labs, cex_axis, uni_blue)
dev.off()


png(filename = here("plots", "fluxes_RNAPmax_noacc.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
plot_RNAP(fluxes2, fluxes_bremer, c(0,3.5), cex_labs, cex_axis, uni_blue)
dev.off()


png(filename = here("plots", "fluxes_RNAPmax_noacc_act.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
plot_RNAP(fluxes3, fluxes_bremer, c(0,3.5), cex_labs, cex_axis, uni_red)
dev.off()


