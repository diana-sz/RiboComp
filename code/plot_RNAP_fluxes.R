library(here)

fluxes_bremer <- read.csv(here("data", "fluxes_bremer.csv"), row.names = 1)
fluxes1 <- read.csv(here("data", "fluxes_RNAPmax_noact.csv"), row.names = 1)
fluxes2 <- read.csv(here("data", "fluxes_RNAPmax_noacc_noact.csv"), row.names = 1)
fluxes3 <- read.csv(here("data", "fluxes_RNAPmax_noacc_act.csv"), row.names = 1)

all_fluxes <- list("fluxes_RNAPmax" = fluxes1,
                   "fluxes_RNAPmax_noacc" = fluxes2,
                   "fluxes_RNAPmax_noacc_act" = fluxes3)

fig_size <- c(18,13)
uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_green <- "#94c154"

plot_RNAP <- function(data, data2, xlim, col_fluxes){
  par(mar = c(4.5,5.5,0.5,0.5))
  ylim <- c(0.00000001, 0.0015)
  plot(vRNAP ~ mu, 
       data, 
       log = "y",
       ylim = ylim, 
       xlim = xlim,
       col = col_fluxes, 
       pch = 19, 
       cex = 0.9,
       xlab = NA, 
       ylab = NA, 
       axes = FALSE)
  par(new = TRUE)
  plot(fluxes ~ mu, 
       data2, 
       log = "y", 
       ylim = ylim, 
       xlim = xlim,
       col = uni_green, 
       pch = 18, 
       cex = 2,
       ylab = expression(paste("RNAP flux [mmol g"^"-1"* "h"^"-1"*"]")),
       xlab = expression(paste("Growth rate [h"^"-1"*"]")),
       cex.lab = 1.8, 
       cex.axis = 1.6)
  par(new = FALSE)
}


for(dataset in names(all_fluxes)){
  png(filename=here("plots", paste0(dataset, ".png")), 
      type = "cairo", 
      units = "cm", 
      width = fig_size[1], 
      height = fig_size[2], 
      res = 300)
  colour <- ifelse(grepl("act", dataset), uni_red, uni_blue)
  plot_RNAP(all_fluxes[[dataset]], 
            fluxes_bremer, 
            c(0,3.5),
            colour)
  dev.off()
}



