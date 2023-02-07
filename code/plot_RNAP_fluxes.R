library(here)

#### read simulation results and put them in one list ##########################
filenames <- c(#"fluxes_RNAPmax_noact",
               #"fluxes_RNAPmax_noacc_noact",
               #"fluxes_RNAPmax_noacc_glc",
               "R_deg",
               "R_deg2",
               "R_deg_hill")

all_fluxes <- list()
for(filename in filenames){
  data <- read.csv(here("data", paste0("fluxes_", filename, ".csv")))
  for(name in unique(data$name)){
    all_fluxes[[name]] <- data[data$name == name,]
  }
}


#### read experimental RNAP fluxes and correct them for degradation of rRNA ####
fluxes_bremer <- read.csv(here("data", "fluxes_bremer.csv"), row.names = 1)
fluxes_bremer$fluxes_corr <- fluxes_bremer$fluxes * c(1.3, 1.13, 1.1, 1.1, 1.1) 


#### define plotting parameters and functions ##################################
fig_size <- c(18,13)
uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_green <- "#94c154"
uni_teal <- "#11897a"

plot_RNAP <- function(data, data2, xlim, col_fluxes){
  par(mar = c(4.5,5.5,0.5,0.5))
  ylim <- c(0.00000001, 0.015)

  # simulated data
  plot(vRNAP ~ growth_rate, data, 
       log = "y",
       ylim = ylim, xlim = xlim,
       col = col_fluxes, 
       pch = 19, cex = 0.9,
       xlab = NA, ylab = NA, 
       axes = FALSE)
  par(new = TRUE)

  # experimental data
  plot(fluxes_corr ~ mu, data2, 
       log = "y", 
       ylim = ylim, xlim = xlim,
       col = uni_teal, 
       pch = 17, cex = 1.8,
       xlab = NA, ylab = NA,
       axes = FALSE)
  par(new = TRUE)
  plot(fluxes ~ mu, data2, 
       log = "y", 
       ylim = ylim, xlim = xlim,
       col = uni_green, 
       pch = 18, cex = 1.8,
       ylab = expression(paste("RNAP flux [mmol g"^"-1"* "h"^"-1"*"]")),
       xlab = expression(paste("Growth rate [h"^"-1"*"]")),
       cex.lab = 1.8, cex.axis = 1.6)
  par(new = FALSE)
}


#### plot data #################################################################
to_plot <- c("R_deg_glc", "R_deg2_glc", 
             "R_deg_hill_noact",
             "R_deg_hill_glc_noacc", "R_deg_hill_glc")
for(dataset in to_plot){
  png(filename=here("plots", paste0("fluxes_", dataset, ".png")), 
      type = "cairo", 
      units = "cm", 
      width = fig_size[1], 
      height = fig_size[2], 
      res = 300)
  colour <- ifelse(grepl("glc", dataset), uni_red, uni_blue)
  plot_RNAP(all_fluxes[[dataset]], 
            fluxes_bremer, 
            c(0,2.5),
            colour)
  dev.off()
}

