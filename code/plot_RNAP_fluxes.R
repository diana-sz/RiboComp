library(here)

#### read simulated and experimental fluxes ####################################
data <- read.csv(here("data", paste0("fluxes_x0.36.csv")))
fluxes_bremer <- read.csv(here("data", "fluxes_bremer.csv"))

#### define plotting parameters and functions ##################################
fig_size <- c(18,14)
uni_blue <- "#0063a6"
uni_green <- "#94c154"
uni_teal <- "#11897a"

plot_RNAP <- function(data, data2, xlim, colour, point_size){
  par(mar = c(4.5,5.5,0.5,0.5))
  ylim <- c(0.00000001, 0.02)

  # simulated data
  plot(vRNAP ~ growth_rate, data, 
       log = "y",
       ylim = ylim, xlim = xlim,
       col = colour, 
       pch = 19, cex = point_size,
       xlab = NA, ylab = NA, 
       axes = FALSE)
  par(new = TRUE)

  # experimental data
  plot(fluxes_corrected ~ mu, data2,
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


#### plot RNAP fluxes ##########################################################
to_plot <- c("RBA_glc",
             "deg_glc",  
             "deg_hill-2_glc",
             "deg_hill-6_glc")
sim_types <- c("", "_noacc")

for(dataset in to_plot){
  png(filename=here("plots", paste0("fluxes_", dataset, ".png")), 
      type = "cairo", 
      units = "cm", 
      width = fig_size[1], 
      height = fig_size[2], 
      res = 300)
  
  for(sim_type in sim_types){
    colour <- ifelse(grepl("noacc", sim_type), uni_blue, "grey60")
    point_size <- ifelse(grepl("noacc", sim_type), 0.9, 0.5)
    dataset_name <- paste0(dataset, sim_type)
    plot_RNAP(data[data$name == dataset_name,], fluxes_bremer, 
              c(0,2.5), colour, point_size)
    par(new = TRUE)
  }
  par(new=FALSE)
  dev.off()
}

#### export data for one growth rate to examine the EGVs
one_mu <- data[data$growth_rate == unique(data$growth_rate)[10], ]
write.csv(one_mu, here("data", paste0("fluxes_one_mu.csv")))


