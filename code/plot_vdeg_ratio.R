library(here)

name <- "deg_hill-6"
data <- read.csv(here("data", paste0("RBA_deg_fluxes.csv")))
data <- data[grepl(name, data$name),]


ratio_opt <- data.frame()
par(mfrow = c(2, 3))
for(condition in paste0(name, c("_LB", "_glcAA", "_glyAA", "_glc", "_gly", "_succ"))){
  # only plot the datasets that include R and RNAP activities
  if(grepl("arch", condition) | grepl("ecoli", condition) | grepl("mito", condition)){
    next
  }
  
  # calculate mean of the EGVs at each growth rate
  # at optimum they are almost the same, but we never find the precise optimum
  data_one <- data[data$name == condition,]
  mean_egv <- aggregate(data_one, 
                        list("prot_fraction" = data_one$prot_fraction), 
                        FUN='mean')
  mean_egv$ratio <- mean_egv$vRNase/mean_egv$vRNAP
  
  #plot(ratio ~ prot_fraction, mean_egv, ylim = c(0,1))
  #abline(v=0.36)

  ratio_opt <- rbind(ratio_opt, mean_egv[mean_egv$prot_fraction == 0.36000, c("growth_rate", "ratio")])
}
par(mfrow = c(1, 1))


fig_size <- c(18,13)
mar <- c(4.5,5.5,0.5,0.5)
png(filename=here("plots", paste0("ratio_RNAP_RNase_", name, ".png")),
    type = "cairo",
    units = "cm",
    width = fig_size[1],
    height = fig_size[2],
    res = 300)
  par(mar = mar)
  plot(ratio ~ growth_rate, ratio_opt,
       xlim = c(0,3.2), ylim = c(0,1),
       xlab = expression(paste("Growth rate [h"^"-1"*"]")),
       ylab = "RNAP flux / R degradation flux",
       cex.axis = 1.6, 
       cex.lab = 1.8, 
       pch = 19,
       cex = 1.2
  )
dev.off()

