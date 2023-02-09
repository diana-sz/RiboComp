library(here)

data <- read.csv(here("data", "R_deg_hill_fluxes.csv"))
data <- data[data$name == "R_deg_hill_LB",]

mean_egv <- aggregate(data, list(data$prot_fraction), FUN='mean')
mean_egv$ratio <- mean_egv$vdeg/mean_egv$vRNAP

fig_size <- c(18,13)

png(filename=here("plots", "ratio_vRNAP_vdeg_hill_LB.png"), 
    type = "cairo", 
    units = "cm", 
    width = fig_size[1], 
    height = fig_size[2], 
    res = 300)
    plot(ratio ~ prot_fraction, mean_egv)
    abline(v=0.36)

dev.off()
