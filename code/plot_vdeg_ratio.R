library(here)

all_data <- read.csv(here("data", "RBA_deg_fluxes.csv"))
fraction_degraded_gausing <- read.csv(here("data", "gausing_RNA_deg.csv"))

names <- c("deg_hill-6", "deg_hill-2", "deg")
media <- c("_LB", "_glcAA", "_glyAA", "_glc", "_gly", "_succ")
legend_media <- c("LB", "Glc+AA", "Gly+AA", "Glc", "Gly", "Succ")


#### define plotting parameters and functions ##################################
fig_size <- c(18,14)
leg_size <- 1.4
mar <- c(4.5,5.5,0.5,0.5)
uni_blue <- "#0063a6"
uni_teal <- "#11897a"
uni_green <- "#94c154"
uni_yellow <- "#f6a800"
colors <- c(uni_yellow, uni_green, uni_teal, uni_blue, "gray26", "gray60")
exp_color <- "grey15"

#### Plot degradation:RNAP flux ratio ###########################################

for(name in names){
  data <- all_data[grepl(name, all_data$name),]
  ratio_opt <- data.frame()

  for(i in 1:length(media)){
    # calculate mean of the EGVs at each growth rate
    # at optimum they are almost the same, but we never find the precise optimum
    data_one <- data[data$name == paste0(name, media[i]),]
    mean_egv <- aggregate(data_one, 
                          list("prot_fraction" = data_one$prot_fraction), 
                          FUN='mean')
    mean_egv$ratio <- mean_egv$vRNase/mean_egv$vRNAP
    ratio_opt <- rbind(ratio_opt, mean_egv[mean_egv$prot_fraction == 0.36000, 
                                           c("growth_rate", "ratio")])
  }

  png(filename=here("plots", paste0("ratio_RNAP_RNase_", name, ".png")),
      type = "cairo",
      units = "cm",
      width = fig_size[1],
      height = fig_size[2],
      res = 300)
  par(mar = mar)
  plot(ratio ~ growth_rate, ratio_opt,
       xlim = c(0,3.2), 
       ylim = c(0,1),
       type = "l",
       col = "grey69",
       xlab = expression(paste("Growth rate [h"^"-1"*"]")),
       ylab = "RNAse flux / RNAP flux",
       cex.axis = 1.6, 
       cex.lab = 1.8)
  points(ratio_opt$growth_rate, ratio_opt$ratio,
         col = colors,
         pch = 19,
         cex = 1.4)
  par(new = TRUE)
  plot(percent_degraded/100 ~ growth_rate, fraction_degraded_gausing,
       type = "l",
       col = "grey69",
       xlim = c(0,3.2), ylim = c(0,1), 
       xlab = NA, ylab = NA, axes = FALSE)
  points(fraction_degraded_gausing$growth_rate,
         fraction_degraded_gausing$percent_degraded/100,
         pch = 17, cex = 1.4,
         col = exp_color)
  par(new = FALSE)

  if(name %in% c("deg")){
    legend("bottomright", 
           legend = c(legend_media, "Experimental"),
           pch = c(rep(19, 6), 17),
           col = c(colors, exp_color),
           bty = "n", 
           cex = leg_size)
  }
  dev.off()
}


