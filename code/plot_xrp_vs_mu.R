library(here)

#### read simulation results and put them in one list ##########################
filenames <- c("RBA_deg", "Kostinski_reproduction")
all_data <- list()
for(filename in filenames){
  data <- read.csv(here("data", paste0(filename, "_growth_rates.csv")))
  data$prot_fraction <- round(as.numeric(data$prot_fraction),3)
  for(name in unique(data$name)){
    all_data[[name]] <- data[data$name == name,]
  }
}

#### save max. mu and optimal xrP ##############################################
stats <- data.frame()
for(name in names(all_data)){
  one <- all_data[[name]]
  row <- c(one$prot_fraction[which.max(one$growth_rate)],
           max(one$growth_rate), 
           one$growth_rate[which(one$prot_fraction == 0.6)])
 stats <- rbind(stats, row)
}
colnames(stats) <- c("xrp_opt", "mu_opt", "mu_0.36")
rownames(stats) <- names(all_data)

write.csv(stats, here("data", "stats.csv"))

#### define plotting parameters and functions ##################################
fig_size <- c(18,14)
leg_size <- 1.4
ylim <- c(0,3)
xlim <- c(0,1)
lwd <- 6
mar <- c(4.5,5.5,0.5,0.5)

uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_teal <- "#11897a"
uni_green <- "#94c154"
uni_yellow <- "#f6a800"
uni_orange <- "#dd4814"

plot_curve <- function(data, xlim, ylim, lty, col){
  plot(data$prot_fraction, data$growth_rate,
       xlim = xlim,
       ylim = ylim,  
       lty = lty,
       col = col,
       type = "l", 
       xlab = "Protein mass fraction of ribosome", 
       ylab = expression(paste("Growth rate [h"^"-1"*"]")),
       cex.axis = 1.6, 
       cex.lab = 1.8, 
       lwd = 6)
}

add_max_mu_lines <- function(stats, dataset_names, colors){
  for(i in 1:length(dataset_names)){
    name <- dataset_names[i]
    lines(c(stats[name, "xrp_opt"], stats[name, "xrp_opt"]),
          c(-1, stats[name, "mu_opt"]),
          col = colors[i], lwd = 2)
  }
}

# add_curve_labels <- function(stats, dataset_names, labels){
#   for(i in 1:length(dataset_names)){
#     name <- dataset_names[i]
#     text(stats[name, "xrp_opt"]+0.08, stats[name, "mu_opt"]+0.12,
#          labels[i], cex = leg_size)
#   }
# }

plot_curves <- function(all_data, stats, dataset_names, colors, ltys, xlim, ylim, max_mu_lines){
  plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ylab = NA, xlab = NA)
  
  if(max_mu_lines == TRUE){
    add_max_mu_lines(stats, dataset_names, colors)
  }
  par(new=TRUE)
  
  for(i in 1:length(dataset_names)){
    dataset <- all_data[[dataset_names[i]]]
    plot_curve(dataset, xlim, ylim, ltys[i], colors[i])
    par(new=TRUE)
  }
  par(new=FALSE)
}


media <- c("_LB", "_glcAA", "_glyAA", "_glc", "_gly", "_succ")
colors <- c(uni_yellow, uni_green, uni_teal, uni_blue, "gray26", "gray60")
legend_media <- c("LB", "Glc+AA", "Gly+AA", "Glc", "Gly", "Succ")

#### RBA standard #############################################################
# RBA <- all_data[["RBA_standard"]]
# png(filename = here("plots","RBA_standard_mus.png"), type="cairo", units="cm", 
#     width=fig_size[1], height=fig_size[2], res=300)
# par(mar = mar)
# plot_curve(RBA, xlim, ylim, 3, "grey48")
# 
# dev.off()


#### RBA reverse ###############################################################
png(filename = here("plots", "RBA_reverse_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
  par(mar = mar)
  
  plot_curves(all_data, stats,
              c("RBA_glc", "RBA_reverse"),
              c("grey78", uni_teal),
              c(3,1), xlim, ylim, FALSE)
  
  legend("bottomright", legend = c("RBA", "RBA (RNA expensive)"), 
         lty = c(1,3), lwd = lwd,
         col = c("grey78", uni_teal), bty = "n", cex = leg_size)
dev.off()


#### E. coli vs. Archaea #######################################################
for(deg_type in c("deg_hill-2", "deg_hill-6", "deg")){
  png(filename = here("plots", paste0(deg_type, "_arch_mus.png")), type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  par(mar = mar)
  
  plot_curves(all_data, stats,
            c(paste0(deg_type, c("_glc", "_arch3", "_arch2"))),
            c(uni_blue, uni_red, uni_red), 
            c(1,3,1), xlim, ylim, TRUE)
  #par(new=TRUE)
  #plot_curve(all_data[["RBA_glc"]], xlim, ylim, lty=3, col="grey78")
  if (deg_type == "deg"){
    legend("topleft", 
           legend = c("E. coli",
                      expression(paste("Archaea - higher k"["deg"])),
                      expression(paste("Archaea - higher k"["deg"]," + Thermococcus parameters"))),
                      lty = c(1,1,3), 
           col = c(uni_blue, uni_red, uni_red), 
           lwd = lwd,
           bty = "n")
  }

    par(new=FALSE)
  dev.off()
}


#### E. coli vs. mitochondria ##################################################
for(deg_type in c("deg_hill-2", "deg_hill-6", "deg")){
  if (deg_type == "deg"){
    mt <- paste0(deg_type, "_mito")
  }else{
    mt <- gsub("hill", "hill_mito", deg_type)
  }
  to_plot <- c(paste0(deg_type, "_glc"), paste0(mt, "_mito"))
  
  png(filename = here("plots", paste0(deg_type, "_mito_mus.png")), type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  par(mar = mar)
  
  plot_curves(all_data, stats, to_plot,
              c(uni_blue, uni_red), 
              c(1,1,1), xlim, ylim, TRUE)
  par(new=TRUE)
  if (deg_type == "deg"){
    legend("topleft", 
           legend = c("E. coli", "Mitochondria"),
           lty = c(1,1,3), 
           col = c(uni_blue, uni_red), 
           lwd = lwd,
           bty = "n")
  }
  par(new=FALSE)
  dev.off()
}


#### Six different media #######################################################
png(filename = here("plots", "RBA_media_mus.png"), 
    type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
  par(mar = c(mar[1:3], 0.5))
  
  ylim <- 4
  plot_curves(all_data, stats,  paste0("RBA", media), 
              colors, c(1,1,1,1,1,1), xlim, c(0,ylim), FALSE)
  par(xpd = NA)
  
  legend(0.48, 0.94, legend = legend_media[1:3],
         lty = 1, lwd = lwd,
         col = colors[1:3],
         bty = "n", cex = leg_size)
  legend(0.78, 0.94, legend = legend_media[4:6],
         lty = 1, lwd = lwd,
         col = colors[4:6],
         bty = "n", cex = leg_size)
  
  par(xpd = FALSE)
dev.off()


#### Six different media #######################################################

for(deg_type in c("deg_hill-2", "deg_hill-6", "deg", "Kostinski")){
  png(filename = here("plots", paste0(deg_type, "_media_mus.png")), 
      type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  par(mar = c(mar[1:3], 0.5))
  
  ylim <- ifelse(deg_type == "Kostinski", 2.2, 3.5)
  
  plot_curves(all_data, stats,
              paste0(deg_type, media), 
              colors, c(1,1,1,1,1,1), xlim, c(0,ylim), TRUE)
  par(xpd = NA)
  
  if(deg_type %in% c("deg")){
    legend(0.48, 0.82, legend = legend_media[1:3],
           lty = 1, lwd = lwd,
           col = colors[1:3],
           bty = "n", cex = leg_size)
    legend(0.78, 0.82, legend =  legend_media[4:6],
           lty = 1, lwd = lwd,
           col = colors[4:6],
           bty = "n", cex = leg_size)
  }
  
  if(deg_type %in% c("Kostinski")){
    legend("topright", 
           legend = legend_media,
           lty = 1, lwd = lwd,
           col = colors,
           bty = "n", cex = leg_size)
  }
  par(xpd = FALSE)
  dev.off()
}


# #### Six different media #######################################################
# for(deg_type in c("deg_hill-6")){
#   png(filename = here("plots", paste0(deg_type, "_media2_mus.png")), type="cairo", units="cm", 
#       width=fig_size[1], height=fig_size[2], res=300)
#   par(mar = c(mar[1:3], 0.5))
#   
#   ylim <- 3.5
#   
#   plot_curves(all_data, stats,
#               paste0(deg_type, c("_LB2", "_glcAA2", "_glyAA2", "_glc2", "_gly2", "_succ2")), 
#               colors, c(1,1,1,1,1,1), xlim, c(0,ylim), TRUE)
#   #par(xpd = NA)
# 
#   par(xpd = FALSE)
#   dev.off()
# }
# 

# #### v_RNAP<=k_init*c_RNAP #####################################################
# k_in_cap <- all_data[["RBA_init_rate"]]
# k_ins <- unique(k_in_cap$kin)
# cols <- c(uni_blue, uni_red, uni_teal, uni_green, uni_yellow)[1:length(k_ins)]
# 
# png(filename = here("plots", "RBA_init_rate_mus.png"), type="cairo", units="cm", 
#     width=fig_size[1], height=fig_size[2], res=300)
# 
# par(mar = mar)
# for(v in 1:length(k_ins)){
#   plot_curve(k_in_cap[k_in_cap$kin == k_ins[v],], xlim, c(0,3.05), 1, cols[v])
#   par(new=TRUE)
# }
# plot_curve(RBA, xlim, c(0,3.05), 3, "grey78")
# 
# par(new=FALSE)
# legend(0.82, 3.14, legend = c("RBA", k_ins), 
#        lty = c(3,rep(1,5)), lwd = lwd,
#        col = c("grey78", cols), bty = "n", cex = 1.18, 
#        title = as.expression(bquote("k"["RNAP"]^"in"*" [h"^"-1"*"]")))
# 
# dev.off()


