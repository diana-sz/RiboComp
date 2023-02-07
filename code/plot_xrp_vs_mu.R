library(here)

#### read simulation results and put them in one list ##########################
filenames <- c("RNAPmax", "R_deg_hill", "R_deg2", "R_deg", 
               "RBA", "RBA_init_rate", "Kostinski_reproduction")

all_data <- list()
for(filename in filenames){
  data <- read.csv(here("data", paste0(filename, "_growth_rates.csv")))
  for(name in unique(data$name)){
    all_data[[name]] <- data[data$name == name,]
  }
}


#### save max. mu and optimal xrP ##############################################
stats <- list()
for(name in names(all_data)){
  one <- all_data[[name]]
  stats[[name]] <- c(max(one$growth_rate), 
                     one$prot_fraction[which.max(one$growth_rate)])
}


#### define plotting parameters and functions ##################################
fig_size <- c(18,13)
leg_size <- 1.4
ylim <- c(0,3)
xlim <- c(0,1)
lwd <- 7
mar <- c(4.5,5.5,0.5,0.5)

uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_teal <- "#11897a"
uni_green <- "#94c154"
uni_yellow <- "#f6a800"

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
       lwd = 7)
}

add_max_mu_lines <- function(stats, dataset_names, colors){
  for(i in 1:length(dataset_names)){
    name <- dataset_names[i]
    lines(c(stats[[name]][2], stats[[name]][2]),
          c(-1, stats[[name]][1]),
          col = colors[i], lwd = 2)
  }
}

add_curve_labels <- function(stats, dataset_names, labels){
  for(i in 1:length(dataset_names)){
    name <- dataset_names[i]
    text(stats[[name]][2]+0.08, stats[[name]][1]+0.12,
         labels[i], cex = leg_size)
  }
}

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


#### RBA standard #############################################################
RBA <- all_data[["RBA_standard"]]
png(filename = here("plots","RBA_standard_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
par(mar = mar)
plot_curve(RBA, xlim, ylim, 3, "grey48")
dev.off()


#### RBA reverse ###############################################################
png(filename = here("plots", "RBA_reverse_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
par(mar = mar)

plot_curves(all_data, stats,
            c("RBA_standard", "RBA_reverse"),
            c("grey78", uni_teal),
            c(3,1), xlim, ylim, FALSE)

legend("bottomright", legend = c("RBA", "RBA (RNA expensive)"), 
       lty = c(1,3), lwd = lwd,
       col = c("grey78", uni_teal), bty = "n", cex = leg_size)
dev.off()


#### RBA degradation E. coli vs. Archaea (3 types of degradation function) #####
for(deg_type in c("RNAPmax", "R_deg_hill", "R_deg2", "R_deg")){
  png(filename = here("plots", paste0(deg_type, "_arch_mus.png")), type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  par(mar = mar)
  
    plot_curves(all_data, stats,
              c(paste0(deg_type, c("_noact", "_arch"))),
              c(uni_blue, uni_red), 
              c(1,1,1), xlim, ylim, TRUE)
    par(new=TRUE)
    plot_curve(all_data[["RBA_standard"]], xlim, ylim, lty=3, col="grey78")
    add_curve_labels(stats, 
                   paste0(deg_type, c("_noact", "_arch")),
                   c("E. coli", "Archaea"))
    par(new=FALSE)
  dev.off()
}



#### RBA degradation media #####################################################
colors <- c(uni_blue, uni_teal,uni_green, uni_yellow, uni_red, "grey28")

for(deg_type in c("R_deg_hill", "R_deg2", "R_deg", "RNAPmax", "PRL")){
  png(filename = here("plots", paste0(deg_type, "_media_mus.png")), type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  par(mar = c(mar[1:3], 0.5))
  
  ylim <- ifelse(deg_type == "PRL", 2.2, 3.5)
  
  plot_curves(all_data, stats,
              paste0(deg_type, c("_LB", "_glcAA", "_glyAA", "_glc", "_gly", "_succ")), 
              colors, c(1,1,1,1,1,1), xlim, c(0,ylim), TRUE)
  par(xpd = NA)
  
  if(deg_type %in% c("R_deg")){
    legend(0.48, 0.82, legend = c("LB", "Glc+AA", "Gly+AA"),
           lty = 1, lwd = lwd,
           col = c(uni_blue, uni_teal, uni_green),
           bty = "n", cex = leg_size)
    legend(0.78, 0.82, legend = c("Glc", "Gly", "Succ"),
           lty = 1, lwd = lwd,
           col = c(uni_yellow, uni_red, "grey28"),
           bty = "n", cex = leg_size)
  }
  
  if(deg_type %in% c("PRL")){
    legend("topright", legend = c("LB", "Glc+AA", "Gly+AA", "Glc", "Gly", "Succ"),
           lty = 1, lwd = lwd,
           col = colors,
           bty = "n", cex = leg_size)
  }
  
  par(xpd = FALSE)
  dev.off()
}



#### RBA vs RNA+RNAPmax ########################################################
# png(filename = here("plots", "RNAPmax_mus.png"), type="cairo", units="cm", 
#     width=fig_size[1], height=fig_size[2], res=300)
# par(mar = mar)
# 
# plot_curves(all_data, stats,
#             c("RNAPmax_noact", "RNAPmax_glc", "RBA_standard"), 
#             c(uni_blue, uni_red, "grey78"),
#             c(1, 1, 3), xlim, ylim, FALSE)
# 
# my.expressions <-c(as.expression(bquote("RBA")),
#                    as.expression(bquote("RBA+RNAP"["max"])),
#                    as.expression(bquote("RBA+RNAP"["max"]*"+activities")))
# legend("bottomright", legend = my.expressions, lty = c(3,1,1), lwd = lwd,
#        col = c("grey78", uni_blue, uni_red), bty = "n", cex = leg_size)
# 
# dev.off()


# 
# #### Archaea ###################################################################
# png(filename = here("plots", "RNAPmax_arch_mus.png"), type="cairo", units="cm", 
#     width=fig_size[1], height=fig_size[2], res=300)
# par(mar = mar)
# plot_curve(all_data[["RNAPmax_arch"]], xlim, c(0,2), 1, uni_blue)
# dev.off()


#### v_RNAP<=k_init*c_RNAP #####################################################
k_in_cap <- all_data[["RBA_init_rate"]]
k_ins <- unique(k_in_cap$kin)
cols <- c(uni_blue, uni_red, uni_teal, uni_green, uni_yellow)[1:length(k_ins)]

png(filename = here("plots", "RBA_init_rate_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)

par(mar = mar)
for(v in 1:length(k_ins)){
  plot_curve(k_in_cap[k_in_cap$kin == k_ins[v],], xlim, c(0,3.05), 1, cols[v])
  par(new=TRUE)
}
plot_curve(RBA, xlim, c(0,3.05), 3, "grey78")

par(new=FALSE)
legend(0.82, 3.14, legend = c("RBA", k_ins), 
       lty = c(3,rep(1,5)), lwd = lwd,
       col = c("grey78", cols), bty = "n", cex = 1.18, 
       title = as.expression(bquote("k"["RNAP"]^"in"*" [h"^"-1"*"]")))

dev.off()


#### Print max growth rates ####################################################
for(name in names(stats)){
  one <- stats[[name]]
  print(paste0("Max. growth rate ", name," (", one[1],") at ", one[2]*100,"% rP"))
}

