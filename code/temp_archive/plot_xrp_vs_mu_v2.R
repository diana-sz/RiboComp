library(here)

#### read simulation results and put them in one list ##########################
filenames <- c("RBA_deg")
all_data <- list()
for(filename in filenames){
  data <- read.csv(here("data", paste0(filename, "_growth_rates_v2_quick.csv")))
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
           one$growth_rate[which(one$prot_fraction == 0.355)])
 stats <- rbind(stats, row)
}
colnames(stats) <- c("xrp_opt", "mu_opt", "mu_0.36")
rownames(stats) <- names(all_data)


#### define plotting parameters and functions ##################################
fig_size <- c(18,13)
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

add_curve_labels <- function(stats, dataset_names, labels){
  for(i in 1:length(dataset_names)){
    name <- dataset_names[i]
    text(stats[name, "xrp_opt"]+0.08, stats[name, "mu_opt"]+0.12,
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





#### Six different media #######################################################
colors <- c(uni_blue, uni_teal,uni_green, uni_yellow, uni_red, "grey28")

for(deg_type in c("deg_hill-2", "deg_hill-6", "deg", "Kostinski"#,   "RNAPmax", 
                  )){
  png(filename = here("plots", paste0(deg_type, "_media_mus_v2.png")), type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  par(mar = c(mar[1:3], 0.5))
  
  ylim <- ifelse(deg_type == "Kostinski", 2.2, 3.5)
  
  plot_curves(all_data, stats,
              paste0(deg_type, c("_LB", "_glcAA", "_glyAA", "_glc", "_gly", "_succ")), 
              colors, c(1,1,1,1,1,1), xlim, c(0,ylim), TRUE)
  par(xpd = NA)
  
  if(deg_type %in% c("deg")){
    legend(0.48, 0.82, legend = c("LB", "Glc+AA", "Gly+AA"),
           lty = 1, lwd = lwd,
           col = c(uni_blue, uni_teal, uni_green),
           bty = "n", cex = leg_size)
    legend(0.78, 0.82, legend = c("Glc", "Gly", "Succ"),
           lty = 1, lwd = lwd,
           col = c(uni_yellow, uni_red, "grey28"),
           bty = "n", cex = leg_size)
  }
  
  if(deg_type %in% c("Kostinski")){
    legend("topright", 
           legend = c("LB", "Glc+AA", "Gly+AA", "Glc", "Gly", "Succ"),
           lty = 1, lwd = lwd,
           col = colors,
           bty = "n", cex = leg_size)
  }
  par(xpd = FALSE)
  dev.off()
}




