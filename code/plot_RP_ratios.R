library(here)

#### read simulation results ###################################################
all_data <- read.csv(here("data", "RBA_mass_fractions.csv"))

names <- c("base_activities",
           "base_activities_var_kel",
           "extended_activities", 
           "extended_hill-2_activities", 
           "extended_hill-6_activities", 
           "extended_hill-6_activities_var_kel")
media <- c("_LB", "_glcAA", "_glyAA", "_glc", "_gly", "_succ")
legend_media <- c("LB", "Glc+AA", "Gly+AA", "Glc", "Gly", "Succ")


#### define plotting parameters and functions ##################################
group_data <- function(data, columns){
  grouped <- aggregate(data[, columns], list(data$prot_fraction), FUN='mean')[, columns]
  return(grouped)
}

leg_size <- 1.4
ylim <- c(0,4.5)
xlim <- c(0,1)
mar <- c(4.5,5.5,0.5,0.5)
fig_size <- c(18,14)
uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_teal <- "#11897a"
uni_green <- "#94c154"
uni_yellow <- "#f6a800"
uni_orange <- "#dd4814"
colors <- c(uni_yellow, uni_green, uni_teal, uni_blue, "gray26", "gray60")


#### Plot RNA:protein ratio ####################################################

# get values at maximum growth rates for all conditions
x_opt <- 0.36000
for(name in names){
    opt_RP_ratio <- data.frame()

    for(dataset in c(paste0(name, media))){
      protein_cols <- c("R", "AF", "RNAP", "ENT", "EAA", "IG")
      if(grepl("extended", dataset)){
        protein_cols <- c(protein_cols, "RNase")}
  
      one <- all_data[all_data$name == dataset, ]
      opt <- group_data(one[one$prot_fraction == x_opt,], c(protein_cols, "growth_rate"))
      opt$rP <- opt$R * x_opt
      opt$RNA <- opt$R * (1-x_opt)
      opt$RP_ratio <- opt$RNA / rowSums(opt[, c(protein_cols[2:length(protein_cols)], "rP")])
      opt_RP_ratio <- rbind(opt_RP_ratio, opt)
    }

    fit <- lm(RP_ratio ~ growth_rate, opt_RP_ratio)
    
    plot_name <- gsub("_activities", "", name)
    png(filename = here("plots", paste0("ratio_RNA_protein_", plot_name, ".png")), 
        type="cairo", units="cm", 
        width=18, height=14, res=300)
    par(mar = mar)
    
    plot(NA, xlab = NA, ylab = NA,
         ylim = c(0,1), xlim = c(0,3.5), axes = FALSE)
    abline(fit, col = "grey69", lwd = 0.5)
    par(new = TRUE)
    plot(RP_ratio ~ growth_rate, opt_RP_ratio, 
         ylim = c(0,1), 
         xlim = c(0,3.5),
         ylab = "RNA/protein mass ratio",
         xlab = expression(paste("Growth rate [h"^"-1"*"]")),
         pch = 19,
         cex = 1.4,
         col = colors,
         cex.axis = 1.6, 
         cex.lab = 1.8)
    par(new = FALSE)
    
    if(name %in% c("base_activities",
                   "base_activities_var_kel",
                   "extended_activities", 
                   "extended_hill-6_activities_var_kel")){
      legend("topleft", 
             legend = legend_media,
             pch = 19,
             col = colors,
             bty = "n", 
             cex = leg_size)
    }
    dev.off()
}



