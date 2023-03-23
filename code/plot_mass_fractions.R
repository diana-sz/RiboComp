library(here)
library(RColorBrewer)

#### read simulation results and put them in one list ##########################
all_data <- list()

filenames <- c("RBA_deg")#,
for(filename in filenames){
  data <- read.csv(here("data", paste0(filename, "_mass_fractions.csv")))
  for(name in unique(data$name)){
    all_data[[name]] <- data[data$name == name,]
  }
}


#### define plotting parameters and functions ##################################
fig_size <- c(18,13)
colours <-  c(brewer.pal(7, "Paired"), "grey96")

group_data <- function(data, columns){
  grouped <- aggregate(data[, columns], list(data$prot_fraction), FUN='mean')[, columns]
  return(grouped)
}


make_plot <- function(data, columns, main, colours){
  leg_size <- 1.35
  cex_labs <- 1.8
  cex_axis <- 1.6
  
  # group data from different EGVs (at max. growth rate they should be almost the same)
  grouped <- group_data(data, columns)
  
  # cumulative sum to draw the polygons
  cumsums <- t(apply(grouped, 1, cumsum))  
  xs <- c(unique(data$prot_fraction), rev(unique(data$prot_fraction)))
  
  par(mar = c(4.5,5.5,1,1))
  plot(NA, 
       xlim = c(0,1), 
       ylim = c(0,1), 
       xaxs = "i", yaxs = "i",
       xlab = "Protein mass fraction in ribosome",
       ylab = as.expression(bquote("Mass fraction per cell")),
       cex.lab = cex_labs,
       cex.axis = cex_axis)
  for(protein in ncol(cumsums):1){
    ys <- c(cumsums[, protein], rep(0, nrow(cumsums)))
    polygon(xs, ys, 
            col = colours[protein], 
            border = NA)
    par(new=TRUE)
    box()
  }
  par(new=FALSE)
  legend(1,0, xjust = 1, yjust = 0,
         legend =  gsub("w", "", rev(colnames(grouped))),
         col = "black", pt.bg = rev(colours[1:ncol(grouped)]),
         pch = 22,  pt.cex = 2,
         cex = leg_size, bty = "n")
}

#for(dataset in names(all_data)){
for(dataset in c("deg_glc")){  
  if(grepl("deg", dataset)){
    protein_cols <- c("R", "AF", "RNAP", "ENT", "EAA", "IG", "RNase")
  }else{
    protein_cols <- c("R", "AF", "RNAP", "ENT", "EAA", "IG")
  }
  # png(filename = here("plots", paste0(dataset, "_mass_fraction.png")), 
  #     type="cairo", units="cm", 
  #     width=fig_size[1], height=fig_size[2], res=300)
  make_plot(all_data[[dataset]], protein_cols, dataset, colours)
  #dev.off()
}


leg_size <- 1.4
ylim <- c(0,4.5)
xlim <- c(0,1)
mar <- c(4.5,5.5,0.5,0.5)


#### Plot RNA:protein ratio ####################################################
names <- c("deg", "deg_hill-2", "deg_hill-6")
for(name in names){

png(filename = here("plots", paste0("ratio_RNA_protein_", name, ".png")), type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  par(mar = mar)
  for(dataset in paste0(name, c("_LB", "_glcAA", "_glyAA", "_glc", "_gly", "_succ"))){
    
    if(grepl("deg", dataset)){
      protein_cols <- c("R", "AF", "RNAP", "ENT", "EAA", "IG", "RNase")
    }else{
      protein_cols <- c("R", "AF", "RNAP", "ENT", "EAA", "IG")
    }
    
    one <- all_data[[dataset]]
    x_opt <- 0.35500 #0.36000
    opt <- group_data(one[one$prot_fraction == x_opt,], c(protein_cols, "growth_rate"))
    opt$rP <- opt$R * x_opt
    opt$RNA <- opt$R * (1-x_opt)
    
    opt$RP_ratio <- opt$RNA / rowSums(opt[, c(protein_cols[2:length(protein_cols)], "rP")])
    
    plot(RP_ratio ~ growth_rate, opt, ylim = c(0,1), xlim = c(0,4),
         ylab = "RNA/protein mass ratio",
         xlab = expression(paste("Growth rate [h"^"-1"*"]")),
         pch = 19,
         cex = 1.2,
         cex.axis = 1.6, 
         cex.lab = 1.8)
    
    par(new=TRUE)
  }
  par(new=FALSE)
  dev.off()
}



