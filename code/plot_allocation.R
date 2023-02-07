library(here)
library(RColorBrewer)

#### read simulation results and put them in one list ##########################
filenames <- c("R_deg_hill", "RBA")
for(filename in filenames){
  data <- read.csv(here("data", paste0(filename, "_allocations.csv")))
  for(name in unique(data$name)){
    all_data[[name]] <- data[data$name == name,]
  }
}


#### define plotting parameters and functions ##################################
fig_size <- c(18,13)

# find min. concentration/allocation across EGVs (=necessary minimum)
# everything above that is unnecessary accumulation
group_data <- function(data, columns){
  grouped <- aggregate(data[, columns], list(data$prot_fraction), FUN='min')[, columns]
  
  # calculate the rest (inactive ribosomes)
  storage <- 1-rowSums(grouped)
  
  # only include storage fraction if it is big enough to plot
  if(max(storage) > 0.01){
    grouped$storage <- storage  
  }
  return(grouped)
}

make_plot <- function(data, columns, main, colours){
  leg_size <- 1.35
  cex_labs <- 1.8
  cex_axis <- 1.6
  
  grouped <- group_data(data, columns)
  cumsums <- t(apply(grouped, 1, cumsum))  # cumulative sum to draw the polygons
  xs <- c(unique(data$prot_fraction), rev(unique(data$prot_fraction)))
  
  par(mar = c(4.5,5.5,1,1))
  plot(NA, 
       xlim = c(0,1), ylim = c(0,1), 
       xaxs = "i", yaxs = "i",
       xlab = "Protein mass fraction in ribosome",
       ylab = as.expression(bquote("Ribosome allocation ["*Phi["R"]^"i"*"]")),
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


#### plot data #################################################################
for(dataset in c("R_deg_hill_glc", "RBA_standard")){ #names(all_data)
  
  syn_cols <- c("wrP", "wAF", "wRNAP", "wENT", "wEAA", "wIG")
  if(grepl("R_deg", dataset)){
    syn_cols <- c(syn_cols, "wRNase")
  }
  colours <-c(brewer.pal(length(syn_cols), "Paired"), "grey96")

  png(filename = here("plots", paste0(dataset, "_allocations.png")), 
      type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  make_plot(all_data[[dataset]], syn_cols, dataset, colours)
  dev.off()
}


