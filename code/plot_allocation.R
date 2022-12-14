library(here)
library(RColorBrewer)

phis1 <- read.csv(here("data", "RBA_standard_phis.csv"), row.names = 1)
phis2 <- read.csv(here("data", "RNAPmax_noact_phis.csv"), row.names = 1)
phis3 <- read.csv(here("data", "RNAPmax_act_phis.csv"), row.names = 1)
phis4 <- read.csv(here("data", "RNAPmax_arch_phis.csv"), row.names = 1)

all_data <- list("RBA_standard_phis" = phis1,
                 "RNAPmax_noact_phis" = phis2,
                 "RNAPmax_act_phis" = phis3,
                 "RNAPmax_arch_phis" = phis4)

metabolite_cols <- c("G", "AA", "NT", "rRNA", "rP")
enzyme_cols <- c("R", "AF", "RNAP", "ENT", "EAA", "IG")
syn_cols <- c("wrP", "wAF", "wRNAP", "wENT", "wEAA", "wIG")

fig_size <- c(18,13)
colours <-  c(brewer.pal(6, "Paired"), "grey96")

# find min. concentration/allocation across EGVs (=necessary minimum)
# everything above that is unnecessary accumulation
group_data <- function(data, columns){
  grouped <- aggregate(data[, columns], list(data$x), FUN='min')[, columns]
  storage <- 1-rowSums(grouped)
  if(max(storage) > 0.01){
    grouped$storage <- storage  # only include if it is big enough
  }
  return(grouped)
}

make_plot <- function(data, columns, main, colours){
  leg_size <- 1.35
  cex_labs <- 1.8
  cex_axis <- 1.6
  
  grouped <- group_data(data, columns)
  cumsums <- t(apply(grouped, 1, cumsum))  # cumulative sum to draw the polygons
  xs <- c(unique(data$x), rev(unique(data$x)))
  
  par(mar = c(4.5,5.5,1,1))
  plot(NA, 
       xlim = c(0,1), 
       ylim = c(0,1), 
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

#### Plot allocations ##########################################################
for(dataset in names(all_data)){
  png(filename = here("plots", paste0(dataset, ".png")), 
      type="cairo", units="cm", 
      width=fig_size[1], height=fig_size[2], res=300)
  make_plot(all_data[[dataset]], syn_cols, dataset, colours)
  dev.off()
}
