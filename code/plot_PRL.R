library(here)

data <- read.csv(here("data", "PRL_reproduction.csv"), row.names = 1)
colnames(data) <- c("Succinate", "Glycerol", "Glucose", "Glycerol+AA", "Glucose+AA", "LB")

colors <- list("Succinate" = "mediumpurple", 
               "Glycerol" = "orangered2" , 
               "Glucose" = "yellowgreen", 
               "Glycerol+AA" = "goldenrod1", 
               "Glucose+AA" = "chocolate3", 
               "LB" = "gray42")
lwd <- 5
fig_size <- c(18,13)
leg_size <- 1.4
cex_labs <- 1.8
cex_axis <- 1.6

png(filename = here("plots", "PRL_reproduction.png"), 
    type="cairo", 
    units="cm", 
    width=fig_size[1], 
    height=fig_size[2], 
    res=300)

par(mar = c(4.5,5.5,0.5,0.5))
for(condition in colnames(data)){
  plot(rownames(data), 
       data[,condition],
       ylim = c(0,2.2), 
       xlim = c(0,1), 
       type = "l",
       xlab = "Protein mass fraction of ribosome",
       ylab = expression(paste("Growth rate [h"^"-1"*"]")),
       cex.axis = cex_axis, 
       cex.lab = cex_labs, 
       lwd = lwd, 
       col = colors[[condition]])
  par(new=TRUE)
}
par(new=FALSE)

legend("topright", 
       legend = rev(colnames(data)), 
       bty = "n", 
       cex = leg_size, 
       col = rev(unlist(colors)), 
       lty = 1, 
       lwd = lwd)
abline(v=0.36, 
       lty = 2, 
       lwd = lwd)

dev.off()
