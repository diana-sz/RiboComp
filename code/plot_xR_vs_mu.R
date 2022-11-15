library(here)

RBA <- read.csv(here("data", "00_RBA.csv"))
RBA_reverse <- read.csv(here("data", "01_RBA_reverse.csv"))
RNAPmax <- read.csv(here("data", "02_RNAPmax.csv"))
RNAPmax_act <- read.csv(here("data", "03_RNAPmax_act.csv"))
RNAPmax_arch <- read.csv(here("data", "04_RNAPmax_arch.csv"))
vmax <- read.csv(here("data", "09_vmax.csv"))

lwd <- 7
fig_size <- c(18,13)
leg_size <- 1.4
cex_labs <- 1.8
cex_axis <- 1.6
ylim <- c(0,4.5)
xlim <- c(0,1)

uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_teal <- "#11897a"
uni_green <- "#94c154"
uni_yellow <- "#f6a800"

xlab <- "Protein mass fraction of ribosome"
ylab <- expression(paste("Growth rate [h"^"-1"*"]"))


#----RBA only--------------------------------------------------------------#####
png(filename = here("plots","RBA.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)

par(mar = c(4.5,5.5,1,0.5))
plot(RBA$x, RBA$mu,
     ylim = ylim, xlim = xlim, type = "l", lty = 3,
     xlab = xlab, ylab = ylab,
     cex.axis = cex_axis, cex.lab = cex_labs, lwd = lwd, col = "grey48")
dev.off()


#----RBA reverse-----------------------------------------------------------#####
png(filename = here("plots", "RBA_reverse.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)

par(mar = c(4.5,5.5,1,0.5))
plot(RBA_reverse$x, RBA_reverse$mu,
     ylim = ylim, xlim = xlim, type = "l", lty = 1,
     xlab = xlab, ylab = ylab,
     cex.axis = cex_axis, cex.lab = cex_labs, 
     lwd = lwd,  col = uni_teal)

par(new = TRUE)
plot(RBA$x, RBA$mu,
     ylim = ylim, xlim = xlim, type = "l",  lty = 3,
     xlab = NA, ylab = NA, axes = FALSE,
     lwd = lwd, col = "grey78")

par(new = FALSE)
legend("bottomright", legend = c("RBA", "RBA (RNA expensive)"), lty = c(1,1), lwd = lwd,
       col = c("grey78", uni_teal), bty = "n", cex = leg_size)
dev.off()


#----RBA vs RNA+RNAPmax----------------------------------------------------#####
png(filename = here("plots", "RNAPmax.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
par(mar = c(4.5,5.5,1,0.5))

plot(RNAPmax$x, RNAPmax$mu,
     ylim = ylim, xlim = xlim, type = "l", lty = 1,
     xlab = xlab, ylab = ylab,
     cex.axis = cex_axis, cex.lab = cex_labs, 
     lwd = lwd, col = uni_blue)
par(new = TRUE)

plot(RNAPmax_act$x, RNAPmax_act$mu,
     ylim = ylim, xlim = xlim, type = "l", lty = 1,
     xlab = NA, ylab = NA, axes = FALSE,
     lwd = lwd, col = uni_red)
par(new = TRUE)

plot(RBA$x, RBA$mu,
     ylim = ylim, xlim = xlim, type = "l", lty = 3,
     xlab = NA, ylab = NA, axes = FALSE,
     lwd = lwd, col = "grey78")
par(new = FALSE)

my.expressions <-c(as.expression(bquote("RBA")),
                   as.expression(bquote("RBA+RNAP"["max"])),
                   as.expression(bquote("RBA+RNAP"["max"]*"+activities")))

legend("bottomright", legend = my.expressions, lty = c(3,1,1), lwd = lwd,
       col = c("grey78", uni_blue, uni_red), bty = "n", cex = leg_size)

dev.off()


#----Archaea---------------------------------------------------------------#####
png(filename = here("plots", "RNAPmax_arch.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)

par(mar = c(4.5,5.5,1,0.5))
plot(RNAPmax_arch$x, RNAPmax_arch$mu,
     ylim = c(0,2), xlim = xlim, type = "l",
     xlab = xlab, ylab = ylab,
     cex.axis = cex_axis, cex.lab = cex_labs, lwd = lwd, col = uni_blue)

dev.off()



#----Vmax------------------------------------------------------------------#####
cols <- c(uni_red, uni_blue, uni_teal, uni_green, uni_yellow)
vmaxes <- unique(vmax$vmax)

png(filename = here("plots", "RBA_vmax.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)

par(mar = c(4.5,5.5,1,0.5))
for(v in 1:length(vmaxes)){
  one_vmax <- vmax[vmax$vmax == vmaxes[v],]
  plot(one_vmax$x, one_vmax$mu,
       ylim = ylim, xlim = xlim, type = "l", lty = 1,
       xlab = xlab, ylab = ylab,
       cex.axis = cex_axis, cex.lab = cex_labs, lwd = lwd, 
       col = cols[v])
  par(new=TRUE)
}
par(new=FALSE)
legend("topright", legend = vmaxes, lty = rep(1,5), lwd = lwd,
       col = cols, bty = "n", cex = 1.2, 
       title = as.expression(bquote("V"["max"]*" [mmol g"^"-1"*"h"^"-1"*"]")))

dev.off()


#----Print max growth rates------------------------------------------------#####

print(paste0("Max. growth rate RBA (", max(RBA$mu),") at RP fraction: ", 
            RBA$x[which.max(RBA$mu)]))
print(paste0("Max. growth rate RBA+RNAPmax (", max(RNAPmax$mu),") at RP fraction: ", 
            RNAPmax$x[which.max(RNAPmax$mu)]))
print(paste0("Max. growth rate RBA+RNAPmax+activities (", max(RNAPmax_act$mu),") at RP fraction: ", 
            RNAPmax_act$x[which.max(RNAPmax_act$mu)]))
print(paste0("Maximum growth rate Archaea (", max(RNAPmax_arch$mu),") at RP fraction: ", 
            RNAPmax_arch$x[which.max(RNAPmax_arch$mu)]))

