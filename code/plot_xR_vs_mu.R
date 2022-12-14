library(here)

RBA <- read.csv(here("data", "RBA_standard_mus.csv"))
RBA_reverse <- read.csv(here("data", "RBA_reverse_mus.csv"))
RNAPmax <- read.csv(here("data", "RNAPmax_noact_mus.csv"))
RNAPmax_act <- read.csv(here("data", "RNAPmax_act_mus.csv"))
RNAPmax_arch <- read.csv(here("data", "RNAPmax_arch_mus.csv"))
k_in_cap <- read.csv(here("data", "init_rate_mus.csv"))

fig_size <- c(18,13)
leg_size <- 1.4
ylim <- c(0,4.5)
xlim <- c(0,1)
lwd <- 7
mar <- c(4.5,5.5,0.5,0.5)

uni_blue <- "#0063a6"
uni_red <- "#a71c49"
uni_teal <- "#11897a"
uni_green <- "#94c154"
uni_yellow <- "#f6a800"

plot_curve <- function(data, xlim, ylim, lty, col){
  plot(data$x, data$mu,
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


#----RBA only--------------------------------------------------------------#####
png(filename = here("plots","RBA_standard_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
par(mar = mar)

plot_curve(RBA, xlim, ylim, 3, "grey48")
dev.off()


#----RBA reverse-----------------------------------------------------------#####
png(filename = here("plots", "RBA_reverse_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
par(mar = mar)

plot_curve(RBA_reverse, xlim, ylim, 1, uni_teal)
par(new = TRUE)
plot_curve(RBA, xlim, ylim, 3, "grey78")
par(new = FALSE)

legend("bottomright", legend = c("RBA", "RBA (RNA expensive)"), 
       lty = c(1,3), lwd = lwd,
       col = c("grey78", uni_teal), bty = "n", cex = leg_size)
dev.off()


#----RBA vs RNA+RNAPmax----------------------------------------------------#####
png(filename = here("plots", "RNAPmax_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)
par(mar = mar)

plot_curve(RNAPmax, xlim, ylim, 1, uni_blue)
par(new = TRUE)
plot_curve(RNAPmax_act, xlim, ylim, 1, uni_red)
par(new = TRUE)
plot_curve(RBA, xlim, ylim, 3, "grey78")
par(new = FALSE)

my.expressions <-c(as.expression(bquote("RBA")),
                   as.expression(bquote("RBA+RNAP"["max"])),
                   as.expression(bquote("RBA+RNAP"["max"]*"+activities")))
legend("bottomright", legend = my.expressions, lty = c(3,1,1), lwd = lwd,
       col = c("grey78", uni_blue, uni_red), bty = "n", cex = leg_size)

dev.off()


#----Archaea---------------------------------------------------------------#####
png(filename = here("plots", "RNAPmax_arch_mus.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)

par(mar = mar)
plot_curve(RNAPmax_arch, xlim, c(0,2), 1, uni_blue)
dev.off()


#----v_RNAP<=k_init*c_RNAP-------------------------------------------------#####
k_ins <- unique(k_in_cap$kin)
cols <- c(uni_blue, uni_red, uni_teal, uni_green, uni_yellow)[1:length(k_ins)]

png(filename = here("plots", "RBA_init_rate.png"), type="cairo", units="cm", 
    width=fig_size[1], height=fig_size[2], res=300)

par(mar = mar)
for(v in 1:length(k_ins)){
  one_k_in <- k_in_cap[k_in_cap$kin == k_ins[v],]
  
  plot_curve(one_k_in, xlim, ylim, 1, cols[v])
  par(new=TRUE)
}
plot_curve(RBA, xlim, ylim, 3, "grey78")

par(new=FALSE)
legend(0.82, 4.61, legend = c("RBA", k_ins), 
       lty = c(3,rep(1,5)), lwd = lwd,
       col = c("grey78", cols), bty = "n", cex = 1.2, 
       title = as.expression(bquote("k"["RNAP"]^"in"*" [h"^"-1"*"]")))

dev.off()


#----Print max growth rates------------------------------------------------#####
print(paste0("Max. growth rate RBA (", max(RBA$mu),") at RP fraction: ", 
            RBA$x[which.max(RBA$mu)]))
print(paste0("Min. growth rate RBA (", min(RBA$mu),") at RP fraction: ", 
             RBA$x[which.min(RBA$mu)]))
print(paste0("Max. growth rate RBA+RNAPmax (", max(RNAPmax$mu),") at RP fraction: ", 
            RNAPmax$x[which.max(RNAPmax$mu)]))
print(paste0("Max. growth rate RBA+RNAPmax+activities (", max(RNAPmax_act$mu),") at RP fraction: ", 
            RNAPmax_act$x[which.max(RNAPmax_act$mu)]))
print(paste0("Maximum growth rate Archaea (", max(RNAPmax_arch$mu),") at RP fraction: ", 
            RNAPmax_arch$x[which.max(RNAPmax_arch$mu)]))

