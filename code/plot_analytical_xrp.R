library(here)
library(viridis)
library(ggplot2)
library(RColorBrewer)

xrp <- 0.36 
tol <- 0.01

kelsr_ecoli <- 0.85*c(12, 16.83, 21) #, 20.17, 21, 22.25)
kelsrnap_ecoli <- 85*c(13.2, 14.4, 15.0)/100  # 18.8, 24.2, 31.0
#phis_rnap <- c(0.18, 0.28, 0.42, 0.52, 0.60, 0.65)
#kelsrnap_ecoli <- kelsrnap_ecoli*phis_rnap

kels <- data.frame(kelsr_ecoli, kelsrnap_ecoli)
data <- read.csv(here("data", "analytical_xrp_act.csv"))

data$xrp_opt <- data$xrp
data$xrp_opt[data$xrp > xrp+tol] <- 0
data$xrp_opt[(data$xrp <= xrp+tol) & (data$xrp >= xrp-tol)] <- 1
data$xrp_opt[data$xrp < xrp-tol] <- 0.2
data$xrp_opt[(data$xrp < 0) | (data$xrp > 1)] <- NA

#data$xrp[(data$xrp > xrp+0.005) | (data$xrp < xrp-0.005)] <- 50

#data <- data[order(data$kel_r),]

data[, c("kel_r", "kel_rnap")] <- data[, c("kel_r", "kel_rnap")]/3600

(v <- ggplot(data) +
    geom_tile(data, mapping=aes(kel_r, kel_rnap, fill = xrp_opt)) +
    #geom_path(data=data[!is.na(data$xrp_opt),], aes(kel_r, xrp_opt), color="white", na.rm=TRUE, size=2) +
    theme_minimal() +
    scale_fill_gradient(low = "white", high = "purple") +
  #scale_fill_distiller() +
    #scale_fill_viridis_c() +
    #scale_fill_viridis_c(option = "plasma") +
    #scale_fill_distiller(palette = "RdPu") +
    #scale_fill_viridis_c(option="viridis", name="xrP*") + 
    #geom_text(aes(x=kel_rs[i], y=kel_rnaps[i], label=labels[i], hjust=-0.2), color="white") +
    geom_point(data=kels, mapping=aes(x=kelsr_ecoli, y=kelsrnap_ecoli),shape=23, fill="red", color="white", size=3) +
    labs(x = expression(paste("Translation rate [AA h"^-1*"]")), y = "y axis") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
)
print(v)


data <- read.csv(here("data", "analytical_xrp_act.csv"))

data$xrp_opt <- data$xrp
data$xrp_opt[data$xrp > xrp+tol] <- 1
data$xrp_opt[(data$xrp <= xrp+tol) & (data$xrp >= xrp-tol)] <- 2
data$xrp_opt[data$xrp < xrp-tol] <- 3
data$xrp_opt[(data$xrp < 0) | (data$xrp > 1)] <- 4
data[, c("kel_r", "kel_rnap")] <- data[, c("kel_r", "kel_rnap")]/3600


par(mar = c(5,5,1,1))
colours <- data.frame(col = c(viridis(15, option = "magma")[c(2,15,4)], "grey90"),
                      #col = c(brewer.pal(11, "RdYlBu")[c(10,7,11)], "grey90"),
                      val = 1:4)
plot(kel_rnap ~ kel_r, data=data,
     xlim = c(0,25), ylim = c(0,25),
     cex = 0.4, pch = 15,
     col = colours[match(data$xrp_opt, colours$val), "col"],
     xlab = expression(paste("Translation rate [AA s"^-1*"]"*"k"["R"]^"el"*" f"["RNAP"]^"act")),
     ylab = expression(paste("Transcription rate [NT s"^-1*"]")),
     cex.lab = 1.3,
     cex.axis = 1.3)
par(new = TRUE)
plot(kelsr_ecoli, kelsrnap_ecoli, 
     xlim = c(0,25), ylim = c(0,25),
     cex = 2, pch = 23, col = "white", bg = "deeppink2",
     axes = FALSE, ylab = NA, xlab = NA)
par(new = FALSE)

