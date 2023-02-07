library(here)

# E. coli translation and transcription rates * activities of R and RNAP
kelsr_ecoli <- 0.85*c(12, 16.83, 21) #, 20.17, 21, 22.25)
kelsrnap_ecoli <- 85*c(13.2, 14.4, 15.0)/100  # 18.8, 24.2, 31.0

# optimal xrPs (protein fraction in ribosome) calculated analytically 
# for varying translation/transcription rates (output of analytical_xrp.py)
data <- read.csv(here("data", "analytical_xrp.csv"), row.names = 1)
data[, c("kel_r", "kel_rnap")] <- data[, c("kel_r", "kel_rnap")]/3600

xrps <- c(0, 0.2, 0.6, 0.8, 1, 0.36)  # xrPs to plot
tol <- 0.02
fig_size <- c(18,13)
uni_red <- "#a71c49"
uni_green <- "#94c154"
axis_lim <- c(0,25)

png(filename=here("plots", "analytical_xrp.png"), 
    type = "cairo", 
    units = "cm", 
    width = fig_size[1], 
    height = fig_size[2], 
    res = 300)

  par(mar = c(4.5,5,0.5,0.5))
  
  plot(NA, xlab = NA, ylab = NA, xlim = axis_lim, ylim = axis_lim, axes = FALSE)
  par(new=TRUE)
  
  for(xrp in xrps){
    data1 <- data[(data$xrp < xrp+tol) & (data$xrp > xrp-tol), ]
    data1<- rbind(c(0, 0.0000001, xrp), data1)
  
    smoothingSpline <- smooth.spline(data1$kel_r, data1$kel_rnap,
                                     w = 1/data1$kel_rnap,
                                     spar = 0.5)
    # plot(kel_rnap ~ kel_r, data=data1,
    #      col = "grey80",
    #      axes = FALSE,
    #      xlim = axis_lim,
    #      ylim = axis_lim,
    #      cex = 0.4,
    #      pch = 15)

    lines(smoothingSpline,
          lwd = ifelse(xrp == 0.36, 4, 1),
          col = ifelse(xrp == 0.36, uni_red, "grey50"))
    par(new=TRUE)

    text(22.5, data1$kel_rnap[length(data1$kel_rnap)]+0.6,
         adj = 0,
         labels = xrp,
         cex = 1.8)
    par(new=TRUE)
  
  }

  plot(kelsr_ecoli, kelsrnap_ecoli,
       pch = 23,
       cex = 1.5,
       col = uni_green,
       bg = uni_green,
       xlim = axis_lim, 
       ylim = axis_lim,
       xlab = NA,
       ylab = NA,
       cex.axis = 1.6)

  text(kelsr_ecoli, kelsrnap_ecoli-1.25, 
       labels = c("Suc", "Gly", "Glc"),
       cex = 1.6,
       adj = 0.5)

  title(#xlab = expression(paste("k"["R"]^"el"*" f"["R"]^"act"*" [AA s"^-1*"]")),
    xlab = expression(paste("Translation rate [AA s"^-1*"]")),
    line = 3.5,
    cex.lab = 1.8)
  title(#ylab = expression(paste("k"["RNAP"]^"el"*" f"["RNAP"]^"act"*" [NT s"^-1*"]")),
    ylab = expression(paste("Transcription rate [NT s"^-1*"]")),
    line = 2.5,
    cex.lab = 1.8)

dev.off()

