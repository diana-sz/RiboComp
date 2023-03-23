library(here)

#### read simulation results and put them in one list ##########################
filenames <- c("sanity_check")

all_data <- list()
for(filename in filenames){
  data <- read.csv(here("data", paste0(filename, "_growth_rates.csv")))
  for(name in unique(data$name)){
    all_data[[name]] <- data[data$name == name,]
  }
}


#### save max. mu and optimal xrP ##############################################
stats <- list()
for(name in names(all_data)){
  one <- all_data[[name]]
  stats[[name]] <- c(max(one$growth_rate), 
                     one$prot_fraction[which.max(one$growth_rate)])
}



#### Print max growth rates ####################################################
for(name in names(stats)){
  one <- stats[[name]]
  print(paste0("Max. growth rate ", name," (", one[1],") at ", one[2]*100,"% rP"))
}


all_data[["R_deg_hill_growth_rates"]] == all_data[["RBA_growth_rates"]]
plot(one$prot_fraction, one$growth_rate)

