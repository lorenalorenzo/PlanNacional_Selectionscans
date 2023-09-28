#!/usr/bin/env Rscript

####Explore missingness distribution per sp####

#install.packages("tidyverse")
library(tidyverse)

#locally
setwd("/Users/lorenalorenzo/")


species<- c("lc", "ll", "lp", "lr")
samples<- c(19, 32, 11, 18)
table<- data.frame(cbind(species, samples))

for (i in 1:4)
{
  #Read the data 
  print("Reading missing data")
  data <- read.csv(paste0(table$species[i], "_missing_gts.csv"), header= T)
  
  print("Calculating missing data")
  miss_data<- data.frame(table(data$n_missing))
  miss_data$miss_prop= as.numeric(levels(miss_data$Var1))/as.numeric(table$samples[i])
  miss_data$freq_prop= ((miss_data$Freq)/NROW(data$n_missing))
  miss_data$freq_cumsum= cumsum(miss_data$freq_prop)

  print("Renaming missing data")
  assign(paste0(table$species[i], "_miss_data"), miss_data)
} 

#cum_miss <- 
ggplot() +
  geom_line(data= lc_miss_data, aes (x=miss_prop, y=freq_cumsum), alpha=0.5, color="blue", size=1.5) +
  geom_point(data= lc_miss_data, aes(x=miss_prop, y=freq_cumsum), alpha=0.5, color="blue", size=4, shape=15) +
  geom_line(data= ll_miss_data, aes (x=miss_prop, y=freq_cumsum), alpha=0.5, color="red", size=1.5) +
  geom_point(data= ll_miss_data, aes(x=miss_prop, y=freq_cumsum), alpha=0.5, color="red", size=4, shape=15) +
  geom_line(data= lp_miss_data, aes (x=miss_prop, y=freq_cumsum), alpha=0.5, color="green", size=1.5) +
  geom_point(data= lp_miss_data, aes(x=miss_prop, y=freq_cumsum), alpha=0.5, color="green", size=4, shape=15) +
  geom_line(data= lr_miss_data, aes (x=miss_prop, y=freq_cumsum), alpha=0.5, color="yellow", size=1.5) +
  geom_point(data= lr_miss_data, aes(x=miss_prop, y=freq_cumsum), alpha=0.5, color="yellow", size=4, shape=15) +
  scale_x_continuous() +
  xlab("Proportion of missing data") +
  ylab("Proportion of SNPs included") +
  theme_minimal()
 