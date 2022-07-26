#!/usr/bin/env Rscript

print("Starting plotting LASSI results code")

#Install dependencies

##install.packages("dplyr")
##install.packages("tidyverse")
##install.packages("ggplot2")

library(dplyr)
library(tidyverse)
library(ggplot2)

#Set working directory
setwd("rehh_test")

#Name variables
species<- c("lc", "ll", "lp", "lr")


for (sp in species)
{
  #Read data
  df<- read.table(paste0(sp, "_salti.lassip.hap.out"), sep="\t", header=T)
  print(paste("Read data for", sp))
  
  #Group by chr
  data_cum <- df %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(pos)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(chr, bp_add)

  data <- df %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = pos + bp_add)

  axis_set <- data %>% 
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))
 
  print("Acumulated by chr")
  
  #Change the "sp_L" column name to another in order to make it easier to manage
  colnames(data)[12]<- "statistic"
  
  print("rename sp_L column to statistic")
  
  #get the 95% and 99% cut-off of the data
  print(paste(sp, round(length(data$statistic) * 0.05), "top 5% values"))
  print(paste(sp, round(length(data$statistic) * 0.01), "top 1% values"))
  
  cut_95 <- sort(data$statistic)[round(length(data$statistic) * 0.95)]
  cut_99 <- sort(data$statistic)[round(length(data$statistic) * 0.99)]
  
  print(paste(sp, cut_95, "5% cut-off" ))
  print(paste(sp, cut_99, "1% cut-off" ))
  
  print("Cut-offs done")
  
  #Extract a dataframe with the 1% outliers
  outliers<- data %>% filter(statistic >= cut_99) %>% select(c(1,2,3,5,12))
  write.table(outliers, file=(paste0(sp, "_99outliers")), sep="\t", row.names = FALSE, quote = FALSE)
  
  #Save as 0-based bed file
  bed <- outliers %>% mutate(start= start - 1) %>% select(c(1,2,3))
  write.table(bed, file=(paste0(sp, "_99outliers.bed")), sep="\t", col.names= FALSE, row.names = FALSE, quote = FALSE)
  
  print("outliers df saved")
}  
  ##########################Plotting results###################################
  
  #Give a name and size to the output plot
  pdf(paste0(sp,"_saltilassi_cutoffs_manhattanplot.pdf"), width=10, height=5)

  
  #Plot with cutoffs
  print(ggplot(data, aes(x = bp_cum, y = statistic, color = as_factor(chr))) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept= cut_95, color="grey", linetype="dashed") +
      geom_hline(yintercept= cut_99, color="red", linetype="dashed") +
      scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
      labs (x="Chromosome", y=expression(Lambda), title=sp ) +
      theme(plot.title = element_text(hjust = 1, vjust = -10), legend.position = "none"))
  
  #Another option if we don't want to use pdf() is to name the plot as an object and then:
  #ggsave(manhattanplot, file=paste0(sp,"_saltilassi_manhattanplot.pdf"), width=10, height=5)
  
  #Save the pdf results        
  dev.off()
  print("Manhattanplot done and saved")


  #Plot statistic distribution
  pdf(paste0(sp,"_L_statistic_distplot.pdf"), width=10, height=5)
  
  print(ggplot(data, aes(x= statistic)) +
          geom_freqpoly(binwidth=10) + 
          ylim(0, 20000) +
          labs (title=sp))
  
  dev.off()
  
  print("Distributionplot done and saved")
}



