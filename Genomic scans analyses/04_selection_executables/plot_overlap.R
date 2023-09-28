#!/usr/bin/env Rscript

#Install dependencies
#install.packages("tidyverse")
library(tidyverse)

#Set working directory
#in the terminal: setwd("/mnt/lustre/scratch/nlsas/home/csic/bie/llf/selection_scan")

#locally
setwd("/Users/lorenalorenzo/")
#Save the output
#sink(file = "plot_overlap_output.txt")


#set X axis with chr_size dataframe
path<- ("files/")  

axis_set <- read.table(paste0(path, "chr_size.txt"), sep="\t", col.names=c("chr", "size")) %>%
  summarize(axis_set, center = size/2)

#Group by chr
data_axis <- axis_set %>%
  mutate(size_cum = lag(cumsum(as.numeric(size)), default = 0)) %>%
  rowwise() %>%
  mutate(center_cum = sum(center, size_cum))

#Name the variables
species<- c("lc", "ll", "lp", "lr")

for (sp in species)
  {
    ###LASSI###
   path<- ("files/")  
  
    #Read the data
    lassi_99outliers <- read.table(paste0(path, sp, "_lassi_99outliers"), sep="\t", header = TRUE) 
    print("Reading lassi results done")
    
    #Group by chr
    data_lassi <- lassi_99outliers %>% 
      inner_join(data_axis, by = "chr") %>% 
      mutate(bp_cum = pos + size_cum)
    #uso end position of window para el cumulative o mejor position (medio de cada ventana)
    print("Acumulated by chr")
    
    ###iHS###
    #Read the data
    ihs_outliers <- read.table(paste0(path, sp, "_ihs_99.9outliers"), sep="\t", header= TRUE, col.names= c("chr", "pos", "ihs", "pval"))
    print("Reading iHS results done")
    
    #Group by chr
    data_ihs <- ihs_outliers %>% 
      inner_join(data_axis, by = "chr") %>% 
      mutate(bp_cum = pos + size_cum)
    
    ###intersect###
    #Read the data
    intersect_outliers <- read.table(paste0(path, sp, "_outliers_ihslassi_intersect.bed"), sep="\t", header= FALSE, col.names= c("chr", "start", "end", "n_wind", "n_snps"))
    print("Reading intersect results done")
    
    #Group by chr
    data_intersect <- intersect_outliers %>% 
      filter(n_snps > 0) %>%
      inner_join(data_axis, by = "chr") %>% 
      mutate(bp_cum = end + size_cum)
    
    
  ##############################Plot overlap########################################  
    
  path<- "plots/"  
    pdf(paste0(path, sp, "_overlap_manhattanplot.pdf"), width=10, height=5)
    
    scale= 6
    
    print(
      ggplot() +
        geom_point(data_ihs, mapping=aes(x= bp_cum, y=pval, alpha = 0.5, color = as_factor(chr)), shape=4) +
        geom_point(data_lassi, mapping=aes(x = bp_cum, y = statistic/scale, alpha = 0.5, color = as_factor(chr))) +
        geom_point(data_intersect, mapping=aes(x = bp_cum, y = -5, alpha = 0.5, color = as_factor(chr), size=n_snps), shape=16) +    
        scale_x_continuous(label = data_axis$chr, breaks = data_axis$center_cum) +
        scale_y_continuous(sec.axis = sec_axis(~./scale, name="Lassi statistic value")) +
        labs (x="Chromosome", y="iHS adjusted p-value", title=sp ) +
        theme(plot.title = element_text(hjust = 1, vjust = -10), legend.position = "none") +
        scale_color_manual(values=rep(c("black","grey"), 18 ))
    )
      
    dev.off()
    
   print("Plot done and saved")
  } 


