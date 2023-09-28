#!/usr/bin/env Rscript

#install.packages
library(tidyverse)
library(rehh)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

#Set directories
path<- "/Users/lorenalorenzo"
files<- "/files/"
plots<- "/plots/"

#Name the variables
species<- c("lc", "ll", "lp", "lr")


for (sp in species)
{

# iHS results -------------------------------------------------------------
  ihs_scan <- read.table(paste0(path, files, sp , "_ihs_scan"), sep="\t", header= TRUE, col.names = c("chr", "pos", "ihs", "pval"))
  print("Read ihs_scan data")
  
  #Group by chr
  data_cum <- ihs_scan %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(pos)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(chr, bp_add)
  
  data_ihs <- ihs_scan %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = pos + bp_add) %>%
    na.omit() %>%
    mutate(abs_ihs= abs(ihs)) %>% #absolute value score of iHS
    mutate (dataset="|iHS|")  #for facet_grid
  
  axis_set <- data_ihs %>% 
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))
  print("Data set")
  
  write.table(data_ihs, file=(paste0(path, files, sp, "data_ihs")), sep="\t", row.names = FALSE, quote = FALSE)
  
  #explore ihs distribution (gaussian)
  # pdf(paste0(path, plots, sp,"_distplot.pdf"), width=10, height=5)
  # distribplot(ihs_scan$IHS, xlab= paste(sp, "iHS"), qqplot = TRUE)
  #dev.off()
  
## define outliers ----------------------------------------------------------------
  
  ihs_outliers<- data_ihs %>% filter(abs_ihs >= 4)
  print(length(ihs_outliers$abs_ihs))
  
  write.table(ihs_outliers, file=(paste0(path, files, sp, "_ihs_4outliers")), sep="\t", row.names = FALSE, quote = FALSE)
  print("outliers df saved")
  
  #Save as 0-based bed file
  bed_4 <- ihs_outliers %>% mutate(start= pos - 1) %>% select(c(1,9,2,3,7))
  write.table(bed_4, file=(paste0(path, files, sp, "_ihs_4outliers.bed")), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  
# saltiLASSI results -----------------------------------------------------------
  print("Starting plotting LASSI results code")

  #Read data
  lassi_scan<- read.table(paste0(path, files, sp, "_salti.lassip.hap.out"), sep="\t", header=T)
  print(paste("Read data for", sp))
  
  #Group by chr
  data_cum <- lassi_scan %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(pos)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(chr, bp_add)
  
  data_lassi <- lassi_scan %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = pos + bp_add) %>%
    mutate(dataset="saltiLASSI")
  
  axis_set <- data_lassi %>% 
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))
  
  print("Acumulated by chr")
  
  
  #Change the "sp_L" column name to another in order to make it easier to manage
  colnames(data_lassi)[12]<- "statistic"
  
  print("rename sp_L column to statistic")
  
  write.table(data_lassi, file=(paste0(path, files, sp, "data_lassi")), sep="\t", row.names = FALSE, quote = FALSE)
  
  #Explore statistic distribution
    #pdf(paste0(path, plots, sp,"_L_statistic_distplot.pdf"), width=10, height=5)
  
    #print(ggplot(data_lassi, aes(x= statistic)) +
         # geom_freqpoly(binwidth=10) + 
         # ylim(0, 20000) +
        #  labs (title=sp))
  
   #dev.off()
  
    #print("distribution plot saved")
  
## define outliers ---------------------------------------------------------
  #cut-off 99%
  print(paste(sp, round(length(data_lassi$statistic) * 0.01), "top 1% values"))
  cut_99 <- sort(data_lassi$statistic)[round(length(data_lassi$statistic) * 0.99)]
  print(paste(sp, cut_99, "1% cut-off" ))
  
  lassi_outliers_99<- data_lassi %>% filter(statistic >= cut_99)
  write.table(lassi_outliers_99, file=(paste0(path, files, sp, "_lassi_1%outliers")), sep="\t", row.names = FALSE, quote = FALSE)
  print("outliers df saved")
  
  bed_99 <- lassi_outliers_99  %>% select(c(1,2,3,12))
  write.table(bed_99, file=(paste0(path, files, sp, "_lassi_1%outliers.bed")), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  #cut-off 95%
  print(paste(sp, round(length(data_lassi$statistic) * 0.05), "top 5% values"))
  cut_95 <- sort(data_lassi$statistic)[round(length(data_lassi$statistic) * 0.95)]
  print(paste(sp, cut_95, "5% cut-off" ))
  
  lassi_outliers_95<- data_lassi %>% filter(statistic >= cut_95)
  write.table(lassi_outliers_95, file=(paste0(path, files, sp, "_lassi_5%outliers")), sep="\t", row.names = FALSE, quote = FALSE)
  print("outliers df saved")
  
  bed_95 <- lassi_outliers_95  %>% select(c(1,2,3,12))
  write.table(bed_95, file=(paste0(path, files, sp, "_lassi_5%outliers.bed")), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  

# plot iHS & lassi --------------------------------------------------------

#Give a name and size to the output plot
pdf(paste0(path, plots, sp,"_ihs_saltilassi_instersect_cutoffs_manhattanplot.pdf"), width=10, height=10)

#Plot together ihs & lassi
  print(
    ggplot() +
        geom_point(data_ihs, mapping=aes(x= bp_cum, y=abs_ihs, alpha = 0.5, color = as_factor(chr)), shape=15) + 
        geom_point(data_lassi, mapping=aes(x = bp_cum, y = statistic, alpha = 0.5, color = as_factor(chr)), shape=19) +    
        geom_hline(yintercept= 4, color="grey", linetype="dashed") +
        geom_hline(yintercept= -4, color="grey", linetype="dashed") +      
        scale_x_continuous(label = data_axis$chr, breaks = data_axis$center_cum) +
        labs (x="Chromosome", title=sp ) +
        theme(plot.title = element_text(hjust = 1, vjust = -10), legend.position = "none") +
        scale_color_manual(values=rep(c("#66C2A5", "#FC8D62"), 18 )) +
        facet_grid(dataset ~ ., scales = "free_y")
        )
#Save the pdf results        
  dev.off()

  print("plot done and saved")

}

## intersect ---------------------------------------------------------------
for (sp in species)
{
#Read the data from bedtools (selection_scan.Rmd)
intersect_outliers_99 <- read.table(paste0(path,files, sp, "_1%outliers_ihslassi_weighted.bed"), sep="\t", header= FALSE, col.names= c("chr", "start", "end", "n_wind", "n_snps", "weighted_statistic")) %>% 
  filter(n_snps > 0) %>%
  mutate(window_size= end - start) %>%
  mutate(weight= n_snps / window_size) 
print("Reading intersect results done")
print(paste(nrow(intersect_outliers_99), "outliers for", sp))
write.table(intersect_outliers_99, file=(paste0(path, files, sp, "_1%outliers_table")), sep="\t", row.names = FALSE, quote = FALSE)

intersect_outliers_95 <- read.table(paste0(path,files, sp, "_5%outliers_ihslassi_weighted.bed"), sep="\t", header= FALSE, col.names= c("chr", "start", "end", "n_wind", "n_snps", "weighted_statistic")) %>% 
  filter(n_snps > 0) %>%
  mutate(window_size= end - start) %>%
  mutate(weight= n_snps / window_size) 
print("Reading intersect results done")
print(paste(nrow(intersect_outliers_95), "outliers for", sp))
write.table(intersect_outliers_95, file=(paste0(path, files, sp, "_5%outliers_table")), sep="\t", row.names = FALSE, quote = FALSE)

#set X axis with chr_size dataframe
axis_set <- read.table(paste0(path, files, "chr_size.txt"), sep="\t", col.names=c("chr", "size")) %>%
  mutate(center = size/2)

#Group by chr
data_axis <- axis_set %>%
  mutate(size_cum = lag(cumsum(as.numeric(size)), default = 0)) %>%
  rowwise() %>%
  mutate(center_cum = sum(center, size_cum))    

data_intersect_99 <- intersect_outliers_99 %>% 
  filter(n_snps > 0) %>%
  inner_join(data_axis, by = "chr") %>% 
  mutate(bp_cum = end + size_cum) %>%
  mutate(dataset="intersect")
write.table(data_intersect_99, file=(paste0(path, files, sp, "_1%outliers_chr")), sep="\t", row.names = FALSE, quote = FALSE)

data_intersect_95 <- intersect_outliers_95 %>% 
  filter(n_snps > 0) %>%
  inner_join(data_axis, by = "chr") %>% 
  mutate(bp_cum = end + size_cum) %>%
  mutate(dataset="intersect")
write.table(data_intersect_95, file=(paste0(path, files, sp, "_5%outliers_chr")), sep="\t", row.names = FALSE, quote = FALSE)
}

## plot results --------------------------------------------------------
for (sp in species)
{
  #Read the data
  data_all <- read.table(paste0(path,files, sp, "_results_table_representation"), sep="\t", header= FALSE, na.strings = ".", col.names = c("chr", "start", "end", "pos", "statistic", "r_chr", "r_start", "r_end", "n_windows", "n_snps", "weighted_statistic", "window_size", "weight", "bp_overlap"))   %>% 
    inner_join(data_axis, by = "chr") %>% 
    mutate(bp_cum = end + size_cum) %>%
    transform(as.numeric(weight)) %>%
    mutate(weight= coalesce(weight, 0)) %>%
    mutate(intersect= ifelse(weight == 0, "lassi outliers", "intersect outliers"))
  
  #Use chr size to get cum_sum
   #set X axis with chr_size dataframe
    axis_set <- read.table(paste0(path, files, "chr_size.txt"), sep="\t", col.names=c("chr", "size")) %>%
                 mutate(center = size/2)
  
    #Group by chr
    data_axis <- axis_set %>%
                  mutate(size_cum = lag(cumsum(as.numeric(size)), default = 0)) %>%
                  rowwise() %>%
                  mutate(center_cum = sum(center, size_cum))      
 
  #set the cut-off
  print(paste(sp, round(length(data_all$statistic) * 0.05), "top 5% values"))
  cut_95 <- sort(data_all$statistic)[round(length(data_all$statistic) * 0.95)]
  print(paste(sp, cut_95, "5% cut-off" ))
  #name the plot
    #pdf(paste0(path, plots, "intersect_manhattanplot.pdf"), width=10, height=5)
  
    #print(

p<-  ggplot() +
      geom_point(data_all, mapping=aes(x = bp_cum, y = statistic, color = intersect, shape= intersect, alpha = weight, size= intersect)) +
      geom_hline(yintercept= cut_95, color="red", linetype="dashed")  +  
      scale_x_continuous(label = data_axis$chr, breaks = data_axis$center_cum) +
      labs (x="Chromosome") +
      scale_shape_manual(values=c(17, 1)) +
      scale_alpha(range = c(0.5, 1), guide= "none") +
      scale_color_manual(values= c("#FC8D62", "grey")) +
      scale_size_manual(values=c(3, 1)) 
  
  assign(paste0(sp, "_plot"), p )
  assign(paste0(sp, "_data"), data_all )
  
  
ggarrange(lc_plot, ll_plot, lp_plot, lr_plot, 
          labels = c("lynx canadensis", "lynx lynx", "lynx pardinus", "lynx rufus"), 
          ncol= 2, nrow= 2,
          common.legend= TRUE,
          legend="bottom")
  
  #save the plot
  dev.off()
  
  print(paste(sp, "plot done"))
}
  


  