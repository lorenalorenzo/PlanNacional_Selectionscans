#!/usr/bin/env Rscript

print ("Starting ihs analysis")

#install.packages("vcfR")
#install.packages("rehh")
library(vcfR)
library(rehh)

#ENABLE command line arguments
sp <- commandArgs(trailingOnly = TRUE)

#sp<- c("lc", "ll", "lp", "lr") ###if I don't run the entire code
#Create a df object empty
df_total <- data.frame()

#Define chromosome variable
chromosomes <- c("A1", "A2", "A3", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "F1", "F2")

#loop to bind each chr_sp_scan df in only one
for (chr in chromosomes)
{
  print("Starting loop")
  path<- "/home/llorenzo/selection_scan/rehh_test/chr_ihs/"
  df<- read.csv(paste0(path, chr, "_", sp, "_", "scan"), sep=" ")
  df_total<- rbind(df_total,df)
  print("Added one more chr in the total df")
}
write.table(df_total, file=(paste0(sp, "_scan")), sep="\t", row.names = FALSE, quote = FALSE)

#Calculate genome-wide iHS values
print ("Calculating ihs") 
wgscan.ihs<- ihh2ihs(df_total, freqbin = 0.05, p.adjust.method = "fdr") 
write.table(wgscan.ihs$ihs, file=(paste0(sp, "_ihs_scan")), sep="\t", row.names = FALSE, quote = FALSE)

#Plotting results
#pdf(paste0(sp, "_", "ihs_manhattanplot.pdf"))
#manhattanplot(wgscan.ihs)
#dev.off()

pdf(paste0(sp, "_", "ihs_manhattanplot_p.adjust.pdf"))
manhattanplot(wgscan.ihs, pval=TRUE)
dev.off()
print ("Saved plots")  