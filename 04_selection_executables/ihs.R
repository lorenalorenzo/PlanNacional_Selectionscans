#!/usr/bin/env Rscript

print ("Starting ihs analysis")

#install.packages("vcfR")
#install.packages("rehh")
library(vcfR)
library(rehh)

#ENABLE command line arguments
chr <- commandArgs(trailingOnly = TRUE)

#path
setwd("$LUSTRE/selection_scan/chr_files/")

#Variables
species<- c("lc", "ll", "lp", "lr")

for (sp in species)
{
  #Read the subsetted vcf    
  print ("Reading vcf")  
  data<- data2haplohh(hap_file=paste0(chr, "_", sp, "_","goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf"), polarize_vcf= "FALSE", allele_coding="01", chr.name=chr, vcf_reader= "vcfR")
  
  #Scan the genome 
  print ("Scanning vcf")  
  data_frame<- scan_hh(data) 
  write.table(data_frame, file=(paste0(chr, "_", sp, "_", "scan")))
  print ("Saved scan data frame")
  
  #Calculate iHS
  print ("Calculating ihs")  
  ihs<- ihh2ihs(data_frame, freqbin = (round((min(data_frame$FREQ_A) * 2), 2)), p.adjust.method = "fdr") 

  #Plotting results
#  pdf(paste0(chr, "_", sp, "_", "ihs_manhattanplot.pdf"))
#  manhattanplot(ihs)
#  dev.off()
#  pdf(paste0(chr, "_", sp, "_", "ihs_manhattanplot_pvalue.pdf"))
#  manhattanplot(ihs, pval=TRUE)
#  dev.off()
#  print ("Saved plots")  
}