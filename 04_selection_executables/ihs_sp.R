#!/usr/bin/env Rscript

print ("Starting ihs analysis")

#If run in CESGA, install packages before running the code in sbatch
  #install.packages("vcfR")
  #install.packages("rehh")
  #The downloaded source packages are in ‘/tmp/Rtmp7TLpRA/downloaded_packages’
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
  path<- "/mnt/lustre/scratch/home/csic/bie/llf/selection_scan/"
  
  #Read the subsetted vcf    
  print ("Reading vcf")  
  data<- data2haplohh(hap_file=paste0(path, chr, "_", sp, "_","goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf"), polarize_vcf= "FALSE", allele_coding="01", chr.name=chr, vcf_reader= "vcfR")
  
  #Scan the genome 
  print ("Scanning vcf")  
  data_frame<- scan_hh(data) 
  write.table(data_frame, file=(paste0(path, chr, "_", sp, "_", "scan")))
  print ("Saved scan data frame")
  
  df<- read.csv(paste0(path, chr, "_", sp, "_", "scan"), sep=" ")
  df_total<- rbind(df_total,df)
  print("Added one more chr in the total df")
}
#Save the sp_scan
write.table(df_total, file=(paste0(path, sp, "_scan")), sep="\t", row.names = FALSE, quote = FALSE)

#Calculate genome-wide iHS values
print ("Calculating ihs") 
wgscan.ihs<- ihh2ihs(df_total, freqbin = 0.05, p.adjust.method = "fdr") 
write.table(wgscan.ihs$ihs, file=(paste0(path, sp, "_ihs_scan")), sep="\t", row.names = FALSE, quote = FALSE)
print ("Saved ihs data frame")