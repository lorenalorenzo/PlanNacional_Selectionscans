#!/usr/bin/env Rscript

print ("Starting xpehh analysis")

#install.packages("vcfR")
#install.packages("rehh")
library(vcfR)
library(rehh)


#Read per sp scan
print("Reading lc scan")
lc_scan<- read.csv("lc_scan", sep=" ")
print("Reading ll scan")
ll_scan<- read.csv("ll_scan", sep=" ")
print("Reading lp scan")
lp_scan<- read.csv("lp_scan", sep=" ")
print("Reading lr scan")
lr_scan<- read.csv("lr_scan", sep=" ")

#Calculating XP-EHH
print ("Calculating xpehh") 

lc_ll_xpehh <- ies2xpehh(scan_pop1 = lc_scan, scan_pop2 = ll_scan, popname1 = "LC", popname2 = "LL", p.adjust.method = "fdr")
lc_lp_xpehh <- ies2xpehh(scan_pop1 = lc_scan, scan_pop2 = lp_scan, popname1 = "LC", popname2 = "LP", p.adjust.method = "fdr")
lc_lr_xpehh <- ies2xpehh(scan_pop1 = lc_scan, scan_pop2 = lr_scan, popname1 = "LC", popname2 = "LR", p.adjust.method = "fdr")
ll_lp_xpehh <- ies2xpehh(scan_pop1 = ll_scan, scan_pop2 = lp_scan, popname1 = "LL", popname2 = "LP", p.adjust.method = "fdr")
ll_lr_xpehh <- ies2xpehh(scan_pop1 = ll_scan, scan_pop2 = lr_scan, popname1 = "LL", popname2 = "LR", p.adjust.method = "fdr")
lp_lr_xpehh <- ies2xpehh(scan_pop1 = lp_scan, scan_pop2 = lr_scan, popname1 = "LP", popname2 = "LR", p.adjust.method = "fdr")

#Plotting results
print("Plotting results")

pdf(paste0("lc_ll_", "xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lc_ll_xpehh, pval=TRUE)
dev.off()

pdf(paste0("lc_lp_", "xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lc_lp_xpehh, pval=TRUE)
dev.off()

pdf(paste0("lc_lr_", "xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lc_lr_xpehh, pval=TRUE)
dev.off()

pdf(paste0("ll_lp_", "xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(ll_lp_xpehh, pval=TRUE)
dev.off()

pdf(paste0("ll_lr_", "xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(ll_lr_xpehh, pval=TRUE)
dev.off()

pdf(paste0("lp_lr_", "xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lp_lr_xpehh, pval=TRUE)
dev.off()

print ("Saved plots")  
