#!/usr/bin/env Rscript

print ("Starting xpehh analysis")

#install.packages("vcfR")
#install.packages("rehh")
library(vcfR)
library(rehh)

#dependencies
path<- "/mnt/lustre/scratch/home/csic/bie/llf/selection_scan/"
#Name variables
species<- c("lc", "ll", "lp", "lr")

#Read per sp scan
  for(sp in species) 
    {
    print(paste("Reading", sp, "scan"))
    scan<- read.csv(paste0(path, sp, "_scan"), sep="\t")
    assign(paste0(sp, "_scan"), scan)
  }

#Calculating XP-EHH
print ("Calculating xpehh") 

lc_ll_xpehh <- ies2xpehh(scan_pop1 = lc_scan, scan_pop2 = ll_scan, popname1 = "LC", popname2 = "LL", p.adjust.method = "fdr")
write.table(lc_ll_xpehh, file=(paste0(path, "lc_ll_xpehh_scan")), sep="\t", row.names = FALSE, quote = FALSE)

lc_lp_xpehh <- ies2xpehh(scan_pop1 = lc_scan, scan_pop2 = lp_scan, popname1 = "LC", popname2 = "LP", p.adjust.method = "fdr")
write.table(lc_lp_xpehh, file=(paste0(path, "lc_lp_xpehh_scan")), sep="\t", row.names = FALSE, quote = FALSE)


lc_lr_xpehh <- ies2xpehh(scan_pop1 = lc_scan, scan_pop2 = lr_scan, popname1 = "LC", popname2 = "LR", p.adjust.method = "fdr")
write.table(lc_lr_xpehh, file=(paste0(path, "lc_lr_xpehh_scan")), sep="\t", row.names = FALSE, quote = FALSE)


ll_lp_xpehh <- ies2xpehh(scan_pop1 = ll_scan, scan_pop2 = lp_scan, popname1 = "LL", popname2 = "LP", p.adjust.method = "fdr")
write.table(ll_lp_xpehh, file=(paste0(path, "ll_lp_xpehh_scan")), sep="\t", row.names = FALSE, quote = FALSE)


ll_lr_xpehh <- ies2xpehh(scan_pop1 = ll_scan, scan_pop2 = lr_scan, popname1 = "LL", popname2 = "LR", p.adjust.method = "fdr")
write.table(ll_lr_xpehh, file=(paste0(path, "ll_lr_xpehh_scan")), sep="\t", row.names = FALSE, quote = FALSE)


lp_lr_xpehh <- ies2xpehh(scan_pop1 = lp_scan, scan_pop2 = lr_scan, popname1 = "LP", popname2 = "LR", p.adjust.method = "fdr")
write.table(lp_lr_xpehh, file=(paste0(path, "lp_lr_xpehh_scan")), sep="\t", row.names = FALSE, quote = FALSE)

#Plotting results
print("Plotting results")

pdf(paste0(path, "lc_ll_xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lc_ll_xpehh, pval=TRUE)
dev.off()

pdf(paste0(path, "lc_lp_xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lc_lp_xpehh, pval=TRUE)
dev.off()

pdf(paste0(path, "lc_lr_xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lc_lr_xpehh, pval=TRUE)
dev.off()

pdf(paste0(path, "ll_lp_xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(ll_lp_xpehh, pval=TRUE)
dev.off()

pdf(paste0(path, "ll_lr_xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(ll_lr_xpehh, pval=TRUE)
dev.off()

pdf(paste0(path, "lp_lr_xpehh_manhattanplot_p.adjust.pdf"))
manhattanplot(lp_lr_xpehh, pval=TRUE)
dev.off()

print ("Saved plots")  