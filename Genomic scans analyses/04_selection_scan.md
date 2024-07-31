Positive selection genomic scan on lynxes
================
Lorena Lorenzo
31-Jul-2024

## 1. Genomic scan methods

The aim of this chapter is to make an approach into positive selection
among lynxes species through genomic scan methods. After an exhaustive
revision of the state of art, we decided to focus on two footprints of
positive selection across genome:

1.  Extended Haplotype Homozigosity

2.  Haplotype Frequency Spectrum Distortion

### 1.1 Data set-up

Once I have filtered, phased and polarized data (see
variant_filtering.md , phasing.md or polarization.md), I can start with
selection analyses. For this step I am going to use both phased and
non-phased vcfs (phased with iHS analyses and non-phased with
saltiLASSI). Those files are:

`${sp}_goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf`

`${sp}_goodsamples_filtered_polarized_variants_header_cat_ref.vcf`

Moreover, after phasing, allele count is not recalculated so that I used
the following code for updating INFO field.

``` bash
species=(lc ll lp lr)

for sp in ${species[@]}
do
cat ${sp}_goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf | fill-an-ac | bgzip -c > ${sp}_goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf.gz
done
```

Both methods of selection need to be run per chr

``` bash

#Separate per chr sp vcf
species=(lc ll lp lr)

for sp in ${species[@]}
do
  echo "$sp"
  CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
  for chr in ${CHR[@]:0:18}
    do
      echo "$chr"
      grep -E "^(#|${chr})" \
      $INPUT_FILE \
      > $OUTPUT_FILE
      echo "Done $chr for $sp vcf"
    done
done
```

#### 1.1.1. Control check: exploring missingness per sp

Code available in missing_distribution.R and missing_distribution.sh

``` bash

for sp in ${species[@]}
  do
    echo "$sp"
    vcftools \
    --vcf ${sp}_goodsamples_cat_ref.filter8.vcf \
    --missing-site \
    --out ${sp}_missing_site_rate
  done
  
  bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' lr_goodsamples_cat_ref.filter8.vcf > lr_genotypes_cat_ref.filter8.vcf
```

### 1.2 iHS analyses

For iHS analyses we decided to use the “rehh” package available for R.
The entire code is available in ihs_sp.R. To run this in CESGA, I had to
create a run_ihs.sh script that specify to run the Rscript.

``` r
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

#Define paths
input_path<- "/mnt/lustre/scratch/nlsas/home/csic/bie/llf/selection_scan/chr_files/"
output_path<- "/mnt/lustre/scratch/nlsas/home/csic/bie/llf/selection_scan/iHS/"

#loop to bind each chr_sp_scan df in only one
for (chr in chromosomes)
{
  print("Starting loop")
 
  
  #Read the subsetted vcf    
  print ("Reading vcf")  
  data<- data2haplohh(hap_file=paste0(input_path, chr, "_", sp, "_","goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf"), polarize_vcf= "FALSE", allele_coding="01", chr.name=chr, vcf_reader= "vcfR")
  
  #Scan the genome 
  print ("Scanning vcf")  
  data_frame<- scan_hh(data) 
  write.table(data_frame, file=(paste0(output_path, chr, "_", sp, "_", "scan")))
  print ("Saved scan data frame")
  
  df<- read.csv(paste0(output_path, chr, "_", sp, "_", "scan"), sep=" ")
  df_total<- rbind(df_total,df)
  print("Added one more chr in the total df")
}

#Save the sp_scan
write.table(df_total, file=(paste0(output_path, sp, "_scan")), sep="\t", row.names = FALSE, quote = FALSE)

#Calculate genome-wide iHS values
print ("Calculating ihs") 
wgscan.ihs<- ihh2ihs(df_total, (round((min(df_total$FREQ_A) * 2), 2))) 
write.table(wgscan.ihs$ihs, file=(paste0(output_path, sp, "_ihs_scan")), sep="\t", row.names = FALSE, quote = FALSE)
print ("Saved ihs data frame")
```

### 1.2 XP-EHH analyses

As with iHS, I used rehh R package. The entire code is available in
xpehh_sp.R. To run this in CESGA, I had to create a run_xpehh.sh script
that specify to run the Rscript.

``` r
#!/usr/bin/env Rscript

print ("Starting xpehh analysis")

#install.packages("vcfR")
#install.packages("rehh")
library(vcfR)
library(rehh)

#Define paths
path<- "/mnt/lustre/scratch/nlsas/home/csic/bie/llf/selection_scan/iHS/"

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
```

### 1.3. salti-LASSI

Here I am going to use saltiLASSI to analyse selection from haplotype
frequency spectrum (HFS). I used windows of 101 SNPS with a 50 slide.
The entire code is available in lassi.sh.

``` bash

species=(lc ll lp lr)

for sp in ${species[@]}
  do
    echo "$sp"
    #make a sp_ind.txt (pop file) needed for lassip
    sed "s/$/\t${sp}/" <(grep "#CHR" $LUSTRE/${sp}_goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf | tr "\t" "\n" | grep ${sp})  > $LUSTRE/selection_scan/${sp}_ind.txt
    sbatch lassi.sh $sp
  done  
  
species=(lc ll lp lr)

for sp in ${species[@]}
  do
  sbatch lassi.sh $sp
  done
```

## 2. Results

### 2.1 Representation of iHS and salti-LASSI results

For graphics we are going to use R. First thing is to run packages
needed and set directories.

``` r
#install.packages
library(tidyverse)
library(rehh)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(grid)
library(ggrepel)

#Set directories
path<- "/Users/lorenalorenzo/PlanNacional_Selectionscans"
files<- "/files/"
plots<- "/plots/"

#Name the variables
species<- c("lc", "ll", "lp", "lr")

names <- c("lc" = "Lynx canadensis",
            "lr" = "Lynx rufus",
            "ll" = "Lynx lynx",
            "lp" = "Lynx pardinus")

chr <- c( "A1", "A2", "A3", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "F1", "F2")

#colors <- c("lc"= "#AB97B7",
#            "ll"= "#F2AB5A",
#            "lp"="#EB595F",
#            "lr"="#71C178")
#print(as.character(names[species]))
#print(as.character(colors[species]))
#colors <- c("lr"= "#ED90A4",
#            "lp"= "#ABB150",
#            "ll"="#00C1B2",
#            "lc"="#ACA2EC")

colors <- c("lr"= "#FFB3B5",
            "lp"= "#BBCF85",
            "ll"="#61D8D6",
            "lc"="#D4BBFC")

#set X axis with chr_size dataframe (for plotting)
    axis_set <- read.table(paste0(path, files, "chr_size.txt"), sep="\t", col.names=c("chr", "size")) %>%
                 mutate(center = size/2)
    #Group by chr
    data_axis <- axis_set %>%
                  mutate(size_cum = lag(cumsum(as.numeric(size)), default = 0)) %>%
                  rowwise() %>%
                  mutate(center_cum = sum(center, size_cum))     
```

In the following code I prepare iHS and saltiLASSI results for
representation with cutoffs.

``` r
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
  #cut-off 95%
  print(paste(sp, round(length(data_lassi$statistic) * 0.05), "top 5% values"))
  cut_95 <- sort(data_lassi$statistic)[round(length(data_lassi$statistic) * 0.95)]
  print(paste(sp, cut_95, "5% cut-off" ))
  
  lassi_outliers_95<- data_lassi %>% filter(statistic >= cut_95)
  write.table(lassi_outliers_95, file=(paste0(path, files, sp, "_lassi_5%outliers")), sep="\t", row.names = FALSE, quote = FALSE)
  print("outliers df saved")
  
  bed_95 <- lassi_outliers_95  %>% select(c(1,2,3,12))
  write.table(bed_95, file=(paste0(path, files, sp, "_lassi_5%outliers.bed")), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}  
```

    ## [1] "Read ihs_scan data"
    ## [1] "Data set"
    ## [1] 2992
    ## [1] "outliers df saved"
    ## [1] "Starting plotting LASSI results code"
    ## [1] "Read data for lc"
    ## [1] "Acumulated by chr"
    ## [1] "rename sp_L column to statistic"
    ## [1] "lc 1424 top 5% values"
    ## [1] "lc 5.76722 5% cut-off"
    ## [1] "outliers df saved"
    ## [1] "Read ihs_scan data"
    ## [1] "Data set"
    ## [1] 693
    ## [1] "outliers df saved"
    ## [1] "Starting plotting LASSI results code"
    ## [1] "Read data for ll"
    ## [1] "Acumulated by chr"
    ## [1] "rename sp_L column to statistic"
    ## [1] "ll 3080 top 5% values"
    ## [1] "ll 23.5728 5% cut-off"
    ## [1] "outliers df saved"
    ## [1] "Read ihs_scan data"
    ## [1] "Data set"
    ## [1] 1625
    ## [1] "outliers df saved"
    ## [1] "Starting plotting LASSI results code"
    ## [1] "Read data for lp"
    ## [1] "Acumulated by chr"
    ## [1] "rename sp_L column to statistic"
    ## [1] "lp 865 top 5% values"
    ## [1] "lp 6.53881 5% cut-off"
    ## [1] "outliers df saved"
    ## [1] "Read ihs_scan data"
    ## [1] "Data set"
    ## [1] 3880
    ## [1] "outliers df saved"
    ## [1] "Starting plotting LASSI results code"
    ## [1] "Read data for lr"
    ## [1] "Acumulated by chr"
    ## [1] "rename sp_L column to statistic"
    ## [1] "lr 8763 top 5% values"
    ## [1] "lr 10.0232 5% cut-off"
    ## [1] "outliers df saved"

![](figs/unnamed-chunk-9-1.png)<!-- -->

### 2.2 Intersect between methods

Once saved the outliers bed from iHS and saltiLASSI, we need the overlap
of both methods. For doing that, here we are using bedtools. Firstly, I
merge every overlapping outlier window from lassi. Then, we intersect
this windows with the outlier SNPs from iHS to get putative selected
genomic REGIONS. Use 5% as cutoff for lassi because after that we are
going to rank 10 best regions.

``` bash
module load bedtools

species=(lc ll lp lr)

for sp in ${species[@]}
  do
    echo "$sp"
       #merge every overlapping outlier window in LASSI 
       bedtools merge \
       -i ${sp}_lassi_5%outliers.bed \
       -c 1 \
       -o count \
       > ${sp}_merged_tmp
       echo "$sp merged"
       
       #intersect with iHS snps so the result are GENOMIC REGIONS under putative selection
       bedtools intersect \
       -a ${sp}_merged_tmp \
       -b ${sp}_ihs_4outliers.bed \
       -c \
       > ${sp}_intersected_tmp
       echo "$sp intersected"
       
       #add snp density        
       bedtools intersect \
       -a ${sp}_intersected_tmp \
       -b $LUSTRE/vcfs/${sp}_goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf \
       -c \
     > ${sp}_snps_density_tmp
     
       #calculate lassi values for those regions
       bedtools map \
       -a ${sp}_snps_density_tmp\
       -b ${sp}_lassi_5%outliers.bed \
       -F 1 \
       -c 4 \
       -o mean,max,sum \
       > ${sp}_lassi_tmp
       echo "$sp lassi statistics calculated"
       
       #calculate iHS values for those regions
       bedtools map \
       -a ${sp}_lassi_tmp \
       -b ${sp}_ihs_4outliers.bed \
       -F 1 \
       -c 5 \
       -o mean,max,sum \
       > ${sp}_lassi_ihs_tmp
       echo "$sp ihs statistics calculated"
       
       #remove 0's (no SNPS)
       awk -F '\t' '$5 != 0' ${sp}_lassi_ihs_tmp > ${sp}_lassi_ihs_regions_tmp
       echo "$sp regions with 0 snps intersected removed"
       
       #print a header (column names)
       echo -e "chr start end windows snps snps_density lassi_mean lassi_max lassi_sum ihs_mean ihs_max ihs_sum" | cat - ${sp}_lassi_ihs_regions_tmp | tr " " "\t" > ${sp}_lassi_ihs_regions 
       
  ######Cross bed results with annotation file (gff3)    
      #match selected regions with genome annotation
       bedtools intersect \
         -a ${sp}_lassi_ihs_regions_tmp \
         -b $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.97.gff3 \
         -wa -wb \
         > ${sp}_annotated_tmp 
        echo "$sp annotated"
        
      #filter only genes
      awk '$15=="gene"' ${sp}_annotated_tmp > ${sp}_gene_tmp
      echo "$sp gene filtered"
    
    #get only interesting columns
    cut -f-12,16,17,21- ${sp}_gene_tmp | tr ' ' '\t' > ${sp}_filtered_tmp
    
    #ensembl id and gene name as column data
    cut -d';' -f1,2 ${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($16 ~ /Name=[[:alnum:]]/) $16=$16; else $16="NA"; print $0}' > ${sp}_names_tmp 
    
    #cut extra-info from ens code and gene name columns
    paste <(cut -d" " -f-14 ${sp}_names_tmp) <(cut -d" " -f15 ${sp}_names_tmp | cut -d":" -f2) <(cut -d" " -f16 ${sp}_names_tmp | cut -d"=" -f2) |  tr "\t" " "  >> ${sp}_columns_tmp
    
    #print a header (column names)
    echo -e "chr start end windows snps snps_density lassi_mean lassi_max lassi_sum ihs_mean ihs_max ihs_sum gene_start gene_end ensembl_id gene_name" | cat - ${sp}_columns_tmp > ${sp}_genomic_regions_annotated
      
        #remove temporal files
       rm *tmp
    done
```

Let’s get some information about the genomic regions defined, such as
size range

![](figs/unnamed-chunk-11-1.png)<!-- -->![](figs/unnamed-chunk-11-2.png)<!-- -->![](figs/unnamed-chunk-11-3.png)<!-- -->![](figs/unnamed-chunk-11-4.png)<!-- -->

### 2.3 Representation of the intersection

For representation purposes, I want to merge lassi outliers windows data
with intersected outliers regions data.

``` bash
module load bedtools

species=(lc ll lp lr)

#Merge info about lassi results and the final table of genomics regions under selection
for sp in ${species[@]}
  do
    echo "$sp"
    bedtools intersect \
    -a <(grep -v "chr" saltiLASSI/${sp}_salti.lassip.hap.out | cut -f1,2,3,5,12) \
    -b <(grep -v "chr" ${sp}_genomic_regions | cut -f1,2,3,4,5,7,10) \
    -wao \
    > ${sp}_lassi_windows_genomic_regions
  done  
```

    ## [1] "lc 1424 top 5% values"
    ## [1] "lc 5.76722 5% cut-off"
    ## [1] "ll 3080 top 5% values"
    ## [1] "ll 23.5728 5% cut-off"
    ## [1] "lp 865 top 5% values"
    ## [1] "lp 6.53881 5% cut-off"
    ## [1] "lr 8763 top 5% values"
    ## [1] "lr 10.0232 5% cut-off"

![](figs/unnamed-chunk-13-1.png)<!-- -->

CLARIFICATION: In the plot we can see \* below the cutoff. That’s
because another overlapping window passed the cutoff value, so now the
region includes windows that by its own didn’t were outliers.

Now I want to explore regions defined in a chr

![](figs/unnamed-chunk-14-1.png)<!-- -->![](figs/unnamed-chunk-14-2.png)<!-- -->![](figs/unnamed-chunk-14-3.png)<!-- -->![](figs/unnamed-chunk-14-4.png)<!-- -->

### 2.3 Representation of XP-EHH

    ## [1] "43 outliers lc_ll_xpehh_scan"
    ## [1] "xpehh outliers df saved"
    ## [1] "38 outliers lc_lp_xpehh_scan"
    ## [1] "xpehh outliers df saved"
    ## [1] "54 outliers lc_lr_xpehh_scan"
    ## [1] "xpehh outliers df saved"
    ## [1] "10 outliers ll_lp_xpehh_scan"
    ## [1] "xpehh outliers df saved"
    ## [1] "64 outliers ll_lr_xpehh_scan"
    ## [1] "xpehh outliers df saved"
    ## [1] "11 outliers lp_lr_xpehh_scan"
    ## [1] "xpehh outliers df saved"

![](figs/unnamed-chunk-15-1.png)<!-- -->

Let’s see the annotation on those regions

``` bash

file_list=$(ls files/*xpehh*outliers)

for i in ${file_list[@]}
  do
    echo "$i"
    awk 'NR>1 {print $1"\t"$2-1"\t"$2}' ${i} | sed 's/chr//g' > ${i}.bed
    
    #merge every overlapping outlier in xpehh
    bedtools merge \
      -d 100000 \
      -i ${i}.bed \
      > ${i}_genomic_regions
      echo "$i merged"

######Cross bed results with annotation file (gff3)    
    #match selected regions with genome annotation
    bedtools intersect \
     -a ${i}_genomic_regions \
     -b files/Felis_catus.Felis_catus_9.0.97.gff3 \
     -wa -wb \
      > ${i}_annotated 
     echo "$i annotated"
  done   
```

``` bash
# intersection w/ intrasp. outliers (ihs&lassi)
    bedtools intersect \
    -wa \
    -a <(grep -v "chr" files/ll_genomic_regions_annotated | tr " " "\t") \
    -b ${i}.bed \
    > ll_lp_xpehh_intersected_ll_lassi_ihs_outliers
```

## 3. Interpretation

### 3.1 Looking for enrichment pathways

We tried to find out pathways or functions specifically being under
selection pressure among the putative selected regions (those
overpassing the filter for both methods). Unfortunately, GO enrichment
analyses done in PANTHER results in no enrichment pathway in no species.

``` r
##Dependencies
#BiocManager::install("biomaRt")
library(topGO)
library(biomaRt)

#####Get annotation from ensembl.org#####

#read from ensembl.org every felcat ensembl ID:
ensembl <- useMart("ensembl", dataset = "fcatus_gene_ensembl")  
###extract the GOterms for every ensembl_id
##ensembl_to_go <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id"), ##mart = ensembl)

# Get attributes
attributes <- c("ensembl_gene_id", "external_gene_name", "go_id", "gene_biotype")

# Retrieve annotations
ensembl_to_go <- getBM(attributes = attributes, mart = ensembl) %>%
  filter(gene_biotype== "protein_coding")

#make a list with every GO terms per ENSEMBL_ ID. 
go_list <- split(ensembl_to_go$go_id, ensembl_to_go$ensembl_gene_id)

  #comment: Without filtering for coding genes, there is a total of 29550 ensembl_id corresponding to coding genes (19588) + non-coding genes (9468) + pseudogenes (494) according to https://www.ensembl.org/Felis_catus/Info/Annotation. We are getting 19564, 24 less than reported by the assembly (don't really know why)


for (sp in species)
{

#####read my gene data set#####  
    df_genes<- read.table(paste0(path, files, sp, "_genomic_regions_annotated"), 
                       sep=" ", header=T)
    genes <- df_genes$ensembl_id #this is a character object with the ensembl_ids ("ENSFCAG00000008235" "ENSFCAG00000008236" "ENSFCAG00000029529" ...)
    assign(paste0(sp, "_genes"), genes)

#cross the felcat annotation (ensembl_id with its GOterms) with my set of genes, to get a logical factor of true (if the gene is in the set) and false (if its not)
allgenes = factor(as.integer(names(go_list) %in% genes))
names(allgenes) <- names(go_list) #ensure names of allgenes are names of go_list (that is, the ensembl_id)

assign(paste0(sp, "_allgenes"), allgenes)

#create the empty results df before the loop
final_overrep <- data.frame()

for (gocat in c("BP", "MF", "CC"))
{
  
#####Creating the topGOdata object#####
godata <- new("topGOdata", ontology = gocat, allGenes = allgenes, 
              annotationFun = annFUN.gene2GO, gene2GO = go_list, nodeSize= 5)

assign(paste0(sp, gocat, "_godata"), godata)
#####run the overrepresentation test#####
over_test <- runTest(godata, statistic = "fisher")

assign(paste0(sp, "_test"), over_test)

result_table <- GenTable(godata, Fisher=over_test, topNodes=over_test@geneData[2], numChar=1000) %>% 
        as_tibble() %>% 
        mutate(p.adj = round(p.adjust(as.numeric(gsub("<", "", Fisher)), method="BH"), 15), ontology= gocat) %>% 
        filter(p.adj<0.05 & Significant >1) 

final_overrep  <- rbind (final_overrep, result_table ) 

}
assign(paste0(sp, "_results"), final_overrep ) %>% mutate(species= names[sp])

}

write.table(lc_results, file ="lc_functional_enrichment.csv", sep =  "\t", row.names = FALSE, quote = FALSE)
write.table(ll_results, file ="ll_functional_enrichment.csv", sep =  "\t", row.names = FALSE, quote = FALSE)
write.table(lp_results, file ="lp_functional_enrichment.csv", sep =  "\t", row.names = FALSE, quote = FALSE)
write.table(lr_results, file ="lr_functional_enrichment.csv", sep =  "\t", row.names = FALSE, quote = FALSE)

#to print in the Rmd
print(lc_results)
```

    ## # A tibble: 5 × 8
    ##   GO.ID      Term        Annotated Significant Expected Fisher    p.adj ontology
    ##   <chr>      <chr>           <int>       <int>    <dbl> <chr>     <dbl> <chr>   
    ## 1 GO:0007186 G protein-…      1342          44    18.3  1.1e-… 1.24e- 6 BP      
    ## 2 GO:0050911 detection …       451          24     6.14 1.2e-… 1.24e- 6 BP      
    ## 3 GO:0004984 olfactory …       687          38     9.27 8.0e-… 1.62e-11 MF      
    ## 4 GO:0004930 G protein-…      1069          42    14.4  3.4e-… 3.45e- 8 MF      
    ## 5 GO:0005886 plasma mem…      3851          82    53.4  1.4e-… 3.02e- 4 CC

``` r
print(ll_results)
```

    ## # A tibble: 5 × 8
    ##   GO.ID      Term          Annotated Significant Expected Fisher  p.adj ontology
    ##   <chr>      <chr>             <int>       <int>    <dbl> <chr>   <dbl> <chr>   
    ## 1 GO:1902236 negative reg…        16           2     0.03 0.000… 0.0126 BP      
    ## 2 GO:0001755 neural crest…        37           2     0.07 0.002… 0.0126 BP      
    ## 3 GO:0031648 protein dest…        40           2     0.08 0.002… 0.0126 BP      
    ## 4 GO:0043632 modification…       501           3     0.99 0.007… 0.0126 BP      
    ## 5 GO:0045927 positive reg…       166           2     0.33 0.007… 0.0126 BP

``` r
print(lp_results)
```

    ## # A tibble: 6 × 8
    ##   GO.ID      Term        Annotated Significant Expected Fisher    p.adj ontology
    ##   <chr>      <chr>           <int>       <int>    <dbl> <chr>     <dbl> <chr>   
    ## 1 GO:1902236 negative r…        16           3     0.14 0.000… 4.28e- 2 BP      
    ## 2 GO:0006334 nucleosome…        77           5     0.69 0.000… 4.28e- 2 BP      
    ## 3 GO:0030527 structural…        90          14     0.87 1.7e-… 2.48e-11 MF      
    ## 4 GO:0046982 protein he…       276          13     2.68 2.9e-… 2.12e- 4 MF      
    ## 5 GO:0003677 DNA binding      1808          35    17.6  8.9e-… 4.33e- 3 MF      
    ## 6 GO:0000786 nucleosome        118          14     1.14 7.8e-… 1.18e- 9 CC

``` r
print(lr_results)
```

    ## # A tibble: 4 × 8
    ##   GO.ID      Term          Annotated Significant Expected Fisher  p.adj ontology
    ##   <chr>      <chr>             <int>       <int>    <dbl> <chr>   <dbl> <chr>   
    ## 1 GO:0019441 tryptophan c…         5           3     0.08 4.6e-… 0.0117 BP      
    ## 2 GO:0006067 ethanol meta…         7           3     0.12 0.000… 0.0204 BP      
    ## 3 GO:0032230 positive reg…         8           3     0.13 0.000… 0.0213 BP      
    ## 4 GO:0004022 alcohol dehy…         6           3     0.11 0.000… 0.0299 MF

``` r
# Define a custom color palette for ontology
ontology_colors <- c("BP" = "#b3e2cd",
                     "MF" = "#fdcdac",
                     "CC" = "#cbd5e8")

ontology_names <- c("BP" = "Biological Process",
                    "MF" = "Molecular Function",
                    "CC" = "Cellular Component")

# Create the barplot for lc_results
for (sp in species)
{
  # Read the results data frame
  results <- read.csv(paste0(sp, "_functional_enrichment.csv"), sep = "\t") %>%
              mutate(Gene_Ratio = Significant / Expected) %>%
              mutate(qscore = -log(p.adj, base=10)) %>%
              arrange(ontology, desc(qscore)) %>%
              mutate(Term = factor(Term, levels = rev(unique(Term))))

  # Create the barplot for all results
  enrichment <- ggplot(results, aes(x=qscore, y=Term, fill=ontology)) +
    geom_bar(stat="identity") +
    labs(x="qscore", y="GO term", title=as.character(names[sp])) +
    scale_fill_manual(name = "Ontology",
                      values = ontology_colors,
                      labels = ontology_names) +
    theme (plot.title = element_text(face = "italic", size= 20),
           axis.text.y = element_text(size = 15))  

  assign(paste0(sp, "_enrichment"), enrichment )

}
## customize
q<- ggarrange(lc_enrichment, ll_enrichment, lp_enrichment, lr_enrichment, 
     ncol= 1, nrow= 4, common.legend= TRUE, legend="bottom") 

print(q)
```

![](figs/unnamed-chunk-18-1.png)<!-- -->

``` r
ggsave(filename = paste0(path, plots, "enrichment.png"), plot= q, width= 15, height = 20)
```

### 3.2 Top regions

We’ve decided to look for the most important regions under selection in
each species. For doing so, we had to decided what parameter describes
better which regions are “the most important ones”. For so, we used
lassi max.value.

Clarification: When trying the functional enrichment in
<https://pantherdb.org>, I didn’t get any result because the database
didn’t find the genes I was providing, so that I did it with topGO as
seen above.

From the resulting table I get the top 10 putative regions per sp.

``` r
#Read data
for (sp in species)
{
data<- read.table(paste0(path, files, sp, "_lassi_ihs_regions"), sep="\t", header=T)
assign(paste0(sp, "_results"), data )

data_top10<- data %>%
  arrange(desc(lassi_max)) %>%
  dplyr::slice(1:10)
#write.table(data_top10, file=(paste0(path, files, sp, "_top10")), sep="\t",row.names = FALSE, quote = FALSE)

assign(paste0(sp, "_top10"), data_top10 )
}
```

In order to explore the distribution of the lassi statistic in each of
the 10 top regions, I firstly get the statistic value in every window of
the top region and then represent it.

``` bash
module load bedtools

species=(lc ll lp lr)

for sp in ${species[@]}
  do
    echo "$sp"
    #convert into a bed file
    awk '{printf("%s\t%d\t%d\t%d\t%d\n", $1, $2, $3, $5, $12)}' saltiLASSI/${sp}_salti.lassip.hap.out  > ${sp}_tmp

    #intersect btw lassi values and genomic regions with 100% coincidence (only windows embedded in the genomic regions are reported)
       bedtools intersect -wa -wb \
       -a <(grep -v "chr" ${sp}_tmp) \
       -b <(grep -v "chr" ${sp}_top10 | cut -f1-3) \
       -f 1.0 \
       > ${sp}_lassi_top10  
       echo "$sp intersected"
    
  #remove temporal files
  rm *tmp
  done
```

Now I want to represent each top10 genomic regions and the correspondent
lassi values. Taking into account I previously calculated the genes
under the genomic outlier regions, we can add this info to the plot

![](figs/unnamed-chunk-21-1.png)<!-- -->![](figs/unnamed-chunk-21-2.png)<!-- -->![](figs/unnamed-chunk-21-3.png)<!-- -->![](figs/unnamed-chunk-21-4.png)<!-- -->

I see there is a conspicuous pattern between species of selection
outliers in chr D3:

![](figs/unnamed-chunk-22-1.png)<!-- -->

## 4. Repeated selection ocurring in several sp. in the genus

This approach is used to identify genes used repeatedly for adaptation
based on genomic scans results.

``` bash
bedtools intersect -wa -wb \
       -a <(cut -f 1-3 files/lc_lassi_ihs_regions) \
       -b <(cut -f 1-3 files/ll_lassi_ihs_regions) <(cut -f 1-3 files/lp_lassi_ihs_regions) <(cut -f 1-3 files/lr_lassi_ihs_regions) \
       -names ll lp lr \
       -sorted \
       > files/lc_repeated_regions
       
bedtools intersect -wa -wb \
       -a <(cut -f 1-3 files/ll_lassi_ihs_regions) \
       -b <(cut -f 1-3 files/lc_lassi_ihs_regions) <(cut -f 1-3 files/lp_lassi_ihs_regions) <(cut -f 1-3 files/lr_lassi_ihs_regions) \
       -names lc lp lr \
       -sorted \
       > files/ll_repeated_regions      
       
bedtools intersect -wa -wb \
       -a <(cut -f 1-3 files/lp_lassi_ihs_regions) \
       -b <(cut -f 1-3 files/lc_lassi_ihs_regions) <(cut -f 1-3 files/ll_lassi_ihs_regions) <(cut -f 1-3 files/lr_lassi_ihs_regions) \
       -names lc ll lr \
       -sorted \
       > files/lp_repeated_regions  
       
bedtools intersect -wa -wb \
       -a <(cut -f 1-3 files/lr_lassi_ihs_regions) \
       -b <(cut -f 1-3 files/lc_lassi_ihs_regions) <(cut -f 1-3 files/ll_lassi_ihs_regions) <(cut -f 1-3 files/lp_lassi_ihs_regions) \
       -names lc ll lp \
       -sorted \
       > files/lr_repeated_regions  
       
######Cross bed results with annotation file (gff3)    
    #match selected regions with genome annotation
species=(lc ll lp lr)

for sp in ${species[@]}
  do
    echo "$sp"
    bedtools intersect \
     -a files/${sp}_repeated_regions \
     -b files/Felis_catus.Felis_catus_9.0.97.gff3 \
     -wa -wb \
      > files/${sp}_annotated_tmp 
     echo "$sp annotated"
     #filter only genes
     awk '$10=="gene"' files/${sp}_annotated_tmp  > files/${sp}_gene_tmp
     echo "$sp gene filtered"   
     #get only interesting columns
    cut -f 1-8,11-12,16 files/${sp}_gene_tmp > files/${sp}_filtered_tmp
    #ensembl id and gene name as column data
    cut -d';' -f1,2 files/${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($12 ~ /Name=[[:alnum:]]/) $12=$12; else $12="NA"; print $0}' > files/${sp}_rep_genes_names_tmp 
    
    #cut extra-info from ens code and gene name columns
    paste <(cut -d" " -f-10 files/${sp}_rep_genes_names_tmp) <(cut -d" " -f11 files/${sp}_rep_genes_names_tmp | cut -d":" -f2) <(cut -d" " -f12 files/${sp}_rep_genes_names_tmp | cut -d"=" -f2) |  tr "\t" " "  >> files/${sp}_rep_columns_tmp
    
    #print a header (column names)
    echo -e "chr start end sp chr_2 start_2 end_2 gene_chr gene_start gene_end ensembl_id gene_name" | cat - files/${sp}_rep_columns_tmp > files/${sp}_repetitive_genomic_regions_annotated
      
        #remove temporal files
       rm files/*tmp
    done
  done   
```

VENN DIAGRAM:

In order to represent the repetitive selection between species in
relation to gene name, I will get the list of genes for each species and
compare those lists:

``` r
library(ggvenn)


# Initialize an empty list to store gene sets
gene_lists <- list()
ensembl_lists<- list()

# Read the gene data sets for each species and store in the list
for (sp in species) {
  df_genes <- read.table(paste0(path, files, sp, "_genomic_regions_annotated"), 
                         sep=" ", header=TRUE)
  ens_ids<- na.omit(df_genes$ensembl_id)
  ensembl_lists[[sp]]<- ens_ids
  genes <- na.omit(df_genes$gene_name)
  gene_lists[[sp]] <- genes
}

# Assign meaningful names to the gene sets
names(gene_lists) <- c("Lynx canadensis", "Lynx lynx", "Lynx pardinus", "Lynx rufus")

# Generate the Venn diagram using ggvenn
p <- ggvenn(
        gene_lists,
        fill_color = c(as.character(colors["lc"]), as.character(colors["ll"]), as.character(colors["lp"]), as.character(colors["lr"])),
        stroke_size = 0,
        show_percentage = FALSE,
        set_name_size = 6,
        text_size = 6,
)

ggsave(filename = paste0(path, plots, "genes_venn_diagram.pdf"), plot= p, width= 15, height = 10)
```
