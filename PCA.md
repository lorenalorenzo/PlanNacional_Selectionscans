# PCA

------------------------------------------------------------------------

### 

```{bash}
###PCA
#Repeat plink per sp.
species=(lc ll lp lr)

for i in ${species[@]}
  do
    VCF=/home/llorenzo/selection_scan/rehh_test/${i}_goodsamples_filtered_phased_polarized_ACfiltered_header_cat_ref.vcf
    # perform linkage pruning - i.e. identify prune sites
    plink_1.9 --vcf $VCF --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.1 --out ${i}_goodsamples_filtered_phased_polarized_cat_ref
    #run PCA in the subset
    plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract ${i}_goodsamples_filtered_phased_polarized_cat_ref.prune.in \
    --make-bed --pca --out ${i}_goodsamples_filtered_phased_polarized_cat_ref
  done


```

```{r}
library(tidyverse)
library(viridis)


species<- c("lc", "ll", "lp", "lr")
path <- "~/"

##PREPHASED PCAs##
for (sp in species)
{
  #Read pre-phased data
  # import eigen vec and val
pca <- read_table(paste0(path, sp, "_goodsamples_cat_ref.eigenvec"),
                   col_names = FALSE)
eigenval <- scan(paste0(path, sp, "_goodsamples_cat_ref.eigenval"))  

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

print("named cols")

# first convert to percentage variance explained
pve <- data.frame(PC = 1:(ncol(pca)-1), pve = eigenval/sum(eigenval)*100)
# then make a plot
pdf(paste0(path, sp, "prephased_pve.pdf"))
print(ggplot(pve, aes(PC, pve)) +
        geom_bar(stat = "identity") +
        ylab("Percentage variance explained") + 
        theme_light())
dev.off()
print("pve plot done")

# plot PC1 PC2
pdf(paste0(path, sp, "prephased_pca.pdf"))

print(ggplot() +
   geom_point(pca, mapping=aes(x= PC1, y= PC2, color=ind) , size=4, alpha = 0.5) +
    geom_text(pca, mapping=aes(x= PC1, y= PC2, label=ind), size=3, check_overlap = TRUE) +
    coord_equal() +
    theme() +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")))
dev.off()
print("pca plot done")
}

##PHASED PCAs##
for (sp in species)
{
#Read phased data
pca <- read_table(paste0(path, sp, "_goodsamples_filtered_phased_polarized_cat_ref.eigenvec"),
                   col_names = FALSE)
eigenval <- scan(paste0(path, sp, "_goodsamples_filtered_phased_polarized_cat_ref.eigenval"))

print("read data")
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

print("named cols")

# first convert to percentage variance explained
pve <- data.frame(PC = 1:(ncol(pca)-1), pve = eigenval/sum(eigenval)*100)
# then make a plot
pdf(paste0(path, sp, "phased_pve.pdf"))
print(ggplot(pve, aes(PC, pve)) +
        geom_bar(stat = "identity") +
        ylab("Percentage variance explained") + 
        theme_light())
dev.off()
print("pve plot done")

# plot PC1 PC2
pdf(paste0(path, sp, "phased_pca.pdf"))

print(ggplot() +
   geom_point(pca, mapping=aes(x= PC1, y= PC2, color=ind) , size=4, alpha = 0.5) +
    geom_text(pca, mapping=aes(x= PC1, y= PC2, label=ind), size=3, check_overlap = TRUE) +
    coord_equal() +
    theme() +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")))
dev.off()
print("pca plot done")
}
```

```{bash}
bcftools view -e 'INFO/AF=1.00 | INFO/AF=0.00' *_lp_whatshap_phased_duplicated.vcf | wc -l
```

```{bash}

#LOOKING FOR THE ERROR

species=(lc lp lr)
CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))

for sp in ${species[@]}
  do
    for chr in ${CHR[@]:0:19} 
        do
          grep -v "#" ${chr}_${sp}_whatshap_phased_duplicated.vcf  >> ${sp}_tryal_nvariants.vcf
      done
      wc -l ${sp}_tryal_nvariants.vcf
      cat <(grep "#" A1_${sp}_whatshap_phased_duplicated.vcf) ${sp}_tryal_nvariants.vcf > ${sp}_tryal_nvariants_header.vcf    
      bcftools view -e 'INFO/AF=1.00 | INFO/AF=0.00' ${sp}_tryal_nvariants_header.vcf | wc -l
    done
    

for sp in ${species[@]}
  do
    for chr in ${CHR[@]:0:19} 
        do
          grep -v "#" ${chr}_${sp}_shapeit_phased_duplicated_cat_ref.vcf  >> ${sp}_tryal_shapeit_nvariants.vcf
      done
      wc -l ${sp}_tryal_shapeit_nvariants.vcf
      cat <(grep "#" $LUSTRE/phasing/intermediate_files/A1_${sp}_whatshap_phased_duplicated.vcf) ${sp}_tryal_shapeit_nvariants.vcf > ${sp}_tryal_shapeit_nvariants_header.vcf    
      bcftools view -e 'INFO/AF=1.00 | INFO/AF=0.00' ${sp}_tryal_shapeit_nvariants_header.vcf | wc -l
    done  
   
```
