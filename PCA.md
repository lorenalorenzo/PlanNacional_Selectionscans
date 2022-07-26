---
Title: PCA of filtered VCF
Author: Lorena Lorenzo Fernández
Date: 22 July 2021
---

# Once I have completed all of the filtering steps of my VCF, I can proceed to run a Principal Component Analysis,
# in order to explore the overall population structure and visualize my data. For this I will be using plink 1.9 and my final filtered vcf.
# For understanding, see PCA.md on Enrico's Github (https://github.com/Enricobazzi/Selection_Eurasian_Lynx/blob/master/PCA.md)

VCF=/home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter8.vcf

# perform linkage pruning - i.e. identify prune sites
plink_1.9 --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out goodsamples_cat_ref

#run PCA in the subset
plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract goodsamples_cat_ref.prune.in \
--make-bed --pca --out goodsamples_cat_ref

####################Open R and plot####################

#Separate per sp vcf
species=(lc ll lp lr)

for i in ${species[@]}
  do
    #Name per sp samples
    samples=($(grep -m1 "#CHR" /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter8.vcf | tr '\t' '\n' | grep "${i}"))
    /opt/gatk-4.1.0.0/gatk SelectVariants \
      -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
      -V /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter8.vcf \
      $(for j in ${samples[@]}; do echo "-sn ${j}"; done) \
      -O /home/llorenzo/vcf_filtering/${i}_goodsamples_cat_ref.filter8.vcf
  done

#Repeat plink per sp.
for i in ${species[@]}
  do
    VCF=/home/llorenzo/vcf_filtering/${i}_goodsamples_cat_ref.filter8.vcf
    # perform linkage pruning - i.e. identify prune sites
    plink_1.9 --vcf $VCF --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.1 --out ${i}_goodsamples_cat_ref
    #run PCA in the subset
    plink_1.9 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract ${i}_goodsamples_cat_ref.prune.in \
    --make-bed --pca --out ${i}_goodsamples_cat_ref
  done
