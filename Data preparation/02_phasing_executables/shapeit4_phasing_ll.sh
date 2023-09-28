#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -c 1
#SBATCH --mem=60GB


CHR=($(grep -v "#" $LUSTRE/ll_goodsamples_cat_ref.filter8.vcf | cut -f1 | grep -vE "(KZ|AANG)" | uniq))

module load cesga/2020
module load samtools
module load gcccore/6.4.0 bcftools/1.9
module load gcccore/system shapeit4/4.2.1


for i in ${CHR[@]}
do
    bgzip -c $LUSTRE/${i}_ll_prephased_cat_ref.vcf > $LUSTRE/${i}_ll_prephased_cat_ref.vcf.gz
    bcftools index $LUSTRE/${i}_ll_prephased_cat_ref.vcf.gz
    shapeit4.2 --input $LUSTRE/${i}_ll_prephased_cat_ref.vcf.gz \
                --map $LUSTRE/genetic_maps/chr${i}.gmap \
                --region ${i} \
                --use-PS 0.0001 \
                --output $LUSTRE/${i}_ll_shapeit_phased_cat_ref.vcf \
                --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
done
