#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -c 1
#SBATCH --mem=60GB

sp=($(echo $1))

CHR=($(grep -v "#" $LUSTRE/${sp}_goodsamples_cat_ref.filter8.vcf | cut -f1 | grep -vE "(KZ|AANG)" | uniq))

module load cesga/2020
module load samtools
module load gcccore/6.4.0 bcftools/1.9
module load gcccore/system shapeit4/4.2.1

for j in ${CHR[@]}
do
    bgzip -c $LUSTRE/${j}_${sp}_whatshap_phased_duplicated.vcf > $LUSTRE/${j}_${sp}_whatshap_phased_duplicated.vcf.gz
    bcftools index $LUSTRE/${j}_${sp}_whatshap_phased_duplicated.vcf.gz
    shapeit4.2 --input $LUSTRE/${j}_${sp}_whatshap_phased_duplicated.vcf.gz \
                --map $LUSTRE/genetic_maps/chr${j}.gmap \
                --region ${j} \
                --use-PS 0.0001 \
                --output $LUSTRE/${j}_${sp}_shapeit_phased_duplicated_cat_ref.vcf \
                --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
done
