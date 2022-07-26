#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

module load bedtools

#Now we are going to see the intersect btw sumatran tiger vcf and lynx vcf
echo "Starting lynx_tiger intersect"
bedtools intersect -a $LUSTRE/goodsamples_cat_ref.filter8.vcf -b $LUSTRE/Samaha_et_al_Sumatran_tiger.vcf > $LUSTRE/lynx_tiger_intersect_cat_ref.vcf
echo "ended"

echo "Starting lynx_cheetah intersect"
bedtools intersect -a $LUSTRE/goodsamples_cat_ref.filter8.vcf -b $LUSTRE/Samaha_et_al_cheetah.vcf > $LUSTRE/lynx_cheetah_intersect_cat_ref.vcf
echo "ended"
