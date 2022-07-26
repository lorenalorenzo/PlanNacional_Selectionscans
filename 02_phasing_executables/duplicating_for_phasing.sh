#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -c 1
#SBATCH --mem=20GB

#We cannot use shapeit4 bc of poor number of samples so we are trying to use duplications. First we need to separate GT information from Header and the rest of the information

chr=($(echo $1))

  grep -v "#" $LUSTRE/${chr}_lp_whatshap_phased_cat_ref.vcf  | rev | cut -f 1-11 | rev > $LUSTRE/${chr}_lp_GT.vcf
  grep -v "#" $LUSTRE/${chr}_lc_whatshap_phased_cat_ref.vcf  | rev | cut -f 1-19 | rev > $LUSTRE/${chr}_lc_GT.vcf
  grep -v "#" $LUSTRE/${chr}_lr_whatshap_phased_cat_ref.vcf  | rev | cut -f 1-18 | rev > $LUSTRE/${chr}_lr_GT.vcf
  paste <(grep "#CHR" $LUSTRE/${chr}_lp_whatshap_phased_cat_ref.vcf ) <(paste <(grep "#CHR" $LUSTRE/${chr}_lp_whatshap_phased_cat_ref.vcf  \
  | tr '\t' '\n' | grep "_") <(yes "2" | head -11) | tr '\t' '_' | tr '\n' '\t') | cut -f 1-31 > $LUSTRE/${chr}_lp_samplesnames.vcf
  paste <(grep "#CHR" $LUSTRE/${chr}_lc_whatshap_phased_cat_ref.vcf ) <(paste <(grep "#CHR" $LUSTRE/${chr}_lc_whatshap_phased_cat_ref.vcf  \
  | tr '\t' '\n' | grep "_") <(yes "2" | head -19) | tr '\t' '_' | tr '\n' '\t') | cut -f 1-47 > $LUSTRE/${chr}_lc_samplesnames.vcf
  paste <(grep "#CHR" $LUSTRE/${chr}_lr_whatshap_phased_cat_ref.vcf ) <(paste <(grep "#CHR" $LUSTRE/${chr}_lr_whatshap_phased_cat_ref.vcf  \
  | tr '\t' '\n' | grep "_") <(yes "2" | head -18) | tr '\t' '_' | tr '\n' '\t') | cut -f 1-45 > $LUSTRE/${chr}_lr_samplesnames.vcf

  sp=(lc ll lp lr)

  for j in ${sp[@]}
  do
    grep "##" $LUSTRE/${chr}_${j}_whatshap_phased_cat_ref.vcf  > $LUSTRE/${chr}_${j}_whatshap_phased_duplicated.vcf
    grep -v "#" $LUSTRE/${chr}_${j}_whatshap_phased_cat_ref.vcf  | cut -f 1-9 > $LUSTRE/${chr}_${j}_info.vcf
    cat $LUSTRE/${chr}_${j}_samplesnames.vcf >> $LUSTRE/${chr}_${j}_whatshap_phased_duplicated.vcf
    paste $LUSTRE/${chr}_${j}_info.vcf $LUSTRE/${chr}_${j}_GT.vcf $LUSTRE/${chr}_${j}_GT.vcf >> $LUSTRE/${chr}_${j}_whatshap_phased_duplicated.vcf
  done
