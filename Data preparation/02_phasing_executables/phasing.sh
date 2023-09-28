#!/bin/bash

#Dependencies needed
module load cesga/2020
module load samtools
module load whatshap/1.1
module load gcccore/system shapeit4/4.2.1

#Define variables
CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
species=(lc ll lp lr)

#First step: pre-phase with Whatshap
for sp in ${species[@]}
  do
    echo "$sp"
    for chr in ${CHR[@]:0:20} 
      do
        echo "Phasing $chr in $sp"
        whatshap phase \
          -o $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_cat_ref.vcf \
          --tag=PS \
          --reference=$STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
          $LUSTRE/phasing/${chr}_${sp}_tmp.vcf \
          $(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/lorena_bams/*${sp}*_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam)
      done
  done
  
#Duplicate lc,lr and lp (required for using shapeit4)
for chr in ${CHR[@]:0:20} 
  do
    grep -v "#" $LUSTRE/phasing/${chr}_lp_whatshap_prephased_cat_ref.vcf  | rev | cut -f 1-11 | rev > $LUSTRE/phasing/${chr}_lp_GT.vcf
    grep -v "#" $LUSTRE/phasing/${chr}_lc_whatshap_prephased_cat_ref.vcf  | rev | cut -f 1-19 | rev > $LUSTRE/phasing/${chr}_lc_GT.vcf
    grep -v "#" $LUSTRE/phasing/${chr}_lr_whatshap_prephased_cat_ref.vcf  | rev | cut -f 1-18 | rev > $LUSTRE/phasing/${chr}_lr_GT.vcf
    paste <(grep "#CHR" $LUSTRE/phasing/${chr}_lp_whatshap_prephased_cat_ref.vcf ) <(paste <(grep "#CHR" $LUSTRE/phasing/${chr}_lp_whatshap_prephased_cat_ref.vcf  \
    | tr '\t' '\n' | grep "_") <(yes "2" | head -11) | tr '\t' '_' | tr '\n' '\t') | cut -f 1-31 > $LUSTRE/phasing/${chr}_lp_samplesnames.vcf
    paste <(grep "#CHR" $LUSTRE/phasing/${chr}_lc_whatshap_prephased_cat_ref.vcf ) <(paste <(grep "#CHR" $LUSTRE/phasing/${chr}_lc_whatshap_prephased_cat_ref.vcf  \
    | tr '\t' '\n' | grep "_") <(yes "2" | head -19) | tr '\t' '_' | tr '\n' '\t') | cut -f 1-47 > $LUSTRE/phasing/${chr}_lc_samplesnames.vcf
    paste <(grep "#CHR" $LUSTRE/phasing/${chr}_lr_whatshap_prephased_cat_ref.vcf ) <(paste <(grep "#CHR" $LUSTRE/phasing/${chr}_lr_whatshap_prephased_cat_ref.vcf  \
    | tr '\t' '\n' | grep "_") <(yes "2" | head -18) | tr '\t' '_' | tr '\n' '\t') | cut -f 1-45 > $LUSTRE/phasing/${chr}_lr_samplesnames.vcf
    
    for sp in lc lp lr
      do
        grep "##" $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_cat_ref.vcf  > $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_duplicated.vcf
        grep -v "#" $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_cat_ref.vcf  | cut -f 1-9 > $LUSTRE/phasing/${chr}_${sp}_info.vcf
        cat $LUSTRE/phasing/${chr}_${sp}_samplesnames.vcf >> $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_duplicated.vcf
        paste $LUSTRE/phasing/${chr}_${sp}_info.vcf $LUSTRE/phasing/${chr}_${sp}_GT.vcf $LUSTRE/phasing/${chr}_${sp}_GT.vcf >> $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_duplicated.vcf
    
        rm $LUSTRE/phasing/${chr}_${sp}_GT.vcf
        rm $LUSTRE/phasing/${chr}_${sp}_samplesnames.vcf
        rm $LUSTRE/phasing/${chr}_${sp}_info.vcf
      done
  done
  
#Second step: phase data with Shapeit4
for chr in ${CHR[@]:0:20} 
  do
    for sp in lc lp lr
      do  
        bgzip -c $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_duplicated.vcf > $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_duplicated.vcf.gz
        
        module load cesga/2018  gcccore/6.4.0
        module load gcccore/6.4.0 bcftools/1.9
        bcftools index $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_duplicated.vcf.gz
        
        module load cesga/2020
        module load gcccore/system shapeit4/4.2.1
        shapeit4.2 \
          --input $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_duplicated.vcf.gz \
          --map $LUSTRE/genetic_maps/chr${chr}.gmap \
          --region ${chr} \
          --use-PS 0.0001 \
          --output $LUSTRE/phasing/${chr}_${sp}_shapeit_phased_duplicated_cat_ref.vcf \
          --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
      done
      
      sp=ll
        bgzip -c $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_cat_ref.vcf > $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_cat_ref.vcf.gz
        
        module load cesga/2018  gcccore/6.4.0
        module load gcccore/6.4.0 bcftools/1.9
        bcftools index $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_cat_ref.vcf.gz
    
        module load cesga/2020
        module load gcccore/system shapeit4/4.2.1
        shapeit4.2 \
          --input $LUSTRE/phasing/${chr}_${sp}_whatshap_prephased_cat_ref.vcf.gz \
          --map $LUSTRE/genetic_maps/chr${chr}.gmap \
          --region ${chr} \
          --use-PS 0.0001 \
          --output $LUSTRE/phasing/${chr}_${sp}_shapeit_phased_cat_ref.vcf \
          --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m

  done

#Eliminate duplicated data and get a vcf per sp
#Vcf per species WITHOUT header
for chr in ${CHR[@]:0:19}
  do 
    grep -v "#" $LUSTRE/phasing/phased_vcf/${chr}_lc_shapeit_phased_duplicated_cat_ref.vcf | cut -f 1-28 >> lc_shapeit_tmp.vcf
    grep -v "#" $LUSTRE/phasing/phased_vcf/${chr}_lp_shapeit_phased_duplicated_cat_ref.vcf | cut -f 1-20 >> lp_shapeit_tmp.vcf
    grep -v "#" $LUSTRE/phasing/phased_vcf/${chr}_lr_shapeit_phased_duplicated_cat_ref.vcf | cut -f 1-27 >> lr_shapeit_tmp.vcf
    grep -v "#" $LUSTRE/phasing/phased_vcf/${chr}_ll_shapeit_phased_cat_ref.vcf >> ll_shapeit_tmp.vcf
  done
#Add a header
cat <(grep -E "^(##fileformat|##INFO|##FORMAT|#CHR)" A1_lc_shapeit_phased_duplicated_cat_ref.vcf | cut -f 1-28 ) lc_shapeit_tmp.vcf > lc_shapeit_phased_cat_ref.vcf
cat <(grep -E "^(##fileformat|##INFO|##FORMAT|#CHR)" A1_lp_shapeit_phased_duplicated_cat_ref.vcf | cut -f 1-20 ) lp_shapeit_tmp.vcf > lp_shapeit_phased_cat_ref.vcf
cat <(grep -E "^(##fileformat|##INFO|##FORMAT|#CHR)" A1_lr_shapeit_phased_duplicated_cat_ref.vcf | cut -f 1-27 ) lr_shapeit_tmp.vcf > lr_shapeit_phased_cat_ref.vcf
cat <(grep -E "^(##fileformat|##INFO|##FORMAT|#CHR)" A1_ll_shapeit_phased_cat_ref.vcf) ll_shapeit_tmp.vcf > ll_shapeit_phased_cat_ref.vcf

#Merge sp vcf in only one vcf (note that previously I tested that each vcf had the same number of variants and remember I eliminated variants from MT and scaffolds KZ/AANG )
for sp in ${species[@]}
  do
    bgzip -c $LUSTRE/phasing/phased_vcf/${sp}_shapeit_phased_cat_ref.vcf > ${sp}_shapeit_phased_cat_ref.vcf.gz
    bcftools index $LUSTRE/phasing/phased_vcf/${sp}_shapeit_phased_cat_ref.vcf.gz
  done 
  
bcftools merge lc_shapeit_phased_cat_ref.vcf.gz ll_shapeit_phased_cat_ref.vcf.gz lp_shapeit_phased_cat_ref.vcf.gz lr_shapeit_phased_cat_ref.vcf.gz  -o goodsamples_phased_cat_ref.vcf.gz
    