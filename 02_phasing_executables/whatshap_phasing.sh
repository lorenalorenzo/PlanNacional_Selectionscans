#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

sp=($(echo $1))
chr=($(echo $2))

module load cesga/2020
module load whatshap/1.1


whatshap phase \
  -o $LUSTRE/${chr}_${sp}_whatshap_phased_cat_ref.vcf \
  --reference=$STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
  $LUSTRE/${chr}_${sp}_tmp.vcf \
  $(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/lorena_bams/*${sp}*_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam)
