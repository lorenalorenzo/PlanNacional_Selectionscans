#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --output=lassi.out

#load dependencies
module load lassip

#name variables
species=(lc ll lp lr)
CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))

#make a sp_ind.txt (pop file) needed for lassip
for sp in ${species[@]}
  do
  sed "s/$/\t${sp}/" <(grep "#CHR" $LUSTRE/${sp}_goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf | tr "\t" "\n" | grep "c_")  \
    > $LUSTRE/selection_scan/${sp}_ind.txt
  done

#run lassi
for sp in ${species[@]}
  do
    echo "$sp"
    for chr in ${CHR[@]:0:19} 
      do
        echo "trying salti-lassi in $chr"    
        lassip \
          --vcf $LUSTRE/selection_scan/${chr}_${sp}_goodsamples_filtered_phased_polarized_variants_header_cat_ref.vcf \
          --calc-spec \
          --hapstats \
          --salti \
          --winsize 101 \
          --winstep 10 \
          --pop $LUSTRE/selection_scan/${sp}_ind.txt \
          --out $LUSTRE/selection_scan/${chr}
      done
       lassip \
       --spectra $LUSTRE/selection_scan/*${sp}.lassip.hap.spectra.gz \
       --salti \
       --out $LUSTRE/selection_scan/${sp}_salti
  done