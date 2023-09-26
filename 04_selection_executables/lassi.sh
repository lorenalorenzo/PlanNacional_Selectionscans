#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=10GB

#name variables
sp=($(echo $1))
CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))

#load dependencies
module load lassip

#run lassi
for chr in ${CHR[@]:0:18} 
      do
        echo "trying salti-lassi in $chr"    
        lassip \
          --vcf $LUSTRE/selection_scan/chr_files/${chr}_${sp}_goodsamples_filtered_polarized_variants_header_cat_ref.vcf \
          --unphased \
          --calc-spec \
          --hapstats \
          --salti \
          --winsize 101 \
          --winstep 50 \
          --pop $LUSTRE/selection_scan/${sp}_ind.txt \
          --out $LUSTRE/selection_scan/saltiLASSI/${chr}
      done

lassip \
       --spectra $LUSTRE/selection_scan/saltiLASSI/*${sp}.lassip.hap.spectra.gz \
       --salti \
       --out $LUSTRE/selection_scan/saltiLASSI/${sp}_salti
