#!/bin/bash


for sp in lc ll lp lr
 do
 
  # generate a CSV table with the number of missing genotypes at each SNP for each sp
    echo "chromosome, position, n_missing" > $LUSTRE/vcfs/${sp}_missing_gts.csv

  # chromosome and position from any of the VCFs
    grep -v "#" $LUSTRE/vcfs/${sp}_goodsamples_cat_ref.filter8.vcf | cut -f1-2 > $LUSTRE/vcfs/${sp}_chr_pos.tmp

  # adding missing gts to table
    echo "calculate missing gts in ${sp}"
    nsamples=($(grep -m1 "#CHR" $LUSTRE/vcfs/${sp}_goodsamples_cat_ref.filter8.vcf | tr '\t' '\n' | grep "_" | wc -l))
    grep -v "#" $LUSTRE/vcfs/${sp}_goodsamples_cat_ref.filter8.vcf | cut -f8 | cut -d';' -f3 | 
     cut -d'=' -f2 | awk -v nsam="${nsamples}" '{print ((nsam*2)-$1)/2}' > $LUSTRE/vcfs/${sp}_miss.tmp
    
    echo "adding $sp missing gts to table" 
    paste $LUSTRE/vcfs/${sp}_chr_pos.tmp $LUSTRE/vcfs/${sp}_miss.tmp > tmp && mv tmp ${sp}_chr_pos.tmp
  
  # combine table
    cat $LUSTRE/vcfs/${sp}_chr_pos.tmp | tr '\t' ',' >> $LUSTRE/vcfs/${sp}_missing_gts.csv
  
  done

# remove tmp files
rm *.tmp

#extract the proportion of missing data per sample
paste \
<(bcftools query -f '[%SAMPLE\t]\n' ll_goodsamples_cat_ref.filter8.vcf | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' ll_goodsamples_cat_ref.filter8.vcf | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)
