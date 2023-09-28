######################################################################
############### False heterozygosity filter ##########################
######################################################################

# to avoid repeating a SelectVariants to divide the vcf per species,
# we are going to subtract the missingness to the vcf per species from the previous step
species=(lc ll lp lr)

for i in ${species[@]}
  do
    #Eliminate those variants
    bedtools subtract -header \
          -a /home/llorenzo/vcf_filtering/missing_filter/${i}_goodsamples_cat_ref.filter6.2.vcf \
          -b /home/llorenzo/vcf_filtering/missing_filter/${i}_goodsamples_cat_ref.filter7.1.vcf \
              > /home/llorenzo/vcf_filtering/het_filter/${i}_goodsamples_cat_ref.filter7.2.vcf
  done


# with this code we obtain a file where we have the number of het per variant.
for i in ${species[@]}
  do
    paste <(bcftools view /home/llorenzo/vcf_filtering/het_filter/${i}_goodsamples_cat_ref.filter7.2.vcf |\
      awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
        !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
      \
      <(bcftools query -f '[\t%SAMPLE=%GT]\n' /home/llorenzo/vcf_filtering/het_filter/${i}_goodsamples_cat_ref.filter7.2.vcf |\
        awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
        > /home/llorenzo/vcf_filtering/het_filter/${i}_number_het
  done

# now we need to know which variants have MORE THAN 80% OF HET.
awk -F'\t' '{ if ($6 > 15) { print $1, $2-1, $2 } } ' /home/llorenzo/vcf_filtering/het_filter/lc_number_het \
  | tr ' ' '\t' | tail -n +2 \
  > /home/llorenzo/vcf_filtering/het_filter/lc_het_filter8.bed
awk -F'\t' '{ if ($6 > 26) { print $1, $2-1, $2 } } ' /home/llorenzo/vcf_filtering/het_filter/ll_number_het \
  | tr ' ' '\t' | tail -n +2 \
  > /home/llorenzo/vcf_filtering/het_filter/ll_het_filter8.bed
awk -F'\t' '{ if ($6 > 9) { print $1, $2-1, $2 } } ' /home/llorenzo/vcf_filtering/het_filter/lp_number_het \
  | tr ' ' '\t' | tail -n +2 \
  > /home/llorenzo/vcf_filtering/het_filter/lp_het_filter8.bed
awk -F'\t' '{ if ($6 > 14) { print $1, $2-1, $2 } } ' /home/llorenzo/vcf_filtering/het_filter/lr_number_het \
  | tr ' ' '\t' | tail -n +2 \
  > /home/llorenzo/vcf_filtering/het_filter/lr_het_filter8.bed

#Copy the file we are going to modify with final name
cp /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter7.vcf /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter8.vcf

species=(lc ll lp lr)

#Eliminate variants with >80% of het
for i in ${species[@]}
  do
    bedtools subtract -header \
          -a /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter8.vcf \
          -b /home/llorenzo/vcf_filtering/het_filter/${i}_het_filter8.bed \
              > /home/llorenzo/vcf_filtering/tmp
        mv /home/llorenzo/vcf_filtering/tmp /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter8.vcf
  done
