#################################################################
######## Extract max_depth non callable variants from vcf #######
#################################################################

#copy the file we want to modify and put the final file name we want.
cp /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter5.vcf /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.vcf

# Subtract non-callable variants from allsamples vcf
species=(lc ll lp lr)


for i in ${species[@]}
 do
  echo ${i}
  bedtools subtract -header \
      -a /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.vcf \
      -b /home/llorenzo/vcf_filtering/depth_filter/${i}_goodsamples_cat_ref_noncallable_regions.vcf \
      > /home/llorenzo/vcf_filtering/tmp
  mv /home/llorenzo/vcf_filtering/tmp /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.vcf
 done
