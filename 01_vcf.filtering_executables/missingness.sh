######################################################################
######################## Missingness filter ##########################
######################################################################

#Separate per sp vcf
species=(lc ll lp lr)

for i in ${species[@]}
  do
    #Name per sp samples
    samples=($(grep -m1 "#CHR" /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.2.vcf | tr '\t' '\n' | grep "${i}"))
    /opt/gatk-4.1.0.0/gatk SelectVariants \
      -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
      -V /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.2.vcf \
      $(for j in ${samples[@]}; do echo "-sn ${j}"; done) \
      -O /home/llorenzo/vcf_filtering/missing_filter/${i}_goodsamples_cat_ref.filter6.2.vcf
    #Search for variants with >70% of genotypes missing
    bcftools filter -i "F_MISSING>=0.7" -Ov /home/llorenzo/vcf_filtering/missing_filter/${i}_goodsamples_cat_ref.filter6.2.vcf \
          > /home/llorenzo/vcf_filtering/missing_filter/${i}_goodsamples_cat_ref.filter7.1.vcf
    #Number of variants with >70% of genotypes missing
    echo "For ${i} we have $(grep -v "#" /home/llorenzo/vcf_filtering/missing_filter/${i}_goodsamples_cat_ref.filter7.1.vcf| wc -l) variants missing"
  done

#Copy the file we are going to modify with final name
cp /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.2.vcf /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter7.vcf

species=(lc ll lp lr)

for i in ${species[@]}
  do
    echo "subtract missing variants for ${i}"
    #Eliminate those variants
    bedtools subtract -header \
          -a /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter7.vcf \
          -b /home/llorenzo/vcf_filtering/missing_filter/${i}_goodsamples_cat_ref.filter7.1.vcf \
              > /home/llorenzo/vcf_filtering/tmp
      mv /home/llorenzo/vcf_filtering/tmp /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter7.vcf
  done
