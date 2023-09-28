##########################################
######## Separate a VCF by species #######
##########################################

screen -S vcf_per_species

#First we are going to see the vcf header
grep -m1 "#CHR" allsamples_cat_ref.filter5.vcf

#Rename wrong sample using bcftools reheader.
bcftools reheader -s <(echo 'c_lr_xx_0011 c_lc_xx_0011') /home/llorenzo/vcf/allsamples_cat_ref.filter5.vcf -o /home/llorenzo/vcf/renamed_allsamples_cat_ref.filter5.vcf

#Eliminate BAD_SAMPLES (only c_lc_yu_0007)
/opt/gatk-4.1.0.0/gatk SelectVariants \
  -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
  -V /home/llorenzo/vcf/renamed_allsamples_cat_ref.filter5.vcf \
  -xl-sn "c_lc_yu_0007" \
  -O /home/llorenzo/vcf/goodsamples_cat_ref.filter5.vcf

#Separate per sp vcf
species=(lc ll lp lr)

for i in ${species[@]}
  do
    samples=($(grep -m1 "#CHR" /home/llorenzo/vcf/goodsamples_cat_ref.filter5.vcf | tr '\t' '\n' | grep "${i}"))
    /opt/gatk-4.1.0.0/gatk SelectVariants \
      -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
      -V /home/llorenzo/vcf/goodsamples_cat_ref.filter5.vcf \
      $(for j in ${samples[@]}; do echo "-sn ${j}"; done) \
      -O /home/llorenzo/vcf/${i}_goodsamples_cat_ref_filter5.vcf
  done
