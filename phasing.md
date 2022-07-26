# Title: Inferring haplotypes from genotypes: Phasing data

## Author: Lorena Lorenzo Fernández

## Date: 11 March, 2021

Here I present each step process I run trough while achieving my purpose. [The complete and summarized code for phasing is available in phasing.sh.]{.underline}

We are going to use whatshap for blocking and shapeit for haplotype phasing (with per species vcf)

``` bash
#Whatshap common use:

whatshap phase --tag=PS -o ${i}_prephased_cat_ref.vcf --reference=reference.fasta ${i}_goodsamples_cat_ref.filter8.vcf input.bam

#First we need to change the lr_0011 sample name in order to match with vcf header (bf using whatshap)
module load picard

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
            I=c_lr_xx_0011_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
            O=c_lc_xx_0011_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
            RGID=L1-Lynx_S9 \
            RGLB=c_lc_xx_0011_lib \
            RGPL=Illumina \
            RGSM=c_lc_xx_0011 \
            RGPU=05

#And then add an index to the new bam file:

java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
      I=c_lc_xx_0011_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
```

Once did this, I tried to do the phasing in lp_vcf

``` bash
module load cesga/2020
module load whatshap/1.1
#not needed yet: module load shapeit4/4.2.1
species=(lc ll lp lr)

for i in ${species[@]}
do
whatshap phase
  -o $LUSTRE/lp_prephased_cat_ref.vcf \
  --tag=PS \
  --reference=$STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
  $LUSTRE/lp_goodsamples_cat_ref.filter8.vcf \
  $STORE2/lynx_genome/lynx_data/CatRef_bams/*lp*_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
```

This was so time consuming (2,5 days). See phasing_lp.sh. Because of that we decided to do it per chr (See per_chr_sp_vcf.sh).

    sp=(lc ll lp lr)

    for i in ${sp[@]}
    do
      echo "$i"
      CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
      for j in ${CHR[@]:0:20}
        do
          echo "$j"
          grep -E "^(#|${j})"  ${i}_goodsamples_cat_ref.filter8.vcf > ${j}_${i}_tmp.vcf
          echo "Done $j for $i vcf"
        done
    done

So, once we have a per chr_sp_vcf, we pre-phase it.

    sp=(lc ll lp lr)

    for i in ${sp[@]}
    do
      echo "$i"
      CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
      for j in ${CHR[@]:0:20}
        do
        echo "Work for $i in $j"
        sbatch whatshap_phasing.sh $i $j
        done
    done

Having done the pre-phasing step with Whatshap, we are now going to phase with shapeit4. Following the manual, we need to add a **genetic map**. If we don't have a genetic map, by default the program uses a 1cM per 1Mb of recombination rate. However, we have that the average recombination rate in cats is 1.9cM. Talking with the author of shapeit4, the only way to apply our calculated 1.9cM rate is making a genetic map with this data. So we will need a per chr map with information about variants (pos), chromosome (chr) and cM.

``` bash
#BUILDING A GENETIC MAP WITH CONSTANT 1.9CM
  CHR=($(grep -v "#" goodsamples_cat_ref.filter8.vcf | cut -f1 | grep -vE "(KZ|AANG)" | uniq))

    for i in ${CHR[@]}
    do
      grep -v "#" goodsamples_cat_ref.filter8.vcf | grep ${i} | cut -f1,2 | \
      awk '{ print $2, $1 }' | \
      awk {'if ( NR==1 ) print $1, $2, 0; else print $1, $2, $1-p, ($1-p)*0.0000019; p=$1'} | \
      awk 'BEGIN{print "pos", "chr", "cM"} {sum+=$4} {print $1, $2, sum}' | tr ' ' '\t' > /home/llorenzo/chr${i}.gmap
    done
```

With this code we are getting a genetic map with our variant information (from our filtered vcf). As a little explanation we have taken the variant and chromosome information from the vcf. Then, we have calculated the distance between variants: position (n) - position (n-1). With the distance, we multiplicated to 1.9cM and divide to 1Mb, which is the same as multiplicating to 0.0000019. Finally we did a cumulative sum of the cM calculated and add the header.

The next step is phasing data with shapeit4.

``` bash
    module load cesga/2020
    module load gcccore/system shapeit4/4.2.1

    shapeit4.2 --input $LUSTRE/A1_lc_prephased_cat_ref.vcf.gz --map $LUSTRE/genetic_maps/chrA1.gmap --region A1 --use-PS 0.0001 --output $LUSTRE/A1_lc_phased_cat_ref.vcf
```

The problem found by this point is that the program needs more than 20 individuals to infer phase. So, we decided to investigate if we can use long reads from PN2017. Finally, this wasn't an option because the individuals at long reads are not the same I have in my dataset.

Meanwhile, I am going to try to phase ll (lynx lynx) with shapeit4 (See shapeit4_phasing.sh) and the rest with whatshap (See whatshap_phasing.sh) without the long reads while waiting email response

``` bash
sp=(lc ll lp lr)

for i in ${sp[@]}
do
  echo "$i"
  CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
  for j in ${CHR[@]:0:20}
    do
    echo "Work for $i in $j"
    sbatch whatshap_phasing.sh $i $j
    done
done
```

Finally, comparing Whatshap and shapeit4 results for lynx lynx I found that Whatshap only infers phase when is sure about it (leaving aprox 50% of the data unphased). For that reason, and following the advice of the whatshap developer, I am going to duplicate data to make shapeit4 work. I decided to duplicate the whatshap_prephased results, after being sure that duplication don't affect whasthap prephasing results.

``` bash
CHR=($(grep -v "#" $LUSTRE/ll_goodsamples_cat_ref.filter8.vcf | cut -f1 | grep -vE "(KZ|AANG)" | uniq))

for i in ${CHR[@]}
do
echo "Duplicating data for $i"
sbatch duplicating_for_phasing.sh $i
done
```

Once duplicated, we are going to test if shapeit4 works.

``` bash
sp=(lc lp lr)

for i in ${sp[@]}
do
  echo "Shapeit4 for $i"
  sbatch shapeit4_phasing_duplicates.sh $i
done  
```

After that, what we want is a per sp vcf and eliminate duplicates

```{bash}
for chr in ${CHR[@]:0:20} 
  do
    bgzip -c $LUSTRE/phasing/${chr}_ll_prephased_cat_ref.vcf > $LUSTRE/phasing/${chr}_ll_prephased_cat_ref.vcf.gz
    bcftools index $LUSTRE/phasing/${chr}_ll_prephased_cat_ref.vcf.gz
    
    shapeit4.2 \
      --input $LUSTRE/phasing/${chr}_ll_prephased_cat_ref.vcf.gz \
      --map $LUSTRE/genetic_maps/chr${chr}.gmap \
      --region ${chr} \
      --use-PS 0.0001 \
      --output $LUSTRE/phasing/${chr}_ll_shapeit_phased_cat_ref.vcf \
      --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
  done
```
