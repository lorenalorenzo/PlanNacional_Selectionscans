# Title: Inferring ancestral state: Polarizing data

## Author: Lorena Lorenzo Fernández

## Date: 18 October, 2021

Some selection analysis need ancestral and derived state (polarized data) so we explored the different methods to infer polarization.

The simplest way to do it is assuming that the reference (in our case Cat) is the ancestral state for each variant. However, this is a big assumption, and the polarization results will be clumsy. So, in order to obtain better results, we are going to try inferring ancestral state with 3 outgroups: puma, tiger and cat (finally lion, serval and cat).

Before than deepen in parsimony issues, we need a synteny between our outgroups in order to infer ancestry. That means that for each variant position in our vcf, we need the equivalency in CAT,TIGER and PUMA. As we have lynxes map to cat, we have the cat synteny (col REF in our vcf) but not for the rest of sp.

In order to get synteny btw this sp. we have explored several ways:

1.  A **synteny "calculator" program**. There are several programs that explore synteny but only in terms of genes, so this way is discarded.

2.  Use **synteny calculated by F.Abascal** for the lynx pardinus genome project. We cannot use this data because it has been calculated for previous versions of Cat reference genome and doesn't include puma.

3.  **Align reference genomes** of puma-cat and tiger-cat and get the equivalency of positions. To align we tried mashmap with each fasta (reference and query). The output had information about cat chromosome, cat fragment size, cat start and end position and then the Puma (or tiger) chr, size, start-end position and percentage of identity (and more info). We pretended to use this in order to exactly know what position we need to look for in the fasta of each sp. E.g: In our vcf we have that in chrA2 (which is referred to Cat chrA2) in position 50000 the REF is A and the ALT is C. So, if it works as we want, mashmap would say that position 50000 of ChrA2 in cat is the same as position 400000 in chrB1 in Tiger, and exactly position in the fasta file is a T. So, we would have that LYNX (ALT) have a C, CAT have an A (REF) and TIGER have a T. Nonetheless it doesn't work this way and the alignment don't count the gaps, so we can't extract the exact base position by this way.

4.  Use the **liftOver GATK tool**. The liftOver command allows you to change the reference to a vcf. It is useful for versions change, so that you can change your data from Cat_4.0 to Cat_9.0 for example. In our case, we want to change the reference from one sp to another, so we have to be careful and prove that this is assumable. The liftOver tool needs a chain format file which is similar to mashmap output but containing gap info. Unfortunately, this chain info is available for model organism but not for tiger neither puma, so we have to generate it. Looking for this, we found a program called Crossmap that does something similar to liftOver but it has in his workflow explained the way to make your own chain file. Presumably we have to align our data (nucmer), then filter it (delta-filter) and then change to chain format (crossmap-delta-to-chain). Once the chain is generated, we will proceed with liftOver (as CrossMap is way more confusing). NOTE: For testing the reliability of this, we are also going to generate a chain file for Cat to lynx canadensis reference and compare the vcf output in comparison to the lynx lynx assembly mapped to lynx canadensis reference. Workflow recommended from crossmap Github: run nucmer (See nucmer_alignment.sh) on reference and query assemblies to align, then filter using the MUMmer utility delta-filter (recommend parameters for a lenient mapping between genome assemblies: -i 99.5 -l 1000).

5.  Use **vcf generated for Samatha et al. 2021** and intersect with mi vcf. By intersecting the vcf from the paper with my vcf I obtain 341,342 variants shared with tiger and 1,068,771 variants shared with cheetah. If we have into consideration that we had 20,578,917 SNPs, over than 19million variants are not covered by this method and would be erroneous to say they are ancestral.

6.  **Align samples from 2 outgroup sp. (serval and lion) to Cat v.9 ref. genome** and take the variant position (synteny). As in Li et al 2019 they used several felid sp aligned, we asked Murphy for the data. Raw data (FASTQ) is available in NCBI, so I had to download it (See downloading_sra.sh). Once downloaded, I map the fastq data to the reference (See serval/lion_CatRef_mapping.sh)

So, after all we used the option 6.

**LION-SERVAL MAPPING:**

I downloaded raw data for lion and serval samples, available in NCBI: <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR836361>

<https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6071636>

Once I have the R1 and R2 fastq for each one, it was mapped to Felis_Catus v.9 reference genome (as it is the one used for my vcf). With this steps, I got a bam file for each species. Now, with the bam file I pretend to get a fasta, with the base called in each position. For this, I am going to use mpileup and pu2fa (available in genomics-b server), as used before by the group (<https://github.com/mlucenaperez/contemporary_analysis/blob/master/*rufus_ancestral_genome.Rmd>).

So, the first thing I need is to move the bam files from CESGA to genomics-b server.

```{bash}
#copy lion bam to genomic server
scp csbiellf@ft2.cesga.es:/mnt/netapp1/Store_csebdjgl/reference_genomes/lion_serval_polarization/SRR836361/lion_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.ba* /home/llorenzo

#copy serval bam to genomic server
scp csbiellf@ft2.cesga.es:/mnt/netapp1/Store_csebdjgl/reference_genomes/lion_serval_polarization/SRR6071636/serval_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.ba* /home/llorenzo
```

Then, I will run samtools mpileup and the output is used by pu2fa to get the final fasta file. Finally, I get the intersect between the fasta and the variants of the vcf, getting the base sinteny in lion and serval.

```{bash}
#Run in a screen in genomics-b server
#Define general paths
REF=/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
CHR=($(cat /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
VCF=/home/llorenzo/goodsamples_cat_ref.filter8.vcf

############################ LION ############################
#Define paths
INBAM=/home/llorenzo/lion_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
OUTFASTA=/home/llorenzo/lion_cat_ref_sorted_rg_rmdup_sorted_indelrealigner_pu2fa.fa
INTERSECTFASTA=/home/llorenzo/lion_cat_ref_intersect_vcf.fa
INTERSECTBED=/home/llorenzo/lion_cat_ref_intersect_vcf.bed

#getting fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}
do
echo "$i"
samtools mpileup -s -q30 -f $REF $INBAM -r $i | /home/GRUPOS/grupolince/reference_genomes/Chrom-Compare-master/pu2fa -c $i -C 100 >> $OUTFASTA
done
#We were getting lots of N, the problem was that I didn't define a maximum coverage and by default is 10 (really low) so we put 100

#Intersect fasta with vcf
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA #we only get 222076 N
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab
############################ SERVAL ############################
#Define paths
INBAM=/home/llorenzo/serval_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
OUTFASTA=/home/llorenzo/serval_cat_ref_sorted_rg_rmdup_sorted_indelrealigner_pu2fa.fa
INTERSECTFASTA=/home/llorenzo/serval_cat_ref_intersect_vcf.fa
INTERSECTBED=/home/llorenzo/serval_cat_ref_intersect_vcf.bed

#getting fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}
do
echo "$i"
samtools mpileup -s -q30 -f $REF $INBAM -r $i | /home/GRUPOS/grupolince/reference_genomes/Chrom-Compare-master/pu2fa -c $i -C 100 >> $OUTFASTA
done
#The same problem than before (max cov)

#Intersect fasta with vcf
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA  #we get only 136599 N
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab

#####################BED SINTENY#########################
cut -f2 serval_cat_ref_intersect_vcf.bed  | paste lion_cat_ref_intersect_vcf.bed - | awk -F"\t|:|-" '{printf ("%s\t%s\t%s\t%s=%s%s\n", $1,$2,$3,"LSC",$4,$5)}' - > lion_serval_sinteny_intersect_vcf.bed

grep -v "#" goodsamples_cat_ref.filter8.vcf | cut -f4 | paste -d'\0' lion_serval_sinteny_intersect_vcf.bed - > lion_serval_cat_sinteny_intersect_vcf.bed
```

**OUTGROUP PARSIMONY**

Once we have synteny inferred btw our outgroup sp. we are going to apply the rules following Maddison et al. 1984. Our outgroup sp. are lion, serval and cat, so the options are: (AAA, AAB, ABA, BAA). Using Dani's code:

```{bash}
#Next, apply parsimony criteria to infer the ancestral state, and print the scaffold, position, new ancestral state and previous ancestral state (i.e. the Lynx rufus base). The logic behind this code is extracted from the results for the previous section. The Iberian and Eurasian lynx base is ignored here (because all sites in my VCF are polymorphic in one or both species, or alternative substitutions).
awk '{                                       
split($0,a,":");
split(a[1],b,"=");
split(b[2],c,"");
if (c[1]==c[2] && c[1]==c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[3]); #LSC* all equal
else if (c[1]==c[2] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[3]); #LS vs C null
else if (c[1]==c[3] && c[1]!=c[2]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[3]); #LC* vs S
else if (c[2]==c[3] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[3]); #SC* vs L

else if ((c[1]=="N") && c[2]==c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[3]); #SC* all equal
else if ((c[1]=="N") && c[2]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[3]); #S vs C null

else printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[3]); #others null
}' lion_serval_cat_sinteny_intersect_vcf.bed > outgroup_parsimony_ancestral_state_lion_serval_cat_variants.bed
```

From the 20,578,917 variants in the vcf 19,463,799 remain consistent (94,6%), which basically means that cat reference is the ancestral allele. So, less than a million variants are inconsistent (we cannot resolve ancestry).

```{bash}
#Inconsistent sites (nulls)
awk '$5!=$6 {print $0}' outgroup_parsimony_ancestral_state_lion_serval_cat_variants.bed > outgroup_parsimony_ancestral_state_lion_serval_cat_null_variants.bed

#Eliminate unpolarisable variants from the vcf (general and per sp)
bedtools subtract \
  -a phasing/goodsamples_phased_cat_ref.vcf \
  -b outgroup_parsimony_ancestral_state_lion_serval_cat_null_variants.bed \
  > goodsamples_filtered_phased_polarized_cat_ref.vcf

for sp in ${species[@]}
  do 
bedtools subtract \
  -a phasing/${sp}_phased_cat_ref.vcf \
  -b outgroup_parsimony_ancestral_state_lion_serval_cat_null_variants.bed \
  > goodsamples_filtered_phased_polarized_cat_ref.vcf

#Test if there are variants with 0 AC...
grep "AC=0" goodsamples_phased_cat_ref.vcf | wc -l
#2.691.650--> FALTA METER LL (reducirá este número...deberia quedaer entorno a 50K)
#51215 (UNA VEZ INCORPORADO LL)
#47687 if we test it in the goodsamples_filtered_phased_polarized_cat_ref.vcf file

#Eliminate variants with 0 AC
grep -v "AC=0" goodsamples_filtered_phased_polarized_cat_ref.vcf > goodsamples_filtered_phased_polarized_ACfiltered_cat_ref.vcf

#Add a header
```
