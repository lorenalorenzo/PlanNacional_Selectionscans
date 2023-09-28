# Title: Calculating per sp. depth and Applying depth filter to VCF

## Author: Lorena Lorenzo Fernández

## Date: 19 April 2021

### 1. Analyzing depth distribution in the genome

The whole genome is too big to analyze by region depth, that's why we only need a representative subset of data. For doing this, I am going to subsample BAM files, previously generating a BED file with 100 random positions of 100000 bp length. Thus, depth calculations will be done only considering this random subset.

To generate the random regions file I will use BEDtools random. I'll then remove the low-mappability and repetitive regions from the file using BEDtools subtract

``` bash
# Create a Genome region file for Bedtools:
# A file with the list of chromosomes as col1 and their length as col2, tab separated
# Basically the first two columns of a FAI file:
# Using an fai index file in conjunction with a FASTA/FASTQ file containing reference sequences enables efficient access to arbitrary regions within those reference sequences. The index file typically has the same filename as the corresponding FASTA/FASTQ file, with .fai appended.

cut -f1,2 /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai > /home/llorenzo/depth_filter/Felis_catus_9.0.dna.toplevel.genome

# Bedtools random to generate file of 100 random segments of 100000 bp
# Output a BED file:

bedtools random -l 100000 -n 100 -g /home/llorenzo/depth_filter/Felis_catus_9.0.dna.toplevel.genome | sort > /home/llorenzo/depth_filter/Felis_catus.100x100kbp.genome.bed

# Using bedtools subtract I can remove low-mappability and repetitive regions:

bedtools subtract -a /home/llorenzo/depth_filter/Felis_catus.100x100kbp.genome.bed -b /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Masked_Regions.bed > /home/llorenzo/depth_filter/Felis_catus.100x100kbp.masked.genome.bed
```

Then, I have to calculate region coverage per sample with samtools depth (loop to do this in every sample). Depth at all positions will be calculated (-a) within the regions randomly selected before (-b).

This will be run in a loop for each lp bam sample.

    # Loop of Samtools depth calculations for each sample_bam
    for i in *lp*.bam
      do
      echo "Calculating depth for $i"
      samtools depth -a -b /home/llorenzo/depth_filter/Felis_catus.100x100kbp.masked.genome.bed \
      "$i" \
      > /home/llorenzo/depth_filter/"$i".100x100kbp.masked.depth
    done

With this, I will analyze coverage distribution (average, stdv, maxDepth, minDepth) and summarize information in a table. (See depth_script.R)

    scp llorenzo@genomics-a.ebd.csic.es:/home/llorenzo/depth_filter/*masked.depth .

Now that I have access to CESGA, I will repeat the previous step for the rest of lynx species. The first thing I have to do is copy Felis_catus.100x100kbp.masked.genome.bed to CESGA (from my laptop)

    scp Felis_catus.100x100kbp.masked.genome.bed csbiellf@ft2.cesga.es:/home/csic/bie/llf

Then I will have to run the previous code in an .sh format (See samtools_depth.sh). Upload the sh in the CESGA server:

    scp /Users/lorenalorenzo/github/selection_scan_lynx/Filtering/depth/executables/samtools_depth.sh csbiellf@ft2.cesga.es:/home/csic/bie/llf

and run it:

    sbatch samtool_depth.sh

As doing the loop in every sample will be high time consuming (more than 3 days), I will loop the sbatch for every sample, in other words, send 81 individual works (as we have 81 samples) to the server. Copy the new sh from the laptop to the server:

    scp /Users/lorenalorenzo/github/selection_scan_lynx/Filtering/depth/executables/samtools_depth_per_sample.sh csbiellf@ft2.cesga.es:/home/csic/bie/llf

Nevertheless, as lp samples are already calculated and lc as well (except the one mentioned above), I do not calculate on this samples to make the process faster and do not duplicate results. So, from 81 we pass to 81-11lp-18lc=52 samples(lr,ll and one lc). The lc sample will be send alone as it is difficult for me to include in the list of \$sample. In adittion, "c_ll_ca_0249" and "c_ll_ca_0253" are bad samples which are not counted in the total samples (so, 81 samples without the bad ones, remaining 52 samples).

So, for the rest:

    sample=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
    | grep -v "lp" | grep -v "lc" | rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

    for i in ${sample[@]}
    do
    echo "Calculating depth for $i"
    sbatch samtools_depth_per_sample.sh $i
    done

And for c_lc_zz_0003:

    lcsample=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/c_lc_zz_0003*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
    | rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

    for i in ${lcsample}
    do
    echo "Calculating depth for $i"
    sbatch samtools_depth_per_sample.sh $lcsample
    done

Download results in my laptop for R analyses:

    scp -r csbiellf@ft2.cesga.es:/mnt/lustre/scratch/home/csic/bie/llf .

At this point we may have evaluated per sample and species depth distribution (See depth_sample_distribution.R and depth_species_distribution.R). We have several ways to use this info:

A. Use min. and max. depth generated from distribution to filter each sample BED.

B. Eliminate samples with low quality of incongruences and then filter min and max depth separately We are not going to follow the A option because that way we will discard much "good info".

### 2. (B option) Eliminate samples with low quality of incongruences and then filter min and max depth separately

The first thing we are going to do is discard one sample: c_lc_yu_0007

#### 2.1 MAXIMUM depth filter

We decided to analyze depth distribution per species in order to avoid duplications but keep good information that may be isn't present in every sample due to quality issues. So that filter is a trade off.

First, we need to join each sample_depth info per species (remember the sample called c_lr_xx_0011 is a lc) and then analyze distribution in R similarly as done before with each sample.

    #in CESGA
    cd $LUSTRE/samples_subset_depth
    cut -f1,2 c_lc_ak_0015.100x100kbp.masked.depth > depth_rows

    species=(lc ll lp lr)

    for i in ${species[@]}
     do
      echo ${i}
      samples=($(ls $LUSTRE/samples_subset_depth/*masked.depth |rev | cut -d '/' -f 1 | grep "${i}" | rev | cut -d '.' -f 1)))
      for j in ${samples[@]}
       do
        echo ${j}
        cut -f3 ${j}.100x100kbp.masked.depth > ${j}.depthcolumn
      done
      paste depth_rows *${i}*.depthcolumn > ${i}_allsamples.depth
    done

Now I have to copy the outputs to my laptop in order to analyze in R:

    scp csbiellf@ft2.cesga.es:/mnt/lustre/scratch/home/csic/bie/llf/samples_subset_depth/\*_allsamples.depth .

I had to write a backlash () to use the pattern command (\*).

Once the per species dataset is in the laptop, I analyze them in R (See depth_species_distribution.R).

Now, we need a BED file with information of coverage per species (as the sum of each sample). That's not direct, neither trivial. As argued with Enrico, the pipeline will be: 1. Get a BED_coverage per sample (run BED_coverage.sh) 2. Get BED intervals with bedops (install with conda) and use --partition. One BED interval per species 3. Join cov_column per sample and sum up 4. Callable with the sp_depth.csv (maxDepth) 5. Filter the vcf with the callable sites

#### 2.1.1 Get a BED_coverage per sample

With bedtools genomecov we are calculating coverage per window in each sample, in order to laterly filter with min and max depth.

NOTE: We have to eliminate from the bamlist the samples we are not going to use, because of quality, depth or anything else: c_ll_ca_0249 c_ll_ca_0253 c_lc_yu_0007

And also have into consideration that: c_lr_xx_0011 is a lc

For that purpose, I moved the samples to a directory called "bad_samples" and rename c_lr_xx_0011 to c_lc_xx_0011 in the BED file (not in the BAM)

    sample=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
    | rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

    for i in ${sample[@]}
    do
    echo "BED for $i"
    sbatch bedtools_genomecov.sh $i
    done

Once we have a per sample bed, we need a per species bed (See species_bed.sh).

NOTE: In case we wanted to know how much is the genome represented in our samples bed:

    awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' <nombre_del_bed>

    species=(lc ll lp lr)

    for i in ${species[@]}
    do
    echo "BED for $i"
    sbatch species_bed.sh $i
    done

Although this previous command seems simple, it is presumably so computing-consuming, so we need to make it easier. For that purpose I am going to divide this work by chromosome, so we are going to get 21 files (20 chrm + rest) per species.

    # Array of Chromosomes and species
    CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
    species=(lc ll lp lr)

    # Create a BED with coordinates for each chromosome and in each species
    for i in ${species[@]}
      do
      for j in ${CHR[@]:0:20}
      do
      echo "bed for $j in $i"
      sbatch chr_sp_bed2.sh $i $j
      done
      done

    # For the rest of the info in the BED file:
    for i in ${species[@]}
      do
      echo "bed for rest in $i"
      sbatch rest_sp_bed.sh $i
      done

Once we have a bed per species chromosome, we wanted to join them so that we have ONE BED PER SPECIES!

    species=(lc ll lp lr)

    #First we join each chromosome
    for i in ${species[@]}
      do
      echo "bed for $i"
      sbatch sp_bed.sh $i
      done

    #Then we join the rest! NOTE THAT AANG AND KZ ARE ORDERED INVERSED THAN IN REF_GENOME (Felis_catus_9)
    for i in ${species[@]}
      do
      echo "Adding rest to ${i}"
      cat rest.${i}.bed >> ${i}.bed
      done

With this, we will have every sample ordered by chrm position but there will be overlapping regions. To control for this we are going to use bedops, which is a is a core tool for finding relationships between two or more genomic datasets. The partition operator splits all overlapping input regions into a set of disjoint segments. One or more input files may be provided; this option will segment regions from all inputs.

Because we don't have bedops in neither genomics nor CESGA server, we will use miniconda2 to create an environment where we will install it in CESGA:

    module load miniconda2
    conda create -n bedops
    conda activate bedops
    conda install -c bioconda bedops

Once created the environment called "bedops" and installed the program, we will do the partition (see bedops.sh)

We found out that bedops were having problems with first lines as they were not sorted: E.g: A1 0 1414 instead of A1 0 242 A1 0 242 A1 0 351 A1 0 351 A1 0 1414 To solve this problem, we are going to re-sort sp.bed for first, second and third field (see sorting_bed.sh) and after that, repeat bedops.sh

    #Send sorting sbatch
    species=(lc ll lp lr)

    for i in ${species[@]}
    do
    echo "sorting bed for $i"
    sbatch sorting_bed.sh $i
    done

    #After that, send BEDOPS sbatch
    species=(lc ll lp lr)

    for i in ${species[@]}
    do
    echo "bedops for $i"
    sbatch bedops.sh $i
    done

With the partitioned.bed we have the intervals for every species but not the COVERAGE. Now we have to add samples coverage for each interval and then sum up the columns (total sp coverage per interval). With bedtools intersect we can use the option sorted to make it less computing consuming. For using that, our inputs must be sorted, sp.bed it is but sample.bed isn't so we have to do it now before the bedtools_intersect.

    #Call species
    species=(lc ll lp lr)

    #Sort species samples
    for i in ${species[@]}
      do
      echo "sorting samples for $i"
      sbatch sorting_samples.sh $i
      done

    #bedtools_intersect
    for i in ${species[@]}
    do
    echo "Bedtools intersect for ${i}"
    sbatch bedtools_intersect.sh $i
    done

With this we will obtain the information of depth per sample in a species: CHR START END SAMPLE1 SAMPLE2 ... A1 0 10 12 0 A1 10 25 4 7 A1 25 50 8 3 We should have as columns as samples per species + 3 (chr, start, stop), but just in case we are going to check it out (REMEMBER: 19 lc, 32 ll, 11 lp, 18 lr):

    head -1 <FILENAME> | wc -w

The next step is sum this depth per interval, i.e sum columns from 4 to end (this end depends on the number of samples per sp)

    FS is FIELD SEPARATOR whereas OFS is OUTPUT FIELD SEPARATOR.
    #lynx_canadiensis
    awk '{FS="\t"; OFS="\t"; print $1, $2, $3, $4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22}' $LUSTRE/lc_sorted_partitioned_samples_coverage.bed > $LUSTRE/lc_sorted_partitioned_total_coverage.bed

    #lynx_lynx
    awk '{FS="\t"; OFS="\t"; print $1, $2, $3, $4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$28+$29+$30+$31+$32+$33+$34+$35}' $LUSTRE/ll_sorted_partitioned_samples_coverage.bed > $LUSTRE/ll_sorted_partitioned_total_coverage.bed

    #lynx_pardinus
    awk '{FS="\t"; OFS="\t"; print $1, $2, $3, $4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14}' $LUSTRE/lp_sorted_partitioned_samples_coverage.bed > $LUSTRE/lp_sorted_partitioned_total_coverage.bed

    #lynx_rufus
    awk '{FS="\t"; OFS="\t"; print $1, $2, $3, $4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21}' $LUSTRE/lr_sorted_partitioned_samples_coverage.bed > $LUSTRE/lr_sorted_partitioned_total_coverage.bed

Now we can assign callable vs. non callable sites. If maximum depth per species (See depth_species_distribution.R) is reached (1.5xmode), then this sites are "non-callable" (See max_depth_callable.sh)

YES/ NON-CALLABLE SITES IN BED: Lc 2283425907/38222070 (1.6%) Ll 2346130076/35895243 (1.5%) Lp 2197394627/68759452 (3%) Lr 2076467577/ (2%)

Once explored the callable sites in the bed file, we can apply the max_depth filter by removing variants non-callable from the vcf. Before this we had to separate the vcf per species (See vcf_per_species.sh)

    sbatch max_depth_callable.sh

    #Call species
    species=(lc ll lp lr)

    #Get NON-callable variants from vcf
    for i in ${species[@]}
     do
      echo ${i}
      grep "NO" /home/llorenzo/${i}_sorted_partitioned_total_coverage_defined.bed |
      bedtools intersect -header \
          -a /home/llorenzo/vcf/${i}_goodsamples_cat_ref_filter5.vcf \
          -b - \
          > /home/llorenzo/${i}_goodsamples_cat_ref_noncallable_regions.vcf
     done

    # check amount of variants lost in the vcf after max_depth filter.
    for i in ${species[@]}
     do
      echo ${i}
      cat ${i}_goodsamples_cat_ref_noncallable_regions.vcf | wc -l
     done

Non callable variants in each species: lc: 331338 ll: 314320 lp: 529649 lr: 524749

Now we need to eliminate those variants from the goodsamples_vcf (See bedtools_subtract_noncallable.sh)

At this point we have 21,956,842 variants

#### 2.2 MINIMUM depth filter

With the minimum filter, we are going to discard genotypes with less than 3 reads from the goodsamples_vcf

    # Filter DP variants -> change to missing data
    vcffilter -g "DP > 2" /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.vcf| vcffixup - \
        > /home/llorenzo/vcf_filtering/depth_filter/goodsamples_cat_ref.filter6.1.vcf

*NOTE: This filter puts a "." to non-callable instead of "./."*

Now I will explore the missingness variants by applying this filter in the goodsamples_vcf.

    bcftools filter -i "F_MISSING = 1" -Ov /home/llorenzo/vcf_filtering/depth_filter/goodsamples_cat_ref.filter6.1.vcf \
        > /home/llorenzo/vcf_filtering/depth_filter/goodsamples_cat_ref.filter6.1.1.vcf

The 6.1.1 file contains the variants completely missing for the all samples. In total we lost 3,662 variants. So we are going to eliminate those variants:

    # Subtract missing variants from allsamples vcf
      bedtools subtract -header \
          -a /home/llorenzo/vcf_filtering/depth_filter/goodsamples_cat_ref.filter6.1.vcf \
          -b /home/llorenzo/vcf_filtering/depth_filter/goodsamples_cat_ref.filter6.1.1.vcf \
          > /home/llorenzo/vcf_filtering/goodsamples_cat_ref.filter6.2.vcf

At this point the depth filter is DONE!!!! He have 21,953,180 variants left
