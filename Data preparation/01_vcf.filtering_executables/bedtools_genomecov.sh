#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

# With this script I want a BED with info about genome coverage (CHR  START_POSITION END_POSITION COVERAGE.
# So the output will be the depth per window of the genome in ONE sample.

module load bedtools

sample=($(echo $1))

echo "BED for ${sample}"

bedtools genomecov \
-ibam $STORE2/lynx_genome/lynx_data/CatRef_bams/${sample}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam -bga \
> $LUSTRE/"${sample}".coverage.bed

# The aim of this is to filter each sample with the min y max depth limits established before (See depth_filtering.md)
# However, we have finally decided not to do it in that way
