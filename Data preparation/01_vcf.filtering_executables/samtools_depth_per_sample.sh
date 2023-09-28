#!/bin/bash
#SBATCH -t 1:30:00
#SBATCH -c 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=lorena.lorenzo.fdez@gmail.com

# With this script I want to calculate depth at each position of every sample_bam (81)
# using samtools depth. Depth at all positions will be calculated (-a) within the
# regions randomly selected before (-b) (see depth_filtering.md for more detail).

module load samtools

sample=($(echo $1))

samtools depth -a -b /home/csic/bie/llf/Felis_catus.100x100kbp.masked.genome.bed \
$STORE2/lynx_genome/lynx_data/CatRef_bams/${sample}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
> $LUSTRE/"${sample}".100x100kbp.masked.depth
