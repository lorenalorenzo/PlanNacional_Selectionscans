#!/bin/bash
#SBATCH -t 3-10:00:00
#SBATCH -c 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=lorena.lorenzo.fdez@gmail.com

#############################################################
START=$(date)
echo "Samtools depth SCRIPT for $1 starting : $START"
#############################################################

############################################
## Samtools Depth launcher (samtools/1.10)##
############################################

# With this script I want to calculate depth at each position of every sample_bam (81)
# using samtools depth. Depth at all positions will be calculated (-a) within the
# regions randomly selected before (-b) (see depth_filtering.md for more detail).
# This will be run in a loop for each sample_bam

module load samtools

#####################################
## Calculating Depth with SAMtools ##
#####################################

# Loop of Samtools depth calculations for each sample_bam
for i in $STORE2/lynx_genome/lynx_data/CatRef_bams/*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
  do
  sample=($(echo $i | rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))
  echo "Calculating depth for $sample"
  samtools depth -a -b /home/csic/bie/llf/Felis_catus.100x100kbp.masked.genome.bed \
  "$i" \
  > $LUSTRE/"$sample".100x100kbp.masked.depth
done

###########################################################
END=$(date)
echo "Samtools depth SCRIPT for $1 ended : $END"
###########################################################
