#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#We want to create unsolapped intervals in each per species bed, in other words, join intervals with the same coverage.
#Name the species:
sp=($(echo $1))

#First we need to load the environment where we have bedops installed and activate it.
module load cesga/2018
module load miniconda2
conda activate bedops


#Calculate intervals with partition option of bedops
echo "BED partition for $sp"
bedops --partition $LUSTRE/${sp}_sorted.bed > $LUSTRE/"${sp}"_sorted_partitioned.bed
echo "Done partition for $sp"
