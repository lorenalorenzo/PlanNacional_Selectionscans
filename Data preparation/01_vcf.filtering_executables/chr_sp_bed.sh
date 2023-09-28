#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#Now that we have a per sample bed file, we want a per species bed file. For doing that,
#we need to join every sample bed of one species and order the content (field 1 and 2)

sp=($(echo $1))
chr=($(echo $2))

echo "Starting bed for "${chr}" in ${sp}"
cat /mnt/lustre/scratch/home/csic/bie/llf/samples_bed/*${sp}*.coverage.bed | grep "${chr}" | sort -k1,1 -k2,2n | uniq \
> /mnt/lustre/scratch/home/csic/bie/llf/"${chr}"."${sp}".bed
echo "Done bed for "${chr}" in "${sp}""
