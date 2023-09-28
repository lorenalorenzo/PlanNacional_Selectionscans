#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#Note that this script pretends to do the same as chr_sp_bed.sh but,
#because of computing issues, this code tries to make it simpler.

#Now that we have a per sample bed file, we want a per species bed file. For doing that,
#we need to join every sample bed of one species and order the content (field 1 and 2)

#First we call our variables mentioned in the sbatch loop
sp=($(echo $1))
chr=($(echo $2))

#This line removes the file notsorted.bed just in case we run the code several times.
#Note that normally when you run a code twice or more, the output will be automatically substituted.
#In this case, as we say ">>" it means we want to add the info in each for loop, not rewrite the file.
#So, if we dont remove each time we start this sh, it will add information to the same file instead of start it again.
rm /mnt/lustre/scratch/home/csic/bie/llf/"${chr}"."${sp}".notsorted.bed

#new variable with the files of the species being analysed in each sbatch loop.
samples_sp=($(ls /mnt/lustre/scratch/home/csic/bie/llf/samples_bed/*${sp}*.coverage.bed))

for i in ${samples_sp[@]}
do
  echo "extracting "${chr}" from "$i""
  grep "${chr}" $i >> /mnt/lustre/scratch/home/csic/bie/llf/"${chr}"."${sp}".notsorted.bed
done

echo "Done bed for "${chr}" in "${sp}" but not sorted"

cat /mnt/lustre/scratch/home/csic/bie/llf/"${chr}"."${sp}".notsorted.bed | sort -k1,1 -k2,2n | uniq \
>  /mnt/lustre/scratch/home/csic/bie/llf/"${chr}"."${sp}".bed

echo "Done bed for "${chr}" in "${sp}" sorted"
