#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#Now that we have the unsolapped intervals we want the per sample depth corresponding to that intervals.
#With this script we want to run a per species intersect with the samples_bed and the partitioned_bed,
#conserving the intervals of the partitioned_bed and the depth data of each sample.

sp=($(echo $1))

#Copy the intervals file with the name we want to our per species sample coverage bed file.
cp $LUSTRE/${sp}_sorted_partitioned.bed $LUSTRE/${sp}_sorted_partitioned_samples_coverage.bed

#Create a variable with the sample names
samples=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bed/samples_bed/*coverage.bed \
| grep "${sp}" | rev | cut -d '/' -f1 | rev | cut -d '.' -f 1 | cut -d '_' -f 1-4 | uniq))

#Load bedtools
module load cesga/2018
module load bedtools

#With bedtools intersect add the coverage info to the partitioned_bed.
#We used -sorted option because it is recommended for really large bed previously sorted
for j in ${samples[@]}
do
echo "Bedtools intersect coverage for ${j}"
  bedtools intersect -wa -wb \
        -a $LUSTRE/${sp}_sorted_partitioned.bed \
        -b $STORE2/lynx_genome/lynx_data/CatRef_bed/samples_bed/${j}_sorted.coverage.bed \
        -sorted \
        | cut -f7 \
        > $LUSTRE/${sp}_tmp
  paste $LUSTRE/${sp}_sorted_partitioned_samples_coverage.bed $LUSTRE/${sp}_tmp > $LUSTRE/${sp}_tmp2
  mv $LUSTRE/${sp}_tmp2 $LUSTRE/${sp}_sorted_partitioned_samples_coverage.bed
  done
echo "Done bed coverage for ${sp}"

#Try code with one sample
#bedtools intersect -wa -wb \
#    -a $LUSTRE/lc_sorted_partitioned.bed \
#    -b $STORE2/lynx_genome/lynx_data/CatRef_bed/samples_bed/c_lc_ak_0015.coverage.bed \
#    -sorted \
#    | cut -f7
#Maybe I have to do it per chr because of computing consuming
#Try code in one chr of one sample
#grep "A1" $LUSTRE/lc_partitioned.bed | bedtools intersect -wa -wb -a - \
#  -b <(grep "A1" $STORE2/lynx_genome/lynx_data/CatRef_bed/samples_bed/c_lc_ak_0015.coverage.bed)
