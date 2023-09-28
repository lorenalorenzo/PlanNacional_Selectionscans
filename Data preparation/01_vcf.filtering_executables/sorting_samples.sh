#!/bin/bash
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#First we call the species
sp=($(echo $1))

#Call the samples of this species (grep sp)
samples=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bed/samples_bed/*coverage.bed \
| grep "${sp}" | rev | cut -d '/' -f1 | rev | cut -d '.' -f 1))

#Sort by sample
for j in ${samples[@]}
do
  echo "sorting $j"
  sort -k1,1 -k2,2n -k3,3n \
  $STORE2/lynx_genome/lynx_data/CatRef_bed/samples_bed/${j}.coverage.bed \
  > $STORE2/lynx_genome/lynx_data/CatRef_bed/samples_bed/${j}_sorted.coverage.bed
done
