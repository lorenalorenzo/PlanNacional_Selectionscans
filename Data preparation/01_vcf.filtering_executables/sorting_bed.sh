#!/bin/bash
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#We want to create unsolapped intervals in each per species bed, in other words, join intervals with the same coverage.
#We found out that bedops were having problems with first lines as they were not sorted:
#E.g: A1 0 1414 instead of  A1 0 242
#     A1 0 242              A1 0 351
#     A1 0 351              A1 0 1414
#To solve this problem, we are going to re-sort sp.bed for first, second and third field.

sp=($(echo $1))

echo "Starting sort to ${sp}"
sort -k1,1 -k2,2n -k3,3n $LUSTRE/${sp}.bed > $LUSTRE/${sp}_sorted.bed
echo "Ended sort to ${sp}"
