#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -c 1

# define species
species=(lc ll lp lr)

# define the table
table=/home/csic/bie/llf/sp_depth.csv

# assign YES vs NO (callable vs non-callable). See max_depth_callable.sh
for i in ${species[@]}
 do
  max=($(grep "${i}" ${table} | cut -d',' -f6))
  echo "${i} max depth is ${max}"
  cat $LUSTRE/${i}_sorted_partitioned_total_coverage.bed |
  awk -v max="${max}" '{FS="\t"; OFS="\t"; if ($4 > max) print $1, $2, $3, $4, "NO"; else print $1, $2, $3, $4, "YES";}' \
  > $LUSTRE/${i}_sorted_partitioned_total_coverage_defined.bed
done
