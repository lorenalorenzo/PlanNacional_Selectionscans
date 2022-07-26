#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#Now that we have a per chr_sp file, we want a per species bed file. For doing that,
#we need to join every chr_sp file.

sp=($(echo $1))
CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))

rm $LUSTRE/"${sp}".bed

for j in ${CHR[@]:0:20}
do
  echo "adding "${j}" to "${sp}""
  cat $LUSTRE/${j}.${sp}.bed >> $LUSTRE/${sp}.bed
done

echo "BED for "${sp}" done"
