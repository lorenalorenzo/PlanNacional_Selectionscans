#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH --mem=80GB


module load cesga/2018
module load miniconda2
conda activate mummer_env

echo "Activated mummer"
nucmer -p $LUSTRE/Felis_catus_aligned_Lynx_canadensis $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa $STORE2/reference_genomes/Ref_Genome_LyCa/lc4.fa
echo "Done alignment cat-lynx"
