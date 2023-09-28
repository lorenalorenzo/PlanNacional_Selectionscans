#!/bin/bash
#SBATCH -t 1-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=lorena.lorenzo.fdez@gmail.com

##############################################
## BWA mem and BAM post-processing launcher ##
##############################################

# This program runs BWA mem and then standard filtering steps for per-sample read alignment
# to the Felis catus reference genome for serval and lion samples. 
# It will be run on the finis terrae II server of CESGA.


# The following softwares are used:
# BWA
# SAMtools
# Picard-tools
# GATK

module load cesga/2018

module load bwa/0.7.17
module load samtools/1.10
module load picard/2.21.8
module load gatk/3.7-0-gcfedb67

###################################
## VARIABLE and PATHS definition ##
###################################
# This is an adaptation of the lynx mapping code to outgroup species (serval and lion), 
# which data have been download from NCBI and we want to align to Cat reference genome in its 9 version 
# to make it comparable to my lynx vcf for polarization.
# SRR6071636 --> serval code
# SRR836361 --> lion code

#Defining paths
THREADS=24
REF=/mnt/netapp1/Store_csebdjgl/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa
INfastq_PATH=/mnt/netapp1/Store_csebdjgl/reference_genomes/lion_serval_polarization/SRR836361
INfastqR1=($(ls $INfastq_PATH | grep "_1"))
INfastqR2=($(ls $INfastq_PATH | grep "_2"))
OUT=/mnt/netapp1/Store_csebdjgl/reference_genomes/lion_serval_polarization/SRR836361
SAMPLE=(lion)

#################################
## Alignment and BAM filtering ##
#################################

# Mapping 
echo " - Mapping lion -"
bwa mem $REF $INfastq_PATH/${INfastqR1} $INfastq_PATH/${INfastqR2} \
-t $THREADS | samtools view -hbS -@ $THREADS - -o $OUT/${SAMPLE}.cat_ref.bam

# Sorting
echo " - Sorting ${SAMPLE} -"
samtools sort -@ $THREADS $OUT/${SAMPLE}.cat_ref.bam -o $OUT/${SAMPLE}.cat_ref.sorted.bam \
&& rm $OUT/${SAMPLE}.cat_ref.bam

# Adding READ Groups
echo " - Adding READ Groups of ${SAMPLE} -"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I=$OUT/${SAMPLE}.cat_ref.sorted.bam \
O=$OUT/${SAMPLE}_cat_ref_sorted_rg.bam \
RGID=SRR836361 RGLB=Genome_of_Insert_Size_400bp \
RGPL=Illumina RGPU=SRS417526 RGSM=${SAMPLE} \
VALIDATION_STRINGENCY=SILENT && rm $OUT/${SAMPLE}.cat_ref.sorted.bam

# Marking Duplicates, Re-Sorting and Indexing
echo " - Marking Duplicates of ${SAMPLE} -"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
METRICS_FILE=$OUT/${SAMPLE}_rmdup.txt \
I=$OUT/${SAMPLE}_cat_ref_sorted_rg.bam \
O=$OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup.bam \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800
rm $OUT/${SAMPLE}_cat_ref_sorted_rg.bam

echo " - Re-Sorting ${SAMPLE} -"
samtools sort $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup.bam \
-@ $THREADS -o $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted.bam
rm $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup.bam

echo " - Indexing ${SAMPLE} -"
samtools index $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted.bam

# Realigning:
# RealignerTargetCreator
echo " - Realigner Target Creator on ${SAMPLE} -"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-nt $THREADS -R $REF -I $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted.bam \
-o $OUT/${SAMPLE}_realignertargetcreator.intervals
# IndelRealigner
echo " - Realigning INDELS of ${SAMPLE} -"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner \
-R $REF -targetIntervals $OUT/${SAMPLE}_realignertargetcreator.intervals \
-I $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted.bam \
-o $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
rm $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted.bam

# Get stats
echo " - Getting stats of ${SAMPLE} -"
samtools flagstat $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
> $OUT/${SAMPLE}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.stats