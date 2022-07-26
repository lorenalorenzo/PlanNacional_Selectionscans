# **Title: MSMC2**

## Author: Lorena Lorenzo Fernández

### Date: 24 January 2022

I need to get **masked files** which basically is a bed file with the regions "callable". *In the Github pipeline they used a masked file per sample, is it important??*

-   Masked regions because of read depth:

I have the \$LUSTRE/sp_sorted_partitioned_total_coverage_defined.bed which includes the information about each site codified as YES (for callable) and NO (non-callable) but this is only for maximum depth filter. For minimum depth filter I eliminated genotypes with less than 3 reads (individual bed coverage), so we can't talk about sites non-callable because of low depth (¿? discuss with Enrico).

-   Masked repetitive regions:

This file is available in \$STORE2/reference_genomes/Felis_catus_Ref/Masked_Regions.bed
