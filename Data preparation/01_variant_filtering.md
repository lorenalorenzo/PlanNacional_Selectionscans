---
Title: Variant filtering steps for selection analysis in Lynx
Author: Lorena Lorenzo Fernández
Date: 11 March, 2021
---

Once we have our genotypes called, we want to filter the SNPs in order to get better information. The first 5 steps of this filters are extensively explained in Enrico's Github: <https://github.com/Enricobazzi/PlanNacional_Demographic_models_Oscar/blob/master/2.Variant_Filtering.md>

## (1,2,3,4,5) Applying General Filters

The script applying filters 1 through 5, which are independent of species and sequencing technology, and can therefore be applied to the whole dataset can be found at 2.Variant_Filtering-executables/General_Filter_1-2-3-4-5.sh

At this point we have 22,940,737 variants.

## (6) Filter by depth

This step is extensively explained in [depth_filtering.md](https://github.com/lorenalorenzo/selection_scan_lynx/blob/main/Filtering/depth/depth_filtering.md)

Basically we are going to filter maximum and minimum depth in different ways. For maximum we eliminate those variants where the maximum is reached (calculated by species) and for minimum we eliminate those genotypes where the minimum is not reached.

After applying depth filter we have 21,953,180 variants.

## (7) Under-represented, excessively missing variants

It's important to filter out variants which are missing completely in one or more species. For this I will divide the VCF into 4 different per-species VCFs and generate 4 lists of variants with no genotype called in any individual (of that species).

With the per-species VCF, I can use BCFtools again to filter out variants with missing genotypes for too many samples. Each species will have a different value of required number of samples with missing genotype (as the sample size is different for each species). By telling BCFtools to include only those variants (-i), the output can be used as a list of positions to exclude from the original VCF.

With this script I want to apply a filter based on data missingness to my dataset of 80 individuals, composed of 11 Lynx pardinus, 32 Lynx lynx, 18 Lynx rufus and 19 Lynx canadiensis. I will remove variants absent in at least 70% of the individuals (8lp, 22ll, 13lr, 13lc). (See missingness.sh)

lc_missing variants=115,057 ll_missing variants=87,115\
lp_missing variants=1,247,407\
lr_missing variants=69,576

After applying missingness filter we have 20,592,706 variants.

## (8) Filtering by het\>80% samples

Finally, to avoid false positive heterozygous calls due to mapping errors, we excluded any site where more than 80% of the samples are heterozygous (as in Kovalaskas et al. 2020).

With the per-species VCF, we need the chr, position and genotype info per sample. After that, we pass a filter for \>80% of heterozygosity (15lc, 26ll, 9lp, 14lr). Finally, we eliminate those variants from each vcf and test the variants left:

lc_het variants=4,177\
ll_het variants=5,258\
lp_het variants=6,197\
lr_het variants=5,770

From the total vcf we now have 20,578,917 variants
