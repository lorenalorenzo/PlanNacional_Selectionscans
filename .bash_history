pwd
file_list <- ls files/*xpehh*outliers
file_list <- ls *xpehh*outliers
file_list <- ls \.*xpehh*outliers
file_list= ls *xpehh*outliers
echo ${file_list[@]}
file_list=(ls *xpehh*outliers)
echo ${file_list[@]}
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' files/${i} | sed 's/chr//g' > files/${i}.bed    
    bedtools merge       -d 1000       -i files/${i}.bed       > ${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a files/${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > ${i}_annotated ;      echo "$i annotated";   done   
echo ${file_list[@]}
file_list=(ls *xpehh*outliers)
echo ${file_list[@]}
file_list=$(ls *xpehh*outliers)
echo ${file_list[@]}
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' files/${i} | sed 's/chr//g' > files/${i}.bed    
    bedtools merge       -d 1000       -i files/${i}.bed       > ${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a files/${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > ${i}_annotated ;      echo "$i annotated";   done   
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' ${i} | sed 's/chr//g' > ${i}.bed    
    bedtools merge       -d 1000       -i files/${i}.bed       > ${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a files/${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > ${i}_annotated ;      echo "$i annotated";   done   
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' ${i} | sed 's/chr//g' > ${i}.bed    
    bedtools merge       -d 1000       -i ${i}.bed       > ${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a ${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > ${i}_annotated ;      echo "$i annotated";   done   
file_list=(ls files/*xpehh*outliers)
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' ${i} | sed 's/chr//g' > files/${i}.bed    
    bedtools merge       -d 1000       -i files/${i}.bed       > files/${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a files/${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated ;      echo "$i annotated";   done   
file_list=$(ls files/*xpehh*outliers)
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' ${i} | sed 's/chr//g' > files/${i}.bed    
    bedtools merge       -d 1000       -i files/${i}.bed       > files/${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a files/${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated ;      echo "$i annotated";   done   
pwd
ls
echo ${file_list[@]}
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' ${i} | sed 's/chr//g' > ${i}.bed    
    bedtools merge       -d 1000       -i ${i}.bed       > ${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a ${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > ${i}_annotated ;      echo "$i annotated";   done   
for i in ${file_list[@]};   do     echo "$i";     awk 'NR>1 {print $1"\t"$2-1"\t"$2}' ${i} | sed 's/chr//g' > ${i}.bed    
    bedtools merge       -d 100000       -i ${i}.bed       > ${i}_genomic_regions;       echo "$i merged"
    bedtools intersect      -a ${i}_genomic_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > ${i}_annotated ;      echo "$i annotated";   done   
      
bedtools intersect        -a files/lc_lassi_ihs_regions        -b files/ll_lassi_ihs_regions files/lp_lassi_ihs_regions files/lr_lassi_ihs_regions        -wo        > files/lc_repeated_regions
cut -f 3 c_lassi_ihs_regions
cut -f 3 files/lc_lassi_ihs_regions
cut -f 1-3 files/lc_lassi_ihs_regions
bedtools intersect        -a <(cut -f 1-3 files/lc_lassi_ihs_regions)        -b <(cut -f 1-3 files/ll_lassi_ihs_regions) <(cut -f 1-3 files/lp_lassi_ihs_regions) <(cut -f 1-3 files/lr_lassi_ihs_regions)        -names ll lp lr        -sorted        > files/lc_repeated_regions
bedtools intersect -wa -wb        -a <(cut -f 1-3 files/lc_lassi_ihs_regions)        -b <(cut -f 1-3 files/ll_lassi_ihs_regions) <(cut -f 1-3 files/lp_lassi_ihs_regions) <(cut -f 1-3 files/lr_lassi_ihs_regions)        -names ll lp lr        -sorted        > files/lc_repeated_regions
bedtools intersect -wa -wb        -a <(cut -f 1-3 files/ll_lassi_ihs_regions)        -b <(cut -f 1-3 files/lc_lassi_ihs_regions) <(cut -f 1-3 files/lp_lassi_ihs_regions) <(cut -f 1-3 files/lr_lassi_ihs_regions)        -names lc lp lr        -sorted        > files/ll_repeated_regions  
bedtools intersect -wa -wb        -a <(cut -f 1-3 files/lp_lassi_ihs_regions)        -b <(cut -f 1-3 files/lc_lassi_ihs_regions) <(cut -f 1-3 files/ll_lassi_ihs_regions) <(cut -f 1-3 files/lr_lassi_ihs_regions)        -names lc ll lr        -sorted        > files/lp_repeated_regions  
       
bedtools intersect -wa -wb        -a <(cut -f 1-3 files/lr_lassi_ihs_regions)        -b <(cut -f 1-3 files/lc_lassi_ihs_regions) <(cut -f 1-3 files/ll_lassi_ihs_regions) <(cut -f 1-3 files/lp_lassi_ihs_regions)        -names lc ll lp        -sorted        > files/lr_repeated_regions   
species=(lc ll lp lr)
for i in ${species[@]};   do     echo "$i";     bedtools intersect      -a files/${i}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_repeated_regions_annotated ;      echo "$i annotated";   done   
for i in ${species[@]};   do     echo "$i";     bedtools intersect      -a files/${i}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated_tmp ;      echo "$i annotated"
     awk '$10=="gene"' ${i}_annotated_tmp  > ${i}_gene_tmp;      echo "$i gene filtered"   
    cut -f 1:8,11:12,16 ${i}_gene_tmp > ${i}_repeated_genes;   done   
for i in ${species[@]};   do     echo "$i";     bedtools intersect      -a files/${i}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated_tmp ;      echo "$i annotated"
     awk '$10=="gene"' ${i}_annotated_tmp  > ${i}_gene_tmp;      echo "$i gene filtered"   
    cut -f1:8,11:12,16 ${i}_gene_tmp > ${i}_repeated_genes;   done   
for i in ${species[@]};   do     echo "$i";     bedtools intersect      -a files/${i}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated_tmp ;      echo "$i annotated"
     awk '$10=="gene"' ${i}_annotated_tmp  > ${i}_gene_tmp;      echo "$i gene filtered"   
    cut -f 1-8,11-12,16 ${i}_gene_tmp > ${i}_repeated_genes;   done   
for i in ${species[@]};   do     echo "$i";     bedtools intersect      -a files/${i}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated_tmp ;      echo "$i annotated"
     awk '$10=="gene"' files/${i}_annotated_tmp  > files/${i}_gene_tmp;      echo "$i gene filtered"   
    cut -f 1-8,11-12,16 files/${i}_gene_tmp > files/${i}_repeated_genes;   done   
grep "Name" files/.DS_Store files/Felis_catus.Felis_catus_9.0.97.gff3 files/PCA files/VCFs files/chr_size.txt files/lc_annotated_outliers files/lc_annotated_tmp files/lc_gene_outliers files/lc_gene_tmp files/lc_genes_selection.csv files/lc_genomic_regions_annotated files/lc_ihs_4outliers files/lc_ihs_4outliers.bed files/lc_ihs_scan files/lc_ihs_top10 files/lc_lassi_5%outliers files/lc_lassi_5%outliers.bed files/lc_lassi_ihs_regions files/lc_lassi_top10 files/lc_lassi_windows_genomic_regions files/lc_ll_xpehh_scan files/lc_ll_xpehh_scan_outliers files/lc_ll_xpehh_scan_outliers.bed files/lc_ll_xpehh_scan_outliers_annotated files/lc_ll_xpehh_scan_outliers_genomic_regions files/lc_lp_xpehh_scan files/lc_lp_xpehh_scan_outliers files/lc_lp_xpehh_scan_outliers.bed files/lc_lp_xpehh_scan_outliers_annotated files/lc_lp_xpehh_scan_outliers_genomic_regions files/lc_lr_xpehh_scan files/lc_lr_xpehh_scan_outliers files/lc_lr_xpehh_scan_outliers.bed files/lc_lr_xpehh_scan_outliers_annotated files/lc_lr_xpehh_scan_outliers_genomic_regions files/lc_missing_gts.csv files/lc_repeated_genes files/lc_repeated_regions files/lc_repeated_regions_annotated files/lc_results_table_representation files/lc_salti.lassip.hap.out files/lc_top10 files/lcdata_ihs files/lcdata_lassi files/ll_annotated_outliers files/ll_annotated_tmp files/ll_gene_outliers files/ll_gene_tmp files/ll_genes_selection.csv files/ll_genomic_regions_annotated files/ll_ihs_4outliers files/ll_ihs_4outliers.bed files/ll_ihs_scan files/ll_ihs_top10 files/ll_lassi_5%outliers files/ll_lassi_5%outliers.bed files/ll_lassi_ihs_regions files/ll_lassi_top10 files/ll_lassi_windows_genomic_regions files/ll_lp_xpehh_scan files/ll_lp_xpehh_scan_outliers files/ll_lp_xpehh_scan_outliers.bed files/ll_lp_xpehh_scan_outliers_annotated files/ll_lp_xpehh_scan_outliers_genomic_regions files/ll_lr_xpehh_scan files/ll_lr_xpehh_scan_outliers files/ll_lr_xpehh_scan_outliers.bed files/ll_lr_xpehh_scan_outliers_annotated files/ll_lr_xpehh_scan_outliers_genomic_regions files/ll_missing_gts.csv files/ll_repeated_genes files/ll_repeated_regions files/ll_repeated_regions_annotated files/ll_results_table_representation files/ll_salti.lassip.hap.out files/ll_top10 files/lldata_ihs files/lldata_lassi files/lp_annotated_outliers files/lp_annotated_tmp files/lp_gene_outliers files/lp_gene_tmp files/lp_genes_selection.csv files/lp_genomic_regions_annotated files/lp_ihs_4outliers files/lp_ihs_4outliers.bed files/lp_ihs_scan files/lp_ihs_top10 files/lp_lassi_5%outliers files/lp_lassi_5%outliers.bed files/lp_lassi_ihs_regions files/lp_lassi_top10 files/lp_lassi_windows_genomic_regions files/lp_lr_xpehh_scan files/lp_lr_xpehh_scan_outliers files/lp_lr_xpehh_scan_outliers.bed files/lp_lr_xpehh_scan_outliers_annotated files/lp_lr_xpehh_scan_outliers_genomic_regions files/lp_missing_gts.csv files/lp_repeated_genes files/lp_repeated_regions files/lp_repeated_regions_annotated files/lp_results_table_representation files/lp_salti.lassip.hap.out files/lp_top10 files/lpdata_ihs files/lpdata_lassi files/lr_annotated_outliers files/lr_annotated_tmp files/lr_gene_outliers files/lr_gene_tmp files/lr_genes_selection.csv files/lr_genomic_regions_annotated files/lr_ihs_4outliers files/lr_ihs_4outliers.bed files/lr_ihs_scan files/lr_ihs_top10 files/lr_lassi_5%outliers files/lr_lassi_5%outliers.bed files/lr_lassi_ihs_regions files/lr_lassi_top10 files/lr_lassi_windows_genomic_regions files/lr_missing_gts.csv files/lr_repeated_genes files/lr_repeated_regions files/lr_repeated_regions_annotated files/lr_results_table_representation files/lr_salti.lassip.hap.out files/lr_top10 files/lrdata_ihs files/lrdata_lassi files/ls.bed files/ls_annotated files/ls_genomic_regions files/missing_ind_rate.imiss files/samples_table.xlsx files/top10_selectedregions_lynxes.numbers files/top10_selectedregions_lynxes.xls eeee
species=(lc ll lp lr)
for sp in ${species[@]};   do
       bedtools intersect          -a <(grep -v "chr" files/${sp}_lassi_ihs_regions)          -b files/Felis_catus.Felis_catus_9.0.97.gff3          -wa -wb          > files/${sp}_annotated_tmp ;         echo "$sp annotated"        
      awk '$15=="gene"' files/${sp}_annotated_tmp > ${sp}_gene_tmp;       echo "$sp gene filtered"     done
cut -f-12,16,17,21- ll_gene_tmp
cut -f-12,16,17,21- ll_gene_tmp | less -S
cut -f-12,16,17,21- ${sp}_gene_tmp | tr ' ' '\t' > ${sp}_filtered_tmp
species=(lc ll lp lr)
for sp in ${species[@]};   do
       bedtools intersect          -a <(grep -v "chr" files/${sp}_lassi_ihs_regions)          -b files/Felis_catus.Felis_catus_9.0.97.gff3          -wa -wb          > files/${sp}_annotated_tmp ;         echo "$sp annotated"        
      awk '$15=="gene"' files/${sp}_annotated_tmp > files/${sp}_gene_tmp;       echo "$sp gene filtered"    
    cut -f-12,16,17,21- files/${sp}_gene_tmp | tr ' ' '\t' > ${sp}_filtered_tmp    
    cut -d';' -f1,2 ${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($16 ~ /Name=[[:alnum:]]/) $16=$16; else $16="NA"; print $0}' > ${sp}_names_tmp     
    paste <(cut -d" " -f-14 ${sp}_names_tmp) <(cut -d" " -f15 ${sp}_names_tmp | cut -d":" -f2) <(cut -d" " -f16 ${sp}_names_tmp | cut -d"=" -f2) |  tr "\t" " "  >> ${sp}_columns_tmp    
    echo -e "chr start end windows snps snps_density lassi_mean lassi_max lassi_sum ihs_mean ihs_max ihs_sum gene_start gene_end ensembl_id gene_name" | cat - ${sp}_columns_tmp > ${sp}_genomic_regions_annotated      
       rm *tmp;     done
cut -f 16 lc_genomic_regions_annotated
cut -f -16 lc_genomic_regions_annotated
cut -f 15 lc_genomic_regions_annotated
cut -d " " -f 16 lc_genomic_regions_annotated
cut -d " " -f 16 lc_genomic_regions_annotated | uniq
cut -d " " -f 16 lc_genomic_regions_annotated | sort | uniq
cut -d " " -f 16 lc_genomic_regions_annotated | sort | uniq | wc -l
cut -d " " -f 16 ll_genomic_regions_annotated | sort | uniq | wc -l
cut -d " " -f 16 lp_genomic_regions_annotated | sort | uniq | wc -l
cut -d " " -f 16 lr_genomic_regions_annotated | sort | uniq | wc -l
remove files/*tmp
rm files/*tmp
for i in ${species[@]};   do     echo "$i";     bedtools intersect      -a files/${i}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated_tmp ;      echo "$i annotated"
     awk '$10=="gene"' files/${i}_annotated_tmp  > files/${i}_gene_tmp;      echo "$i gene filtered"   
    cut -f 1-8,11-12,16 files/${i}_gene_tmp > files/${i}_repeated_genes
    cut -d';' -f1,2 ${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($12 ~ /Name=[[:alnum:]]/) $12=$12; else $12="NA"; print $0}' > files/${sp}_rep_genes_names_tmp ; done
for i in ${species[@]};   do     echo "$i";     bedtools intersect      -a files/${i}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${i}_annotated_tmp ;      echo "$i annotated"
     awk '$10=="gene"' files/${i}_annotated_tmp  > files/${i}_gene_tmp;      echo "$i gene filtered"   
    cut -f 1-8,11-12,16 files/${i}_gene_tmp > files/${i}_filtered_tmp
    cut -d';' -f1,2 files/${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($12 ~ /Name=[[:alnum:]]/) $12=$12; else $12="NA"; print $0}' > files/${sp}_rep_genes_names_tmp ; done
for sp in ${species[@]};   do     echo "$sp";     bedtools intersect      -a files/${sp}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${sp}_annotated_tmp ;      echo "$sp annotated"
     awk '$10=="gene"' files/${sp}_annotated_tmp  > files/${sp}_gene_tmp;      echo "$sp gene filtered"   
    cut -f 1-8,11-12,16 files/${sp}_gene_tmp > files/${sp}_filtered_tmp
    cut -d';' -f1,2 files/${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($12 ~ /Name=[[:alnum:]]/) $12=$12; else $12="NA"; print $0}' > files/${sp}_rep_genes_names_tmp ; done
for sp in ${species[@]};   do     echo "$sp";     bedtools intersect      -a files/${sp}_repeated_regions      -b files/Felis_catus.Felis_catus_9.0.97.gff3      -wa -wb       > files/${sp}_annotated_tmp ;      echo "$sp annotated"
     awk '$10=="gene"' files/${sp}_annotated_tmp  > files/${sp}_gene_tmp;      echo "$sp gene filtered"   
    cut -f 1-8,11-12,16 files/${sp}_gene_tmp > files/${sp}_filtered_tmp
    cut -d';' -f1,2 files/${sp}_filtered_tmp | tr ';' '\t'  | awk '{if ($12 ~ /Name=[[:alnum:]]/) $12=$12; else $12="NA"; print $0}' > files/${sp}_rep_genes_names_tmp     
    paste <(cut -d" " -f-10 files/${sp}_rep_genes_names_tmp) <(cut -d" " -f11 files/${sp}_rep_genes_names_tmp | cut -d":" -f2) <(cut -d" " -f12 files/${sp}_rep_genes_names_tmp | cut -d"=" -f2) |  tr "\t" " "  >> files/${sp}_rep_columns_tmp    
    echo -e "chr start end sp chr_2 start_2 end_2 gene_chr gene_start gene_end ensembl_id gene_name" | cat - files/${sp}_rep_columns_tmp > files/${sp}_repetitive_genomic_regions_annotated      
       rm *tmp;     done
rm files/
rm files/*tmp
cut -f 12 files/lc_repetitive_genomic_regions_annotated | sort | uniq
cut -d " " -f 12 files/lc_repetitive_genomic_regions_annotated | sort | uniq
cut -d " " -f 12 files/lc_repetitive_genomic_regions_annotated | sort | uniq | wc -l
cut -d " " -f 12 files/ll_repetitive_genomic_regions_annotated | sort | uniq | wc -l
cut -d " " -f 12 files/lp_repetitive_genomic_regions_annotated | sort | uniq | wc - l
cut -d " " -f 12 files/lp_repetitive_genomic_regions_annotated | sort | uniq | wc - l
cut -d " " -f 12 files/lp_repetitive_genomic_regions_annotated | sort 
cut -d " " -f 12 files/lr_repetitive_genomic_regions_annotated | sort | uniq | wc - l
cut -d " " -f 12 files/lr_repetitive_genomic_regions_annotated | sort | uniq | wc -l
