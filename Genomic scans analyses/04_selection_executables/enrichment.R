# Script: functional_enrichment.R
# Usage: Rscript functional_enrichment.R <input_file1> <input_file2>
# 
# Description:
# This script performs some operations on two input files provided as command-line arguments.
# 
# Arguments:
# - <input_file1>: Path to the first input file.
# - <input_file2>: Path to the second input file.
# 
# Example:
# To run the script, provide paths to two input files as command-line arguments:
# Rscript my_script.R path/to/input_file1.csv path/to/input_file2.csv
# 
# Author: Your Name
# Date: April 15, 2024

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there are at least two arguments
if (length(args) < 2) {
  stop("You need to provide two input file paths.")
}

# Access the input file paths
input_file1 <- args[1]
input_file2 <- args[2]

# Now you can work with input_file1 and input_file2
# For example, you could read them into R using read.csv(), read.table(), etc.
data1 <- read.csv(input_file1)
data2 <- read.csv(input_file2)

# Perform some operations with the data...


# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiomaRt")
install.packages("topGO")
library(topGO)

# Check if the correct number of command-line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  stop("Usage: Rscript enrichment_analysis_script.R <gene_list_file> <sp_annotation_id>")
}

# Get the gene list file path from the command-line arguments
gene_list_file <- commandArgs(trailingOnly = TRUE)[1]

# Load the gene list
gene_list <- readLines(gene_list_file)

# Provide the info of the species to get the correct annotation
sp_annotation_id <- commandArgs(trailingOnly = TRUE)[2]



# Rest of the script remains the same...
#load the functions
annotation <- function()
        #read from ensembl.org your annotation:
        ensembl <- useMart("ensembl", dataset = "sp_annotation_id")  
        #extract the GOterms for every ensembl_id
        ensembl_to_go <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id"), mart = ensembl)
        #make a list with every GO terms per ENSEMBL_ ID. 
        go_list <- split(ensembl_to_go$go_id, ensembl_to_go$ensembl_gene_id)
  
enrichment <- function(annotation_file, gene_list, GO_category)



#My purpose is to make a function that does an enrichment analysis using the topGo package.
# I need an annotattion file and a list of genes to do this. 
#The annotation file is a file that contains the GO terms for each gene (gene name or ensembl id) and this file can be directly obtained from the biomart website. 
#The list of genes is a list of genes that you want to do the enrichment analysis on. In this case we need a vector of gene names.
#Then, we have to cross the annotation file with the list of genes to get the GO terms for the genes in the list.
#Then, we have to do the enrichment analysis using the topGO package. 
#For that, we need to create the topGO object, that needs the GO terms for the genes in the list, the GO terms for all the genes in the genome and the GO category you want to test.
#Afterthat, we can run the overrepresentation analyses and get the results as a table.
enrichment <- function(annotation_file, gene_list, GO_category)
  {
  #Load the libraries
  library(topGO)
  library(biomaRt)
  library(dplyr)
  #Read the annotation file
  annotation <- read.table(annotation_file, header = TRUE, sep = "\t")
  #Read the gene list
  gene_list <- read.table(gene_list, header = TRUE, sep = "\t")
  #Cross the annotation file with the gene list
  gene_list_annotated <- inner_join(gene_list, annotation, by = c("Gene.name" = "Gene.name"))
  #Create the topGO object
  GOdata <- new("topGOdata",
                description = "GO enrichment analysis",
                ontology = "BP",
                allGenes = gene_list_annotated$Gene.name,
                geneSel = gene_list_annotated$Gene.name,
                annot = annFUN.org,
                mapping = "org.Hs.eg.db",
                ID = "Entrez",
                nodeSize = 10)
  #Run the overrepresentation analysis
  result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  #Get the results as a table
  result_table <- GenTable(GOdata,
                           classicFisher = result,
                           orderBy = "classicFisher",
                           ranksOf = "classicFisher",
                           topNodes = 10)
  #Return the table
  return(result_table)
}