
library(tidyverse)
library(purrr)
library(readr)

setwd("/Users/yujin/Documents/")
expr_dir = "GTEx_Analysis_V6_eQTLInputFiles_geneLevelNormalizedExpressionMatrices"

# Load Gene-Expression Data

GTEx_data = lapply(list.files(path = expr_dir, pattern = "*.expr.txt",full.names = T), 
                   function(filename) {
                     data = read_table2(filename)
                     ID = data$Id
                     data = t(as.matrix(data[,-1]))
                     colnames(data) = ID
                     return(data)
                   })
tissue_names = sapply(list.files(path = expr_dir, pattern = "*.expr.txt"), 
                      function(s){gsub("_Analysis.expr.txt", "", s)})
names(GTEx_data) = tissue_names

# Find genes that are measured across all tissues
intersect_names <- colnames(GTEx_data[[tissue_names[1]]])
for (name in tissue_names) {
  intersect_names <- intersect(intersect_names, colnames(GTEx_data[[name]]))}
for (name in tissue_names) {
  GTEx_data[[name]] <- GTEx_data[[name]][, intersect_names]}

# Select tissues that have more than 150 samples
#n_samples <- unlist(lapply(GTEx_data, nrow))
#GTEx_data <- GTEx_data[n_samples >= 150]
#tissue_names <- tissue_names[n_samples >= 150] 