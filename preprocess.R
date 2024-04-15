
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

# There are in total 16622 genes. 
# We randomly sample 5000 to reduce the size of the data set so that we can include it into our package.  
new_intersect_names = intersect_names[sample.int(length(intersect_names), size = 5000, replace = FALSE)]
for (name in tissue_names){
  GTEx_data[[name]] <- GTEx_data[[name]][, new_intersect_names]
}

save(GTEx_data, file = 'Documents/dlm/data/GTEx_data.rda')
