### Perform gene set enrichment analysis with enrichR BioPlanet
### 23/09/25
### Ariadna

setwd("./")
library(hypeR)
library(msigdbr)
library(enrichR)
library(rio)
library(openxlsx)

# Load gene network results
res <- readRDS("netwk.RDS")

# Format wgcna output
convert_to_module_lists <- function(gene_module_vector) {
  
  # Get genes names
  clean_gene_names <- names(gene_module_vector)
  
  # Convert numeric module labels to ME-style labels
  module_labels <- paste0("ME", gene_module_vector)
  
  # Use split to group cleaned gene names by module
  module_lists <- split(clean_gene_names, module_labels)
  
  return(module_lists)
}

res_list <- convert_to_module_lists(res$colors)
head(res_list)

# Get background
background <- as.vector(unlist(res_list))

# Remove ME0
res_list <- res_list[-1]

# Run enrichR
# BioPlanet
enrichr_bioplanet<- list()
for (module in names(res_list)){
  enrichr_res <- enrichr(genes = as.vector(unlist(res_list[module])), databases = "BioPlanet_2019", background = background)
  enrichr_bioplanet[[module]] <- enrichr_res
}

file <- "./enrichR_bioplanet.xlsx"
wb <- createWorkbook()
for (module in names(enrichr_bioplanet)) {
  addWorksheet(wb, sheetName = module)
  writeData(wb, sheet = module, x = enrichr_bioplanet[[module]][[1]], colNames = TRUE)
}
saveWorkbook(wb, file = file, overwrite = TRUE)
