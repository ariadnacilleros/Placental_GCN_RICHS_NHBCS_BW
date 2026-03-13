### Perform association analysis 
### 10/03/25
### Ariadna

setwd("./")
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(stringr)
library(mediation)

# Read counts data used for module computation
counts <- readRDS("./counts_genenetwork.RDS")
dim(counts) # 673 samples and 15,082 genes
counts[1:5,1:5]

# Read metadata
pds <- read.csv("./final_sample_variables.csv")
head(pds)

# Make sure we have same sample order
identical(colnames(counts), pds$sampleid)
# TRUE

# Read modules
netwk <- readRDS("./netwk.RDS")

# Test association between eigengenes and BW percentile
eigen <- netwk$MEs
identical(rownames(eigen), pds$sampleid)
# TRUE
# Design model
design <- model.matrix(~ BW_percentile2 + Maternal_smoking + Baby_race + Gestational_age + Delivery_method + Infant_gender + Maternal_age + V1 + V2 + V3 + V4 + V5, data = pds)
fit <- lmFit(t(eigen), design)
fit <- eBayes(fit)
fit.results <- topTable(fit, coef = 2, number = Inf)
table(fit.results$adj.P.Val < 0.05)
rownames(fit.results[fit.results$adj.P.Val < 0.05, ])
# 7 associations:"ME23" "ME22" "ME8"  "ME6"  "ME19" "ME1"  "ME15"

# Get SE (t = logFC/se -> se = logFC/t)
fit.results$se <- fit.results$logFC/fit.results$t

# Save R object
saveRDS(object = fit.results, file = "modules_limma_results.RDS")

# Calculate CI
fit_coefs <- fit.results
fit_coefs$gene <- rownames(fit.results)
fit_coefs$lower_ci <- fit.results$logFC - qt(0.975, df=fit$df.residual) * fit.results$se
fit_coefs$upper_ci <- fit.results$logFC + qt(0.975, df=fit$df.residual) * fit.results$se

# Save R object
saveRDS(object = fit_coefs, file = "modules_limma_results_coefs.RDS")

# Perform differential expression analysis within each module
modules_numeric_lables <- netwk$colors
# Filter by modules of interest
modules_of_interest <- modules_numeric_lables[modules_numeric_lables == 23 | modules_numeric_lables == 22 | modules_numeric_lables == 8 | modules_numeric_lables == 6 | modules_numeric_lables == 19 | modules_numeric_lables == 1 | modules_numeric_lables == 15]
# Get counts per module
counts_list <- list()
for (i in unique(modules_of_interest)){
  genes <- names(modules_of_interest[modules_of_interest==i])
  counts_genes <- counts[genes, ]
  name_matrix <- paste("ME", i, sep = "")
  counts_list[[name_matrix]] <- counts_genes
}
# Make sure sample order 
for (i in c(1:length(counts_list))){
  print(i)
  print(identical(colnames(counts_list[[i]]), pds$sampleid))
}
# ALL TRUE
# Run differential expression analysis
# Store results
limma_results_list <- list()
for (i in seq_along(counts_list)) {
  expr <- counts_list[[i]]
  
  # Create design matrix
  design <- model.matrix(~ BW_percentile2 + Maternal_smoking + Baby_race + Gestational_age + Delivery_method + Infant_gender + Maternal_age + V1 + V2 + V3 + V4 + V5, data = pds)
  
  # Fit model
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  
  # Extract results for the response variable (always coef=2 if it's first after intercept)
  top <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  
  limma_results_list[[names(counts_list)[i]]] <- top
}
# Check consistency in # of genes tested per module
cat("All genes tested:", identical(dim(limma_results_list$ME15)[1], dim(counts_list$ME15)[1]), "\n")
cat("All genes tested:", identical(dim(limma_results_list$ME1)[1], dim(counts_list$ME1)[1]), "\n")
cat("All genes tested:", identical(dim(limma_results_list$ME19)[1], dim(counts_list$ME19)[1]), "\n")
cat("All genes tested:", identical(dim(limma_results_list$ME8)[1], dim(counts_list$ME8)[1]), "\n")
cat("All genes tested:", identical(dim(limma_results_list$ME22)[1], dim(counts_list$ME22)[1]), "\n")
cat("All genes tested:", identical(dim(limma_results_list$ME6)[1], dim(counts_list$ME6)[1]), "\n")
cat("All genes tested:", identical(dim(limma_results_list$ME23)[1], dim(counts_list$ME23)[1]), "\n")

# Save R object 
saveRDS(object = limma_results_list, file = "all_limma_results.RDS")

# Check significant results
significant_res <- data.frame()
for (module in names(limma_results_list)){
  print(module)
  res.rbind <- limma_results_list[[module]][limma_results_list[[module]]$adj.P.Val < 0.05, ]
  print(res.rbind)
  if (dim(res.rbind)[1] > 0){
    significant_res <- rbind(significant_res, data.frame(res.rbind, ME = module))
  }
}
cat("Differentially expressed genes:", dim(significant_res)[1], "\n") # 441 DEG 
cat("Downregulated genes:", dim(significant_res[significant_res$logFC < 0, ])[1], "\n") #128 downregulated
cat("Upregulated genes:", dim(significant_res[significant_res$logFC > 0, ])[1], "\n") # 313 upregulated

# Save R object
saveRDS(object = significant_res, file = "sign_limma_results.RDS")
