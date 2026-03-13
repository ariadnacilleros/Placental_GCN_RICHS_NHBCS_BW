### Perform mediation analysis 
### 23/09/25
### Ariadna

setwd("./")
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(stringr)
library(mediation)
set.seed(1234)

# Load R objects needed
significant_res <- readRDS("./sign_limma_results.RDS")
netwk <- readRDS("./netwk.RDS")
pds <- read.csv("./final_sample_variables.csv")
counts <- readRDS("./counts_genenetwork.RDS")

# Perform mediation analysis
mediation_results_f <- data.frame()
mediation_results_r <- data.frame()

for (sign.module in unique(significant_res$ME)){
  print(sign.module)
  
  # Get eigengene of the significant module
  samples.eigen <- rownames(netwk$MEs)
  sign.eigen <- netwk$MEs[,sign.module]
  names(sign.eigen) <- samples.eigen
  
  # Get genes involved in the module
  sign.genes <- rownames(significant_res[significant_res$ME == sign.module, ])
  sign.genes.expr <- counts[rownames(counts) %in% sign.genes, ]
  
  # Set samples' order
  sign.genes.expr <- sign.genes.expr[,match(pds$sampleid, colnames(sign.genes.expr))]
  identical(colnames(sign.genes.expr), pds$sampleid)
  sign.eigen <- sign.eigen[match(pds$sampleid, names(sign.eigen))]
  identical(colnames(sign.genes.expr), names(sign.eigen))
  
  # Loop for mediation per gene
    for (gene in sign.genes){
    print(gene)
    x.value <- sign.eigen
    y.value <- pds
    m.value <- sign.genes.expr[rownames(sign.genes.expr) %in% gene, ]
    
    # Merge all information in one data frame
    print(identical(names(x.value), names(m.value)))
    df.values <- data.frame(
      eigen = x.value,
      gene.express = m.value,
      row.names = names(x.value)
    )
    df.values <- merge(x = df.values, y = y.value, by.x = "row.names", by.y="sampleid")
    colnames(df.values)[1] <- "sampleid"
    
    # Create mediation models: Forward direction
    # First part indirect effect model
    indirect.model1 <- glm(gene.express ~  eigen + Maternal_smoking + Baby_race + Gestational_age + Delivery_method + Infant_gender + Maternal_age + V1 + V2 + V3 + V4 + V5, data = df.values)
    # Second part direct + indirect effect model
    indirect.model2 <- glm(BW_percentile2 ~ eigen + gene.express + Maternal_smoking + Baby_race + Gestational_age + Delivery_method + Infant_gender + Maternal_age + V1 + V2 + V3 + V4 + V5, data = df.values)
    print("Mediation forward running")
    results <- mediate(indirect.model1, indirect.model2, treat='eigen', mediator='gene.express',
                       boot=TRUE, sims=1000)
    summary(results)
    
    # Add results in dataframe
    mediation_results_f <- rbind(mediation_results_f,
                                 data.frame(ME = sign.module, gene = gene, ACME.estimate = results$d0, ACME.CI_Lower = unname(results$d0.ci[1]), ACME.CI_Upper = unname(results$d0.ci[2]), ACME.p_value = results$d0.p,
                                            ADE.estimate = results$z0, ADE.CI_Lower = unname(results$z0.ci[1]), ADE.CI_Upper = unname(results$z0.ci[2]), ADE.p_value = results$z0.p, 
                                            TotalEffect.estimate = results$tau.coef, TotalEffect.CI_Lower = unname(results$tau.ci[1]), TotalEffect.CI_Upper = unname(results$tau.ci[2]), TotalEffect.p_value = results$tau.p,
                                            PropMediated.estimate = results$n0, PropMediated.CI_Lower = unname(results$n0.ci[1]), PropMediated.CI_Upper = unname(results$n0.ci[2]), PropMediated.p_value = results$n0.p))
    
    # Create mediation models: Reverse direction
    # First part indirect effect model
    indirect.model3 <- glm(gene.express ~  BW_percentile2 + Maternal_smoking + Baby_race + Gestational_age + Delivery_method + Infant_gender + Maternal_age + V1 + V2 + V3 + V4 + V5, data = df.values)
    # Second part direct + indirect effect model
    indirect.model4 <- glm(eigen ~ BW_percentile2 + gene.express + Maternal_smoking + Baby_race + Gestational_age + Delivery_method + Infant_gender + Maternal_age + V1 + V2 + V3 + V4 + V5, data = df.values)
    print("Mediation reverse running")
    results2 <- mediate(indirect.model3, indirect.model4, treat='BW_percentile2', mediator='gene.express',
                        boot=TRUE, sims=1000)
    summary(results2)
    
    # Add results in dataframe
    mediation_results_r <- rbind(mediation_results_r,
                                 data.frame(ME = sign.module, gene = gene, ACME.estimate = results2$d0, ACME.CI_Lower = unname(results2$d0.ci[1]), ACME.CI_Upper = unname(results2$d0.ci[2]), ACME.p_value = results2$d0.p,
                                            ADE.estimate = results2$z0, ADE.CI_Lower = unname(results2$z0.ci[1]), ADE.CI_Upper = unname(results2$z0.ci[2]), ADE.p_value = results2$z0.p, 
                                            TotalEffect.estimate = results2$tau.coef, TotalEffect.CI_Lower = unname(results2$tau.ci[1]), TotalEffect.CI_Upper = unname(results2$tau.ci[2]), TotalEffect.p_value = results2$tau.p,
                                            PropMediated.estimate = results2$n0, PropMediated.CI_Lower = unname(results2$n0.ci[1]), PropMediated.CI_Upper = unname(results2$n0.ci[2]), PropMediated.p_value = results2$n0.p))
    
    
  }
}


# Save R object 
saveRDS(object = mediation_results_f, file = "mediation_results_forward.RDS")
saveRDS(object = mediation_results_r, file = "mediation_results_reverse.RDS")

