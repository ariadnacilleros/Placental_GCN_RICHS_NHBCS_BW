### Get Hub genes WGCNA 
### 23/09/25
### Ariadna

setwd("./")
library(WGCNA)
library(openxlsx)

# Read gene network
netwk <- readRDS("./netwk.RDS")
eigengenes <- netwk$MEs
eigengenes[1:5,1:5]

# Read counts
counts <- readRDS("./counts_genenetwork_top8000.RDS")
counts[1:5,1:5]
counts <- t(counts)
counts[1:5,1:5]

# Calculate KME
kmes<-signedKME(datExpr = counts, datME = eigengenes)
kmes[1:5,1:5]

# Order rows by KME in each module
kmes_list <- list()
for (i in colnames(kmes)){
  kmes_list[[i]] <- as.data.frame(kmes[,colnames(kmes) == i], row.names = rownames(kmes))
}
kmes_list_ordered <- lapply(kmes_list, function(df) {
  colname <- names(df)[1]  # get the column name
  df[order(abs(df[[colname]]), decreasing = T), , drop = FALSE]
})

saveRDS(object = kmes_list_ordered,file = "hubgenes.RDS")

# Write xlsx file with sheet per module
kmes_list_ordered <- kmes_list_ordered[match(paste("kME", 0:23, sep=""), names(kmes_list_ordered))]
wb <- createWorkbook()
for (module in names(kmes_list_ordered)) {
  addWorksheet(wb, module)
  writeData(wb, module, kmes_list_ordered[[module]], colNames = T, rowNames = T)
}
saveWorkbook(wb, "hubgenes.xlsx", overwrite = TRUE)
