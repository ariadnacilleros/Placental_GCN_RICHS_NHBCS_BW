### Perform deconvolution
### 24/06/25
### Ariadna 

setwd("./")
library(data.table)
library(sva)
library(MuSiC)
library(dplyr)
library(SingleCellExperiment)

# Read counts 
counts <- readRDS("merged_counts_nozerogenes_nonlowexpr_678samples.RDS")
counts[1:5,1:5]
dim(counts) #28.025 genes & 678 samples

# Read annotation
annot <- readRDS("genes_annotation_allgenes78K.RDS")
head(annot)
dim(annot) #78.932 genes

# Filter annotation by counts
annot <- annot[annot$gene_id %in% rownames(counts), ]
dim(annot) # 28.025 genes

# Remove duplicates by Gene.satble.ID
annot<-annot[!duplicated(annot$Gene.stable.ID),]
dim(annot) # 28.025 genes

# Remove duplicates by HGNC id
table(duplicated(annot$HGNC.symbol)) 
# FALSE  TRUE 
# 20287  7738 
annot<-annot[!duplicated(annot$HGNC.symbol),]
dim(annot) # 20.287 genes

# Read reference 
plcscSet<-readRDS("./Pique2019_LabeledCells_eSet1.rds")
dim(plcscSet) #19194 genes x 79719  cells  
#remove not labeled cells
plcscSet<-plcscSet[, !grepl('UNK', colnames(plcscSet))]
dim(plcscSet) #19194 genes x 79520  cells 
table(rownames(plcscSet)%in%annot$HGNC.symbol)
# FALSE  TRUE 
# 3098 16096 
plcscSet<-plcscSet[rownames(plcscSet)%in%annot$HGNC.symbol,]
dim(plcscSet) 
# Features  Samples 
# 16097    79520 
gc()

# Filter annotation and counts by reference common genes
annot <- annot[annot$HGNC.symbol %in% rownames(plcscSet),]
dim(annot) # 16096 genes
counts <- counts[annot$gene_id, ]
dim(counts) #  16096 genes

# Change counts gene names fro HGNC symbol
identical(rownames(counts), annot$gene_id)
# TRUE
counts_hgnc <- counts
rownames(counts_hgnc) <- annot$HGNC.symbol
dim(counts_hgnc) # 16.096 genes x 678 samples
counts_hgnc[1:5,1:5]

# Read variables
var <- read.csv("samples_variables_initial_set678_pipeline240625.csv")
head(var)

# Match samples order
identical(colnames(counts_hgnc), var$Sample_ID_counts)
#FALSE
samples_order <- match(var$Sample_ID_counts, colnames(counts_hgnc))
var <- var[order(samples_order), ]
identical(colnames(counts_hgnc), var$Sample_ID_counts)
#TRUE

# Save R object
saveRDS(object = counts_hgnc, file = "merged_counts_hgnc_deconvolution.RDS")

# Format covariables to protect for correcting batch effect from cohort variable
colnames(var)
covars <- var[, c(5:8,10:12)]
covars[] <- lapply(covars, function(x) {
  if (is.character(x)) factor(x) else x
})
mod <- model.matrix(~ ., data = covars)

# Adjust by assay (s/p)
var$cohort_assay <- paste(var$cohort, var$assay, sep = "_")
counts_hgnc_assay_cohort <- ComBat_seq(counts = as.matrix(counts_hgnc), batch = var$cohort_assay, covar_mod = mod)
rm(covars)
rm(mod)
rm(counts)
rm(counts_hgnc)
rm(annot)

# Save RDS object
saveRDS(object = counts_hgnc_assay_cohort, file = "merged_counts_hgnc_combat_assay_cohort_deconvolution.RDS")

# Perform deconvolution
plcscSet_SCE <- SingleCellExperiment(
  assays = list(counts = exprs(plcscSet)),
  colData = pData(plcscSet),
  rowData = fData(plcscSet)
)
countsresmu1 <- music_prop(
  bulk.mtx = counts_hgnc_assay_cohort,    # bulk matrix (genes x samples)
  sc.sce = plcscSet_SCE,      # single-cell ExpressionSet
  clusters = 'cellType',   # which column in sc.eset@phenoData identifies clusters
  samples = 'SubjectName'  # which column in sc.eset@phenoData identifies samples
)
gc()
countsresmu1.w <- as.data.frame(countsresmu1$Est.prop.weighted)
summary(countsresmu1.w$STB)
colnames(countsresmu1.w)[colnames(countsresmu1.w)=="Hoffbauer"]<-"Hofbauer"

# Save R object
saveRDS(object = countsresmu1, file = "music_prop_240625.RDS")
