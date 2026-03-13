### Execute WGCNA 
### 22/09/25
### Ariadna

setwd("./")
library(WGCNA)
library(CWGCNA)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read methylation data
counts <- readRDS("./merged_counts_vst_rin_673samples.RDS")
dim(counts) # 673 samples and 15,729 genes
counts[1:5,1:5]

# Read annotation
annot <- readRDS("./genes_annotation_allgenes78K.RDS")

# Remove sex chr genes 
counts <- subset(counts, rownames(counts) %in% c(annot[annot$seqnames != "chrX" & annot$seqnames != "chrY", ]$gene_id))
dim(counts) # 673 samples and 15,172 genes
rm(annot)

# Read metadata
pds <- read.csv("./final_variables_set673_pipeline240625.csv")
head(pds)

# Transform BW percentile
pds$BW_percentile2 <- pds$BW_percentile/100
pds$BW_percentile2 <- pmin(pmax(pds$BW_percentile2, 1e-6), 1 - 1e-6)
pds$BW_percentile2 <- stats::qnorm(pds$BW_percentile2)

# Make sure BW is numeric
colnames(pds)
pds <- pds %>%
  mutate(across(c(Gestational_age, BW_percentile2,  Maternal_age, V1, V2, V3, V4, V5), as.numeric))
pds <- pds %>%
  mutate(across(c(Maternal_smoking, BW_group, Baby_race, Delivery_method, Infant_gender), as.factor))

# Select only covariates of interest
colnames(pds)
pds2 <- pds[c("Sample_ID_counts", "Maternal_smoking","Baby_race","Gestational_age","BW_percentile2","Delivery_method","Infant_gender","Maternal_age",
              "V1","V2","V3","V4","V5")]
colnames(pds2)[1] <- "sampleid"

# Make sure that sample order is the same for both dataframes
identical(colnames(counts), pds2$sampleid)
# FALSE
samples_order <- match(pds2$sampleid, colnames(counts))
pds2 <- pds2[order(samples_order), ]
identical(colnames(counts), pds2$sampleid)
# TRUE

# Read annotation
annot <- readRDS("./genes_annotation_allgenes78K.RDS")
head(annot)
dim(annot) #78.932 genes

# Filter annotation by counts
annot <- annot[annot$gene_id %in% rownames(counts), ]
dim(annot) # 15.172 genes

# Remove duplicates by Gene.satble.ID
annot<-annot[!duplicated(annot$Gene.stable.ID),]
dim(annot) # 15.172 genes

# Remove duplicates by HGNC id
table(duplicated(annot$HGNC.symbol)) 
# FALSE  TRUE 
# 15082    90 
annot<-annot[!duplicated(annot$HGNC.symbol),]
dim(annot) # 15.082 genes

# Filter counts by annotated genes
counts <- counts[annot$gene_id, ]
dim(counts) #  15.082 genes

# Change counts gene names fro HGNC symbol
identical(rownames(counts), annot$gene_id)
# TRUE
counts_hgnc <- counts
rownames(counts_hgnc) <- annot$HGNC.symbol
dim(counts_hgnc) # 15,082 genes x 673 samples
counts_hgnc[1:5,1:5]

# Make sure that sample order is the same for both dataframes
identical(colnames(counts_hgnc), pds2$sampleid)
# TRUE

# Filter the top 8,000 most variable genes
counts_hgnc[1:5,1:5]
top8000 <- counts_hgnc[order(apply(counts_hgnc, FUN = sd, MARGIN = 1), decreasing = T)[1:8000], ]

# Transpose data to have rows as samples and columns as genes
input_mat = t(top8000)
input_mat[1:5,1:5]

# Choose a set of soft-thresholding powers
powers = c(1:20)

# Call the network topology analysis function 
allowWGCNAThreads() 
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  powerVector = powers,
  RsquaredCut = 0.8, 
  verbose = 5
)

# Plot power's results
par(mfrow = c(1,2));
cex1 = 0.8;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.8, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# Selected power 4

# Construct gene network
set.seed(1234)
picked_power = 4
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat, 
                          corType = "bicor",
                          
                          power = picked_power,                

                          deepSplit = TRUE,
                          pamRespectsDendro = F,
                          minModuleSize = 50,
                          maxBlockSize = 8000,
                          
                          reassignThreshold = 0,
                          mergeCutHeight = 0.2,
                          
                          saveTOMs = T,
                          saveTOMFileBase = "TOM",
                          
                          numericLabels = T,
                          verbose = 5, 
                          randomSeed = 123456789)
cor <- temp_cor 

saveRDS(netwk, "./netwk.RDS")
