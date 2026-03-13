### Preprocess counts
### 24/06/25
### Ariadna 

setwd("./")
library(data.table)
library(sva)
library(edgeR)
library(DESeq2)

#Read counts 
counts <- readRDS("./counts_richs_nhbcs.RDS")
head(counts)
dim(counts)
#78932 genes & 718 columns (rownames + samples)

#Make genes as row names 
genes<-counts$Row.names
row.names(counts)<-genes
counts <- counts[,2:ncol(counts)]

#Check dimensions
dim(counts)
#78932 genes & 717 samples

#Read variables with final list of samples
var <- read.csv("./Covariates_merged_final_240625.csv")
head(var)
dim(var) # 706 samples

# Filter samples with missing maternal smoking (18), delivery method (1) and BW percentile (5)
var <- var[!is.na(var$Maternal_smoking), ]
var <- var[!is.na(var$Delivery_method), ]
var <- var[!is.na(var$BW_percentile), ]
dim(var) # 684 samples remaining 

# Remove sex missmatch
sex_missmatch <- read.csv("./YchromSexMissmatch.csv")
var <- var[!(var$sample1 %in% sex_missmatch$x), ]
dim(var) #678 samples

#Filter counts by list of samples without NAs in important variables and duplicates 
counts <- counts[,colnames(counts) %in% var$sample1]
dim(counts) #78.932 genes & 678 samples 

# Create sample id for counts columns
var$Sample_ID_counts <- paste("X", var$Sample_ID, sep="")
head(var)

# Write sample set
write.csv(x = var, file = "samples_variables_initial_set678_pipeline240625.csv", quote = F, row.names = F)

# Change sample names from counts matrix
colnames(counts)
name_map <- setNames(var$Sample_ID_counts, var$sample1)
names(counts) <- ifelse(names(counts) %in% names(name_map), name_map[names(counts)], names(counts))
colnames(counts)

#Save counts 
saveRDS(object = counts, file = "merged_counts_allgenes_678samples.RDS")

#Get genes with 0 in all samples
table(rowSums(counts==0)==ncol(counts))
# FALSE  TRUE 
# 74379  4553
zerogenes<-row.names(counts)[rowSums(counts==0)==ncol(counts)]
length(zerogenes) # 4553

#Write file with zero genes
write.csv(x = zerogenes, file = "zeroGenes.csv", row.names = F)

#Remove zero genes
counts <- counts[!(rownames(counts) %in% zerogenes), ]

#Save counts
saveRDS(object = counts, file = "merged_counts_nozerogenes_678samples.RDS")

#Get low expressed genes (counts < 5 in 30% of the samples) - 0.3*678 = 203
table(rowSums(counts >5) > 203)
# FALSE  TRUE 
# 46354 28025  
keep<-row.names(counts)[rowSums(counts > 5) > 203]
length(keep) # 28.025 genes

#Write file with low expressed genes
write.csv(x = setdiff(rownames(counts), keep), file = "./lowExpressedGenes.csv", row.names = F)

#Remove low expressed genes
counts <- counts[rownames(counts) %in% keep, ]

#Save counts
saveRDS(object = counts, file = "merged_counts_nozerogenes_nonlowexpr_678samples.RDS")

#Check presence of inmportant placenta genes
table('ENSG00000129226.14'%in%rownames(counts)) #CD68
table('ENSG00000177575.13'%in%rownames(counts)) #CD163

# Get annotation of genes
# system("Rscript AtlasPLC_get_genes_annotation_Percentile_pipeline240625.R")

# Get cell types
# system("Rscript AtlasPLC_deconvolution_Percentile_pipeline240625.R")

#Read genes annotation
annot <- readRDS("genes_annotation_allgenes78K.RDS")
head(annot)

#Filter annotation by our genes
annot <- annot[annot$gene_id %in% rownames(counts), ]

#Get protein-coding genes 
protein_coding <- annot[annot$gene_type=="protein_coding", ]$gene_id
length(protein_coding) #15.729 protein coding genes 

#Write file with non_protein-coding genes
length(setdiff(rownames(counts), protein_coding))
write.csv(x = setdiff(rownames(counts), protein_coding), file = "nonproteincodingGenes.csv")

#Filter counts
counts <- counts[rownames(counts) %in% protein_coding, ]

#Save counts
saveRDS(object = counts, file = "merged_counts_nozerogenes_nonlowexpr_nonproteincoding_678samples.RDS")

# Remove samples with high decidual 
decidual <- read.csv(file = "samples_high_decidual_Percentile_pipeline240625.csv")
var <- var[!(var$Sample_ID_counts %in% decidual$Row.names),]
counts <- counts[,var$Sample_ID_counts]

# Write sample file with final list for this pipeline
write.csv(x = var, file = "samples_variables_final_set673_Percentile_pipeline100625.csv", quote = F, row.names = F)

#Save counts
saveRDS(object = counts, file = "merged_counts_nozerogenes_nonlowexpr_nonproteincoding_673samples.RDS")

#Check samples order
identical(colnames(counts), var$Sample_ID_counts)
#TRUE

# Format covariables to protect for correcting batch effect from biological variables
colnames(var)
covars <- var[, c(5:8,10:12)]
covars[] <- lapply(covars, function(x) {
  if (is.character(x)) factor(x) else x
})
mod <- model.matrix(~ ., data = covars)

#Correct batch effect (cohort_assay) and protect biological variables
var$cohort_assay <- paste(var$cohort, var$assay, sep = "_")
counts_assay_cohort <- ComBat_seq(counts = as.matrix(counts), batch = var$cohort_assay, covar_mod = mod)

#Save RDS object
saveRDS(object = counts_assay_cohort, file = "merged_counts_combat_assay_cohort_673samples.RDS")

# Get DESeq object
dds <- DESeqDataSetFromMatrix(countData = (counts_assay_cohort),colData = var, design = ~1)

# Normalization & transformation
dds <- estimateSizeFactors(dds)
vst <- vst(dds)

#Save RDS object
saveRDS(object = assay(vst), file = "merged_counts_vst_673samples.RDS")

# Remove RIN effect and protect biological variables
counts_vst_rin <- limma::removeBatchEffect(assay(vst), covariates = var$RIN, design = mod)

#Save RDS object
saveRDS(object = counts_vst_rin, file = "merged_counts_vst_rin_673samples.RDS")
