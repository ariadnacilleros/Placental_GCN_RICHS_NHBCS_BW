## Placental gene co-expression networks and mediation analysis identifies pathways and genes associated with birthweight

This repository contains the code from the manuscript in review named *Placental gene co-expression networks and mediation analysis identifies pathways and genes associated with birthweight*, in which several placental gene co-expression networks have been calculated to elucidate how placental transcriptome may be shaping birthweight. 
The code consists on:
- Preprocess of the gene counts, including filtering of genes and samples, as well as batch effect correction and normalization - [preprocess_gene_counts.R](https://github.com/ariadnacilleros/Placental_GCN_RICHS_NHBCS_BW/blob/main/preprocess_gene_counts.R)
- Deconvolution of placental cell types using the reference-based method from MuSiC R package - [run_deconvolution.R](https://github.com/ariadnacilleros/Placental_GCN_RICHS_NHBCS_BW/blob/main/run_deconvolution.R)
- Obtention of gene co-expression networks by using WGCNA R package - [run_WGCNA.R](https://github.com/ariadnacilleros/Placental_GCN_RICHS_NHBCS_BW/blob/main/run_WGCNA.R)
- Calculation of the kME values from the genes belonging to each module with WGCNA R package - [get_hubgenes.R](https://github.com/ariadnacilleros/Placental_GCN_RICHS_NHBCS_BW/blob/main/get_hubgenes.R)
- Gene set enrichment of the modules with EnrichR R package and BioPlanet database - [run_enrichR.R](https://github.com/ariadnacilleros/Placental_GCN_RICHS_NHBCS_BW/blob/main/run_enrichR.R)
- Association test between modules and birthweight percentile and differential expression analysis with limma R package - [run_limma.R](https://github.com/ariadnacilleros/Placental_GCN_RICHS_NHBCS_BW/blob/main/run_limma.R)
- Mediation analysis between gene modules, gene and birthweight percentile with mediation R package - [run_mediation.R](https://github.com/ariadnacilleros/Placental_GCN_RICHS_NHBCS_BW/blob/main/run_mediation.R)
