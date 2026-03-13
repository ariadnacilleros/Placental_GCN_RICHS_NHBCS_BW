## Placental gene co-expression networks and mediation analysis identifies pathways and genes associated with birthweight

This repository contains the code from the manuscript in review named *Placental gene co-expression networks and mediation analysis identifies pathways and genes associated with birthweight*, in which several placental gene co-expression networks have been calculated to elucidate how placental transcriptome may be shaping birthweight. 
The code consists on:
- Preprocess of the gene counts, including filtering of genes and samples, as well as batch effect correction and normalization.
- Deconvolution of placental cell types using the reference-based method from MuSiC R package. 
- Obtention of gene co-expression networks by using WGCNA R package.
- Calculation of the kME values from the genes belonging to each module with WGCNA R package.
- Gene set enrichment of the modules with EnrichR R package and BioPlanet database.
- Association test between modules and birthweight percentile and differential expression analysis with limma R package.
- Mediation analysis between gene modules, gene and birthweight percentile with mediation R package.
