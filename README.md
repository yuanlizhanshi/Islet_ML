This the reproducible code for the islet paper.

# Preprocess scRNAseq data

This folder contains the code and data used for the sRNA-seq analysis 

Processed beta cell scRNA-seq dataset is stored in [cowtransfer](https://drctdb.cowtransfer.com/s/ffd1b19e3a6d41)

After ML quality controlled, beta cell scRNA-seq dataset is stored in [cowtransfer](https://drctdb.cowtransfer.com/s/13a4c9b21e8741)

*.R, r code for downstream analysis of scRNA-seq, including quality control,clustering, cell annotation, 

# Pre-trained human islets atlas

The well-annotated scRNA-seq datsets was trained by scVI and scANVI to learn cell representation.

The trained h5ad file is stored in [cowtransfer](https://drctdb.cowtransfer.com/s/9adb968646324d)

# Preprocess scATAC-seq data

This folder contains the code and data used for the scATAC-seq analysis (multiome + scATAC-seq).

Processed beta cell scATAC-seq data is stored in [cowtransfer](https://drctdb.cowtransfer.com/s/e6a91494db9346)

*.smk data, snakemake scripts for cellranger 
*.R, r code for downstream analysis of scATAC-seq, including quality control, clustering, cell annotation and peak calling


# Preprocess ChIP-seq data and HiC data

For ChIP-seq data, we only run the basic upstream analysis, such quality control and mapping. The [bam file of H3H27ac](https://drctdb.cowtransfer.com/s/8590c30adce14e) modification will used for ABC model input.

For HiC data, we could run the basic upstream analysis, such quality control, mapping (This workflow is reference from [Renlab](https://github.com/ren-lab/hic-pipeline) ). The [hic file](https://drctdb.cowtransfer.com/s/f85edaafd12d46) will used for ABC model input.


# Infer and refine GRN form single cell multiomics data

First, split the scATAC bam file based on the cell type information. Then generate the gene expression (TPM) from multiome dataset from [Wang et al. 2023](https://www.nature.com/articles/s41588-023-01397-9)

After prepare the all data, using code in run_ABC.Rmd to generate enhancer promoter interactions with [ABC model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction)

Then using the [motifmachr](https://github.com/GreenleafLab/motifmatchr) to assign the TF to the corresponding cis-regulatory element


# Using XGboost to remove donors with low predicted accuracy rates

Due to the heterogeneity of human data, we found some of the donors with the discrepancy gene expression profile, which exhibit the extremely low predicted accuracy rates (15%). We decided run iterative XGboost to remove these donors until no donor with low predicted accuracy rates, then calculate the differently expressed genes.