library(Seurat)
library(tidyverse)
library(glue)
library(SeuratWrappers)
options(scipen=5)
map(paste0(list.files()[3:15],'/Rawdata'),function(x){
  if (!dir.exists(x)) {
    dir.create(x)
  }
})
map(paste0(list.files()[3:15],'/paper'),function(x){
  if (!dir.exists(x)) {
    dir.create(x)
  }
})
map(paste0(list.files()[3:15],'/Rds'),function(x){
  if (!dir.exists(x)) {
    dir.create(x)
  }
})


islet_top_makers <- c('GCG','TTR','GC','F10',
                      'INS','IAPP','DLK1','ADCYAP1','NPY')

source('function.R')
##sample1-----
all_samples <- str_extract(list.files('./sample1/Rawdata/'),'.*\\d(?=_)') |> unique()

load_sample1 <- function(mtx,barcode,gene,project){
  sparse_mtx <- ReadMtx(mtx = mtx,
          cells = barcode,
          features = gene)
  obj <- CreateSeuratObject(sparse_mtx,project = project)
  return(obj)
}

sce_list <- map(all_samples,function(x){
  project <- x
  mtx <- paste0('./sample1/Rawdata/',x,'_matrix.mtx.gz')
  barcode <- paste0('./sample1/Rawdata/',x,'_barcodes.tsv.gz')
  gene <- paste0('./sample1/Rawdata/',x,'_genes.tsv.gz')
  sce <- load_sample1(mtx = mtx,barcode = barcode,gene = gene,project = project)
  return(sce)
})
saveRDS(sce_list,file = 'sample1/Rds/sample1_sce_raw_list.Rds')


sample1 <- reduce(sce_list,merge)
sample1$sample <- 'sample1'
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-")
VlnPlot_qc(sample1,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)
sample1 <- subset(sample1, subset = nFeature_RNA > 500 &  nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA < 50000 )
p1 <- VlnPlot_qc(sample1,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)
sample <- 'sample1'
saveplot(p1,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample1,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))
#sample1 <- readRDS(glue('{sample}/Rds/{sample}_after_qc.Rds'))


sample1 <- readRDS(file = 'sample1/Rds/sample1_after_qc.Rds')
sample1 <- run_seurat(sample1,cluster = T)

p1 <- show_umap_plot(sample1,label = F)
saveplot(p1,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)
p2 <- FeaturePlot2(sample1,gene = islet_top_makers,ncol = 3)
saveplot(p2,filenames = glue('{sample}_feature_plot2'),width = 6,height = 6)
sample1 <- sample1 %>% NormalizeData() %>% ScaleData()
p3 <- my_DotPlot(sample1,genes = islet_top_makers)
saveplot(p3,filenames = glue('{sample}_dotplot'),width = 3,height = 3)

cluster2cell <- c('Alpha','Beta','Beta','Unknown','Unknown',
                  'Unknown','Beta','Unknown','Unknown','Unknown',
                  'Unknown','Beta','Unknown','Unknown'
)
names(cluster2cell) <- levels(sample1)
sample1 <- RenameIdents(sample1, cluster2cell)
saveRDS(sample1,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))

sample1 <- readRDS(glue('{sample}/Rds/{sample}_umap_final.Rds'))
sample1_islet  <- subset(
  sample1,idents = c('Beta','Alpha')
)

sample1_islet$cell_type <- Idents(sample1_islet)
sample1_islet <- run_seurat(sample1_islet,cluster = F)
p2 <- show_umap_plot(sample1_islet,label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
p3 <- FeaturePlot2(sample1_islet,gene = islet_top_makers,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)

saveRDS(sample1_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final_temp.Rds'))  

sample <- 'sample1'
sample1_islet <- readRDS(glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))
sample1_islet$donor <- paste0('sample1_',str_extract(sample1_islet$orig.ident,'Donor_\\d+'))
sample1_beta_pseudobulk <- generate_pseudobulk(sample1_islet,sample1_islet$donor,mincell = 50)
saveRDS(sample1_beta_pseudobulk,glue('{sample}/Rds/{sample}_beta_pseudobulk.Rds'))




##sample3-----
df4 <- data.table::fread('sample3/Rawdata/GSE124742_All_human.counts.tab.gz')
sce_list3 <- load_sample2(df4,project = 'sample3')
saveRDS(sce_list3,file = 'sample3/Rds/sample3_sce_raw_list.Rds') 

sample <- 'sample3'
sample3 <- sce_list3
sample3$sample <- 'sample3'
sample3[["percent.mt"]] <- PercentageFeatureSet(sample3, pattern = "^MT-")
VlnPlot_qc(sample3,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),ncol = 3)
sample3 <- subset(sample3, subset = nFeature_RNA > 500 &  nFeature_RNA < 10000 & percent.mt < 25)
p2 <- VlnPlot_qc(sample3,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
DefaultAssay(sample3) <- 'RNA'
saveRDS(sample3,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))

sample3 <- readRDS(file = 'sample3/Rds/sample3_after_qc.Rds')
sample3 <- run_seurat(sample3,cluster = T)

p1 <- show_umap_plot(sample3,label = F)
saveplot(p1,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)
p2 <- FeaturePlot2(sample3,gene = islet_top_makers,ncol = 3)
saveplot(p2,filenames = glue('{sample}_feature_plot2'),width = 6,height = 6)
sample3 <- sample3 %>% NormalizeData() %>% ScaleData()
p3 <- my_DotPlot(sample3,genes = islet_top_makers)
saveplot(p3,filenames = glue('{sample}_dotplot'),width = 3,height = 3)
saveRDS(sample3,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))


cluster2cell <- c('Alpha','Alpha','Beta','Beta','Unknown',
                  'Unknown','Unknown','Unknown','Unknown','Alpha',
                  'Unknown','Unknown'
)
names(cluster2cell) <- levels(sample3)
sample3 <- RenameIdents(sample3, cluster2cell)

sample3_islet  <- subset(
  sample3,idents = c('Beta','Alpha')
)
sample3_islet$cell_type <- Idents(sample3_islet)
sample3_islet <- run_seurat(sample3_islet,cluster = F)
p2 <- show_umap_plot(sample3_islet,label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
p3 <- FeaturePlot2(sample3_islet,gene = islet_top_makers,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)

saveRDS(sample3_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))


sample <- 'sample3'
sample3_islet <- readRDS(glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))
sample3_islet$donor <- 'sample3_donor'
sample3_beta_pseudobulk <- generate_pseudobulk(sample3_islet,sample3_islet$donor)
saveRDS(sample3_beta_pseudobulk,glue('{sample}/Rds/{sample}_beta_pseudobulk.Rds'))


##sample4-----
sample4 <- readRDS('sample4/Rds/Shrestha_jci2021.Rds')
sample4_counts <- sample4@assays$RNA@counts
sample4_metadata <- sample4@meta.data
sce_list4 <- CreateSeuratObject(sample4_counts,meta.data = sample4_metadata,project = 'sample4')
saveRDS(sce_list4,file = 'sample4/Rds/sample4_sce_raw_list.Rds')


sample4 <- sce_list4
sample4$sample <- 'sample4'
sample <- 'sample4'
sample4[["percent.mt"]] <- PercentageFeatureSet(sample4, pattern = "^MT-")
p2 <- VlnPlot_qc(sample4,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample4,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))

sample4 <- readRDS(file = 'sample4/Rds/sample4_after_qc.Rds')
sample4 <- run_seurat(sample4,cluster = F)

p1 <- show_umap_plot(sample4,group_by = 'CellTypes',label = F)
saveplot(p1,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)
p2 <- FeaturePlot2(sample4,gene = islet_top_makers,ncol = 3)
saveplot(p2,filenames = glue('{sample}_feature_plot2'),width = 6,height = 6)
saveRDS(sample4,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))
Idents(sample4) <- sample4$CellTypes
sample4_islet  <- subset(
  sample4,idents = c('Beta','Alpha')
)
sample4_islet <- run_seurat(sample4_islet,cluster = F)
Idents(sample4_islet) <- sample4_islet$CellTypes
p2 <- show_umap_plot(sample4_islet,label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
p3 <- FeaturePlot2(sample4_islet,gene = islet_top_makers,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)
saveRDS(sample4_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))

sample <- 'sample4'
sample4_islet <- readRDS(glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))

sample4_islet$donor<- paste0('sample4_',sample4_islet$Source) 

sample4_beta_pseudobulk <- generate_pseudobulk(sample4_islet,sample4_islet$donor)
saveRDS(sample4_beta_pseudobulk,glue('{sample}/Rds/{sample}_beta_pseudobulk.Rds'))



##sample6-----
sample6_files <- list.files('sample6/Rawdata/',full.names = T)

sample6_df <- map_dfc(sample6_files,function(x){
  df <- data.table::fread(x) %>% as.data.frame()
  df <- df[1:23460,]
  rownames(df) <- df$V1
  df <- df[2]
  colnames(df) <- str_extract(x,'(?<=a/)\\w+')
  return(df)
})
sce_list6 <- CreateSeuratObject(as.sparse(sample6_df),assay = 'RNA','sample6')
saveRDS(sce_list6,file = 'sample6/Rds/sample6_sce_raw_list.Rds')


sample6 <- sce_list6
sample6$sample <- 'sample6'
sample <- 'sample6'
sample6[["percent.mt"]] <- PercentageFeatureSet(sample6, pattern = "^MT-")
VlnPlot_qc(sample6,group.by = 'sample',raster = F,features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),ncol = 3)
sample6 <- subset(sample6, subset = nFeature_RNA > 500 &  nFeature_RNA < 10000 & percent.mt < 10)
p2 <- VlnPlot_qc(sample6,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample6,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))

sample6 <- readRDS(file = 'sample6/Rds/sample6_after_qc.Rds')
sample6 <- run_seurat(sample6,cluster = T,resolution = 0.3)


p1 <- show_umap_plot(sample6,label = F)
saveplot(p1,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)
p2 <- FeaturePlot2(sample6,gene = islet_top_makers,ncol = 3)
saveplot(p2,filenames = glue('{sample}_feature_plot2'),width = 6,height = 6)
sample6 <- sample6 %>% NormalizeData() %>% ScaleData()
p3 <- my_DotPlot(sample6,genes = islet_top_makers)
saveplot(p3,filenames = glue('{sample}_dotplot'),width = 3,height = 3)
saveRDS(sample6,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))

cluster2cell <- c('Beta','Alpha','Unknown','Unknown','Alpha',
                  'Alpha','Unknown','Unknown','Unknown','Alpha',
                  'Unknown','Unknown','Unknown','Unknown'
                  )
names(cluster2cell) <- levels(sample6)
sample6 <- RenameIdents(sample6, cluster2cell)
sample6$cell_type <- Idents(sample6)
sample6_islet  <- subset(
  sample6,idents = c('Beta','Alpha')
)
sample6_islet <- run_seurat(sample6_islet,cluster = T)

p2 <- show_umap_plot(sample6_islet,group_by = 'cell_type',label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
p3 <- FeaturePlot2(sample6_islet,gene = islet_top_makers,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)

Idents(sample6_islet) <- sample6_islet$cell_type
saveRDS(sample6_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))

sample <- 'sample6'
sample6_islet <- readRDS(glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))
sample6_islet$donor <- 'sample6_donor'
sample6_beta_pseudobulk <- generate_pseudobulk(sample6_islet,'sample6_donor')
saveRDS(sample6_beta_pseudobulk,glue('{sample}/Rds/{sample}_beta_pseudobulk.Rds'))



##sample10-----
sample10_features <- data.table::fread('sample10/Rawdata/sample10_mtx/sample10_mtx_features.txt.gz') 
sample10_features <- sample10_features %>% left_join(ensembl2symbol,c('ensembl_id' = 'gene_id')) 
sample10_features$gene_name[which(is.na(sample10_features$gene_name))] <- sample10_features$ensembl_id[which(is.na(sample10_features$gene_name))]
data.table::fwrite(sample10_features,file = 'sample10/Rawdata/sample10_mtx/sample10_mtx_features_symbol.txt',sep = '\t')

sample10 <- ReadMtx(
  mtx = 'sample10/Rawdata/sample10_mtx/sample10_mtx.mtx',
  cells = 'sample10/Rawdata/sample10_mtx/sample10_mtx_metadata.txt.gz',skip.cell = 1,
  features = 'sample10/Rawdata/sample10_mtx/sample10_mtx_features_symbol.txt',
  feature.column = 5,skip.feature = 1
)

metadata <- read.delim("sample10/Rawdata/sample10_mtx/sample10_mtx_metadata.txt.gz", row.names=1)
sample10 <- CreateSeuratObject(sample10,meta.data = metadata)
saveRDS(sample10,'sample10/Rds/sample10_sce_raw_list.Rds')

sample10$sample <- 'sample10'
sample <- 'sample10'
sample10[["percent.mt"]] <- PercentageFeatureSet(sample10, pattern = "^MT-")
p2 <- VlnPlot_qc(sample10,group.by = 'sample',raster = F,features = c(glue("nFeature_RNA"),glue("nCount_RNA"),'percent.mt'),ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample10,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))
sample10 <- readRDS(file = 'sample10/Rds/sample10_after_qc.Rds')
sample10 <- run_seurat(sample10,cluster = T)

p1 <- show_umap_plot(sample10,label = F)
saveplot(p1,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)
p2 <- FeaturePlot2(sample10,gene = islet_top_makers,ncol = 3)
saveplot(p2,filenames = glue('{sample}_feature_plot2'),width = 6,height = 6)
sample2 <- sample2 %>% NormalizeData() %>% ScaleData()
p3 <- my_DotPlot(sample2,genes = islet_top_makers)
saveplot(p3,filenames = glue('{sample}_dotplot'),width = 3,height = 3)
saveRDS(sample10,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))


sample10_islet  <- subset(
  sample10,idents = c(4,9),invert =T 
)
sample10_islet <- run_seurat(sample10_islet,cluster = T)

p2 <- show_umap_plot(sample10_islet,label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
p3 <- FeaturePlot2(sample10_islet,gene = islet_top_makers,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)

cluster2cell <- c('Alpha','Beta','Alpha','Beta','Alpha','Beta','Beta','Beta')
names(cluster2cell) <- levels(sample10_islet)
sample10_islet <- RenameIdents(sample10_islet, cluster2cell)
p2 <- show_umap_plot(sample10_islet,label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
saveRDS(sample10_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))


sample <- 'sample10'
sample10_islet <- readRDS(glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))
sample10_islet$donor <- paste0(sample10_islet$sample,'_',sample10_islet$id)
sample10_beta_pseudobulk <- generate_pseudobulk(sample10_islet,sample10_islet$donor)
saveRDS(sample10_beta_pseudobulk,glue('{sample}/Rds/{sample}_beta_pseudobulk.Rds'))




##sample12----
sample12 <- readRDS('sample12/Rawdata/hpap_islet_scRNAseq.rds')
sample12$disease <- sample12$`Diabetes Status` 
sample12 <- subset(sample12,disease %in% c('ND','T2D'))
saveRDS(sample12,'sample12/Rds/sample12_sce_raw_list.Rds')

sample12$sample <- 'sample12'
sample <- 'sample12'
sample12[["percent.mt"]] <- PercentageFeatureSet(sample12, pattern = "^MT-")
sample12$cell_type <- sample12$`Cell Type Grouped`
p2 <- VlnPlot_qc(sample12,group.by = 'sample',raster = F,features = c(glue("nFeature_RNA"),glue("nCount_RNA"),'percent.mt'),ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample12,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))

sample12 <- readRDS('sample12/Rds/sample12_after_qc.Rds')
Idents(sample12) <- sample12$cell_type
islet_makers <- FindAllMarkers(sample12,assay = 'RNA',logfc.threshold = 1)
saveRDS(islet_makers,file = 'sample12/islet_makers.Rds')

sample12 <- run_seurat(sample12)


p1 <- show_umap_plot(sample12,group_by = 'cell_type',label = F)
saveplot(p1,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)
Idents(sample12) <- sample12$cell_type
saveRDS(sample12,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))  

sample12_islet  <- subset(
  sample12,idents = c('Alpha','Beta')
)
sample12_islet <- run_seurat(sample12_islet)

sample12_islet$donor_id <- paste0('sample12','_',str_extract(colnames(sample12_islet),'.*(?=_)'))
sample12_islet$donor <- paste0(sample12_islet$donor_id,'_',sample12_islet$disease)

saveRDS(sample12_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'),compress = FALSE)


library(scCustomize)
p8 <- Cell_Highlight_Plot(sample12_islet,cells_highlight = list(
  'sample12_HPAP-022_ND' = names(which(sample12_islet$donor == 'sample12_HPAP-022_ND'))
)
)
p9 <- Cell_Highlight_Plot(sample12_islet,cells_highlight = list(
  'sample12_HPAP-109_T2D' = names(which(sample12_islet$donor == 'sample12_HPAP-109_T2D'))
)
)
saveplot(p8,filenames = glue('{sample}_sample1_Donor_6_umap'),width = 6,height = 6)
saveplot(p9,filenames = glue('{sample}_sample1_Donor_11_umap'),width = 6,height = 6)





library(COSG)
Idents(sample12_islet) <- sample12_islet$cell_type
marker_cosg <- cosg(
  sample12_islet,
  assay='RNA',
  mu=100,
  remove_lowly_expressed=TRUE,
  expressed_pct=0.8,
  n_genes_user=10)
marker_cosg_marker <- marker_cosg$names %>%
  pivot_longer(cols = everything(),
               names_to = 'cell',
               values_to = 'gene') 



p2 <- show_umap_plot(sample12_islet,group_by = 'cell_type',label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)

islet_top_makers2 <- marker_cosg_marker$gene[10:18]
p3 <- FeaturePlot2(sample12_islet,gene = islet_top_makers,ncol = 3)
p4 <- FeaturePlot2(sample12_islet,gene = islet_top_makers2,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)
saveplot(p4,filenames = glue('{sample}_feature_plot2'),width = 6,height = 6)

##sample14-----

# sample14 <- ReadMtx(
#   mtx = 'sample14/Rawdata/sample14_mtx/sample14_mtx.mtx',
#   cells = 'sample14/Rawdata/sample14_mtx/sample14_mtx_metadata.txt.gz',skip.cell = 1,
#   features = 'sample14/Rawdata/sample14_mtx/sample14_mtx_features.txt.gz',
#   feature.column = 1,skip.feature = 1
# )
# metadata <- read.delim("sample14/Rawdata/sample14_mtx/sample14_mtx_metadata.txt.gz", row.names=1)
# sample14 <- CreateSeuratObject(sample14,meta.data = metadata)
# saveRDS(sample14,'sample14/Rds/sample14_sce_raw_list.Rds')

sample14 <- data.table::fread('sample14/Rawdata/E-MTAB-5061/pancreas_counts_3514sc.txt')
colnames(sample14)[1] <- 'gene'
sample14 <- sample14[which(!duplicated(sample14[[1]])),]
sample14 <- sample14[which(!str_starts(sample14[[1]],'ERCC')),]
sample14 <- as.data.frame(sample14)
rownames(sample14) <- sample14[[1]]
sample14 <- sample14[,-1]
sample14 <- as.sparse(sample14)
sample14_seurat <- CreateSeuratObject(sample14,assay = 'RNA',project = 'sample14')


sample14_2 <- sample14 %>% group_by(gene) %>% summarise(across(everything(),sum))
sample14_2 <- sample14_2[which(!str_starts(sample14_2[[1]],'ERCC')),]
sample14_2 <- as.data.frame(sample14_2)
rownames(sample14_2) <- sample14_2[[1]]
sample14_2 <- sample14_2[,-1]
sample14_2 <- as.sparse(sample14_2)
sample14_seurat <- CreateSeuratObject(sample14_2,assay = 'RNA',project = 'sample14')

metadata <- data.table::fread('sample14/Rawdata/E-MTAB-5061/E-MTAB-5061.sdrf.txt') |> as.data.frame()
cell_id <- intersect(colnames(sample14_2),metadata$`Source Name`)
metadata <- metadata[metadata$`Source Name` %in% cell_id,] 
sample14_seurat <- subset(sample14_seurat,cells = cell_id)
rownames(metadata) <- metadata$`Source Name`
metadata <- metadata[colnames(sample14_seurat),]
identical(colnames(sample14_seurat),rownames(metadata))

sample14_seurat <- AddMetaData(sample14_seurat,metadata = metadata)
saveRDS(sample14_seurat,'sample14/Rds/sample14_sce_raw_list.Rds')

sample14_seurat$sample <- 'sample14'
sample <- 'sample14'
sample14_seurat[["percent.mt"]] <- PercentageFeatureSet(sample14_seurat, pattern = "^MT-")
VlnPlot_qc(sample14_seurat,group.by = 'sample',raster = F,features = c(glue("nFeature_RNA"),glue("nCount_RNA"),'percent.mt'),ncol = 3)
sample14_seurat <- subset(sample14_seurat, subset = nFeature_RNA > 500 & percent.mt < 10 )
sample14_seurat <- subset(sample14_seurat, subset = nFeature_RNA > 500 & percent.mt < 10 &  
                          nCount_RNA < as.numeric(str_extract(colnames(catable(sample14_seurat$nCount_RNA))[6],'\\d+')) & 
                          nCount_RNA > as.numeric(str_extract(colnames(catable(sample14_seurat$nCount_RNA))[1],'\\d+')) )

p2 <- VlnPlot_qc(sample14_seurat,group.by = 'sample',raster = F,features = c(glue("nFeature_RNA"),glue("nCount_RNA"),'percent.mt'),ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample14_seurat,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))


sample14_seurat <- run_seurat(sample14_seurat,cluster = F)
sample14_seurat$donor <- sample14_seurat$`Characteristics [individual]`
library(harmony)
sample14_seurat <- NormalizeData(sample14_seurat) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
sample14_seurat <- sample14_seurat %>% RunHarmony(group.by.vars = "donor") %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% FindNeighbors()

sample14_seurat <- sample14_seurat %>% FindClusters(resolution = 0.8)

p1 <- my_DotPlot(sample14_seurat,genes = all_makers,scale = T)
p2 <- show_umap_plot(sample14_seurat,label = F)
p1+p2
FeaturePlot2(sample14_seurat,gene = islet_top_makers,ncol = 3)


cluster2cell <- c('Unknown','Alpha','Alpha','Unknown','Unknown',
                  'Alpha','Beta','Unknown','Unknown','Unknown',
                  'Beta','Unknown','Unknown','Unknown'
)
names(cluster2cell) <- levels(sample14_seurat)
sample14_seurat <- RenameIdents(sample14_seurat, cluster2cell)
saveRDS(sample14_seurat,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))


#sample14_seurat <- readRDS(glue('{sample}/Rds/{sample}_umap_final.Rds'))

sample14_islet  <- subset(
  sample14_seurat,idents = c('Alpha','Beta')
)
sample14_islet$cell_type <- Idents(sample14_islet)


library(harmony)
sample14_islet <- NormalizeData(sample14_islet) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
sample14_islet <- sample14_islet %>% RunHarmony(group.by.vars = "donor") %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% FindNeighbors()
saveRDS(sample14_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))



p2 <- show_umap_plot(sample14_islet,label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
p3 <- FeaturePlot2(sample14_islet,gene = islet_top_makers,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)
saveRDS(sample14_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))

islet_top_makers <- c(
  'GCG','TTR','GC','F10',
  'INS','IAPP','DLK1','ADCYAP1','NPY'
)
other_cell_makers <- c(
  'SST','PPY','CFTR','KRT19','REG1A',
  'PRSS1','PECAM1','COL1A1','CD68','CD74'
)
all_makers <- c(islet_top_makers,other_cell_makers)

p6 <- my_DotPlot(sample14_islet,group_by = 'cell_type',genes = all_makers,scale = F)
saveplot(p6,filenames = glue('{sample}_dotplot'),width = 6,height = 3)

p4 <- FeaturePlot2(sample14_islet,gene = islet_top_makers,ncol = 3)
saveplot(p4,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)

p5 <- FeaturePlot2(sample14_islet,gene = other_cell_makers[1:9],ncol = 3)
saveplot(p5,filenames = glue('{sample}_feature_plot_othermaker'),width = 6,height = 6)



##sample16----
sample <- 'sample16'
sample16_count <- data.table::fread('sample16/Rawdata/counts.txt')
sample16_count <- sample16_count[,c(-2:-6)]

colnames(sample16_count)[-1] <- str_extract(colnames(sample16_count)[-1],'SRR\\d+')
sample16_count <- as.data.frame(sample16_count)
rownames(sample16_count) <- sample16_count[[1]]
sample16_count <- sample16_count[,-1]
sample16_count <- as.sparse(sample16_count)

sample16_metadata <- data.table::fread('sample16/Rawdata/sample16_metadata.txt')
runinfo <- data.table::fread('sample16/Paper/SraRunInfo.csv') %>% select(Run,SampleName)
sample16_metadata <- sample16_metadata %>% left_join(runinfo,by = c('GSM' = 'SampleName'))
sample16_metadata <- as.data.frame(sample16_metadata)
rownames(sample16_metadata) <- sample16_metadata$Run
sample16_metadata <- sample16_metadata[colnames(sample16_count),]

identical(colnames(sample16_count),rownames(sample16_metadata))
sample16 <- CreateSeuratObject(sample16_count,meta.data = sample16_metadata,assay = 'RNA',project = 'sample16')
saveRDS(sample16,'sample16/Rds/sample16_sce_raw_list.Rds')

sample16$sample <- 'sample16'
sample <- 'sample16'
sample16[["percent.mt"]] <- PercentageFeatureSet(sample16, pattern = "^MT-")
VlnPlot_qc(sample16,group.by = 'sample',raster = F,features = c(glue("nFeature_RNA"),glue("nCount_RNA"),'percent.mt'),ncol = 3)
sample16 <- subset(sample16, subset = nFeature_RNA > 500 & 
                            nCount_RNA < as.numeric(str_extract(colnames(catable(sample16$nCount_RNA))[6],'\\d+')) & 
                            nCount_RNA > as.numeric(str_extract(colnames(catable(sample16$nCount_RNA))[1],'\\d+')) )

p2 <- VlnPlot_qc(sample16,group.by = 'sample',raster = F,features = c(glue("nFeature_RNA"),glue("nCount_RNA"),'percent.mt'),ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample16,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))



sample16 <- run_seurat(sample16,cluster = F)
sample16$donor <- paste0('sample16_',)
sample16$donor_id <- str_replace_all(sample16$donor_id,' ','_')

sample16$donor <- paste0('sample16_',sample16$donor_id,'_',
                              map_chr(sample16$donor_id,function(x){
                                if (str_detect(x,'Non')) {
                                  return('ND')
                                }else{
                                  return('T2D')
                                }
                              })) 

sample16 <- sample16 %>% FindNeighbors()

Idents(sample16) <- sample16$cell_type
sample16 <- sample16 %>% FindSubCluster(graph.name = 'SCT_snn',cluster = 'beta',resolution = 2)
DimPlot(sample16,group.by = 'sub.cluster',label = T)
sample16 <- subset(sample16,sub.cluster %in% c('beta_11','beta10'),invert = T)


sample16 <- sample16 %>% FindSubCluster(graph.name = 'SCT_snn',cluster = 'alpha',resolution = 2)
DimPlot(sample16,group.by = 'cell_type',label = T)
sample16 <- subset(sample16,sub.cluster %in% c('alpha_11','alpha_12'),invert = T)



cluster2cell <- c('Beta','Alpha','PP','delta'
)
names(cluster2cell) <- levels(sample16)
sample16 <- RenameIdents(sample16, cluster2cell)
saveRDS(sample16,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))


sample16_islet  <- subset(
  sample16,idents = c('Alpha','Beta')
)
sample16_islet$cell_type <- Idents(sample16_islet)
sample16_islet <- run_seurat(sample16_islet,cluster = F)

saveRDS(sample16_islet,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))
##sample17----
sample <- 'sample17'

load_sample17 <- function(file,project){
    df <- data.table::fread(file)
    df <- as.data.frame(df)
    rownames(df) <- df[[1]]
    df <- df[,-1]
    df <- as.sparse(df)
    colnames(df) <- paste0(str_extract(file,'(?<=\\d_).*(?=\\.dr)'),'#',colnames(df))
    sce <- CreateSeuratObject(df,assay = 'RNA',project)
    sce$donor <- paste0(project,'_',str_extract(file,'(?<=\\d_).*(?=\\.dr)'))
    return(sce)
}

sample17_list <- map(list.files('sample17/raw_data/',full.names = T),load_sample17,project = 'sample17',.progress = T)
sample17 <- purrr::reduce(sample17_list,merge)
saveRDS(sample17,file = 'sample17/Rds/sample17_sce_raw_list.Rds')

sample17$sample <- 'sample17'
sample <- 'sample17'
sample17[["percent.mt"]] <- PercentageFeatureSet(sample17, pattern = "^MT-")
sample17$percent.mt[which(is.na(sample17$percent.mt))] <- 0

VlnPlot_qc(sample17,group.by = 'donor',raster = T,features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),ncol = 3)

sample17 <- subset(sample17, subset = nFeature_RNA > 200 &  nFeature_RNA < 4000 & percent.mt < 10)
sample17 <- subset(sample17, subset = nFeature_RNA > 200 &  nFeature_RNA < 4000 & percent.mt < 10 &
                            nCount_RNA < as.numeric(str_extract(colnames(catable(sample17$nCount_RNA))[6],'\\d+')) &
                            nCount_RNA > as.numeric(str_extract(colnames(catable(sample17$nCount_RNA))[1],'\\d+')) )


p2 <- VlnPlot_qc(sample17,group.by = 'donor',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc'),width = 6,height = 4)
saveRDS(sample17,file = glue('{sample}/Rds/{sample}_after_qc.Rds'))


sample17 <- run_seurat(sample17,cluster = T,resolution = 0.3)

p1 <- show_umap_plot(sample17,label.size = 10,label = F)
saveplot(p1,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)
p2 <- FeaturePlot2(sample17,gene = islet_top_makers,ncol = 3)
saveplot(p2,filenames = glue('{sample}_feature_plot2'),width = 6,height = 6)
sample17 <- sample17 %>% NormalizeData() %>% ScaleData()
p3 <- my_DotPlot(sample17,genes = islet_top_makers)
saveplot(p3,filenames = glue('{sample}_dotplot'),width = 3,height = 3)
saveRDS(sample17,file =  glue('{sample}/Rds/{sample}_umap_final.Rds'))

sample17_alpha_beta <- subset(sample17,idents = c(0,1,4,5,9,11))
sample17_alpha_beta <- run_seurat(sample17_alpha_beta,cluster = T,resolution = 0.3)
sample17_alpha_beta <- sample17_alpha_beta %>% NormalizeData() %>% ScaleData()
sample17_alpha_beta <- subset(sample17_alpha_beta,idents = 9,invert = T)
sample17_alpha_beta <- SCTransform(sample17_alpha_beta) %>% RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:20,umap.method = 'umap-learn',metric = 'correlation',verbose = FALSE)   %>%
    FindNeighbors() %>% 
    FindClusters(resolution = 1)

p1 <- show_umap_plot(sample17_alpha_beta,label.size = 10,label = F)
p3 <- my_DotPlot(sample17_alpha_beta,genes = islet_top_makers)
p1+p3
cluster2cell <- c('Alpha','Alpha','Beta','Alpha','Alpha',
                  'Alpha','Alpha','Beta','Beta','Alpha',
                  'Beta','Beta','Unknown','Unknown','Beta',
                  'Beta','Beta','Beta','Beta'
)
names(cluster2cell) <- levels(sample17_alpha_beta)
sample17_alpha_beta <- RenameIdents(sample17_alpha_beta, cluster2cell)
sample17_alpha_beta$cell_type <- Idents(sample17_alpha_beta)

sample17_alpha_beta  <- subset(
    sample17_alpha_beta,idents = c('Beta','Alpha')
)
sample17_alpha_beta <- run_seurat(sample17_alpha_beta,cluster = T)

p2 <- show_umap_plot(sample17_alpha_beta,group_by = 'cell_type',label = F)
saveplot(p2,filenames = glue('{sample}_alpha_beta_umap'),width = 4,height = 4)
p3 <- FeaturePlot2(sample17_alpha_beta,gene = islet_top_makers,ncol = 3)
saveplot(p3,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)
Idents(sample17_alpha_beta) <- sample17_alpha_beta$cell_type

saveRDS(sample17_alpha_beta,file =  glue('{sample}/Rds/{sample}_alpha_beta_final.Rds'))

DefaultAssay(sample17_alpha_beta) <- 'RNA'

sample17_islet <- NormalizeData(sample17_alpha_beta) %>%  FindVariableFeatures()
sample17_islet <- RunFastMNN(object.list = SplitObject(sample17_islet, split.by = "donor")) %>%
    RunUMAP(reduction = "mnn", dims = 1:30) %>%
    FindNeighbors(reduction = "mnn", dims = 1:30)

sample17_islet$status <- if_else(str_detect(sample17_islet$donor,'HT'),'ND','T2D')
p5 <- show_umap_plot(sample17_alpha_beta,group_by = 'cell_type',label = F)

sample17_islet <- sample17_islet %>% RunUMAP(reduction = "mnn", dims = 1:30) 
sample17_islet <- sample17_islet %>%  FindNeighbors(reduction = "mnn", dims = 1:30) %>% FindClusters(resolution = 0.3)
show_umap_plot(sample17_islet,label = F)
cluster2cell <- c('Alpha','Alpha','Beta','Beta','Beta')
names(cluster2cell) <- levels(sample17_islet)
sample17_islet <- RenameIdents(sample17_islet, cluster2cell)
sample17_islet$cell_type <- Idents(sample17_islet)
islet_top_makers <- c(
    'GCG','TTR','GC','F10',
    'INS','IAPP','DLK1','ADCYAP1','NPY'
)
other_cell_makers <- c(
    'SST','PPY','CFTR','KRT19','REG1A',
    'PRSS1','PECAM1','COL1A1','CD68','CD74'
)
all_makers <- c(islet_top_makers,other_cell_makers)
p7 <- show_umap_plot(sample17_islet,label = F)
saveplot(p7,filenames = glue('{sample}_alpha_beta_umap_intergation'),width = 4,height = 4)

p8 <- show_umap_plot(sample17_islet,group_by = 'donor',label = T)
saveplot(p8,filenames = glue('{sample}_alpha_beta_umap_intergation_donor'),width = 4,height = 4)

p6 <- my_DotPlot(sample17_islet,group_by = 'cell_type',genes = all_makers,scale = F)
saveplot(p6,filenames = glue('{sample}_dotplot'),width = 6,height = 3)

p4 <- FeaturePlot2(sample17_islet,gene = islet_top_makers,ncol = 3)
saveplot(p4,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)

p5 <- FeaturePlot2(sample17_islet,gene = other_cell_makers[1:9],ncol = 3)
saveplot(p5,filenames = glue('{sample}_feature_plot_othermaker'),width = 6,height = 6)


saveRDS(sample17_islet,glue('sample17/Rds/{sample}_alpha_beta_final.Rds'))


sample17_beta_cell <- subset(sample17_islet,cell_type == 'Beta')
Idents(sample17_beta_cell) <- sample17_beta_cell$status

sample17_islet_DEGs <- FindMarkers(sample17_beta_cell,ident.1 = 'ND',ident.2 = 'T2D') %>%
    rownames_to_column('gene')

