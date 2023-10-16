library(Seurat)
library(tidyverse)
library(glue)
library(SeuratWrappers)
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
islet_top_makers <- c(
  'GCG','TTR','GC','F10',
  'INS','IAPP','DLK1','ADCYAP1','NPY'
)
other_cell_makers <- c(
  'SST','PPY','CFTR','KRT19','REG1A',
  'PRSS1','COL1A1','CD68','CD74'
)
all_makers <- c(islet_top_makers,other_cell_makers)
gene_length <- readRDS('../../support_data/human_ensembl2symbol_with_length.Rds')
get_tpm <- function(count_mtx){
  tpm <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  count_mtx$gene_name <- rownames(count_mtx)
  count_mtx <- count_mtx %>% left_join(gene_length[3:4],by = c('gene_name')) %>% 
    filter(!is.na(length))
  count_tpm <- tpm(count_mtx[,str_which(colnames(count_mtx),'Biobank')],count_mtx$length)
  rownames(count_tpm) <- count_mtx$gene_name
  return(count_tpm)
}
source('../function.R')
# seurat_obj_list <- readRDS('Rds/biobank_multiome_scRNA_seurat_obj_list.Rds')
# seurat_obj <- purrr::reduce(seurat_obj_list,merge)
seurat_obj <- readRDS('Rds/biobank_multiome_scRNA_seurat_obj.Rds')
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj$sample <- 'biobank_multiome'
sample <- 'biobank_multiome'
VlnPlot_qc(seurat_obj,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),raster = FALSE, ncol = 3)
VlnPlot_qc(seurat_obj,features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),raster = FALSE, ncol = 3)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 &  nFeature_RNA < 5000 & percent.mt < 10 &  
                          nCount_RNA < as.numeric(str_extract(colnames(catable(seurat_obj$nCount_RNA))[6],'\\d+')) & 
                          nCount_RNA > as.numeric(str_extract(colnames(catable(seurat_obj$nCount_RNA))[1],'\\d+')) )
#run 2 
p1 <- VlnPlot_qc(seurat_obj,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),raster = FALSE,ncol = 3)
saveplot(p1,filenames = glue('{sample}_final_qc'),width = 4,height = 3)

p2 <- VlnPlot_qc(seurat_obj,features = c("nFeature_RNA", "nCount_RNA",'percent.mt'),raster = FALSE, ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc_donor'),width = 12,height = 3)

seurat_obj <- NormalizeData(seurat_obj) %>%  FindVariableFeatures()
seurat_obj <- RunFastMNN(object.list = SplitObject(seurat_obj, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

saveRDS(seurat_obj,'Rds/biobank_multiome_scRNA_seurat_obj_umap.Rds')
#seurat_obj <- readRDS('Rds/biobank_multiome_scRNA_seurat_obj_umap.Rds')
p1 <- show_umap_plot(seurat_obj,raster=FALSE)

seurat_obj <- seurat_obj %>% FindClusters(resolution = 0.3)

p2 <- show_umap_plot(seurat_obj,label = FALSE,raster=FALSE)
saveplot(p2,filenames = glue('{sample}_raw_umap'),width = 4,height = 4)


islet_makers <- readRDS('../sample12/islet_makers.Rds')
islet_makers_filter <- 

makers <- c(
  'GCG','INS','SST','PPY', # alpha beta delta gamma
  'CFTR','KRT19','REG1A','PRSS1', # Ductal Acinar
  'PECAM1','COL1A1','PTPRC' #Endo, Stellate, Immune
)


p3 <- my_DotPlot(seurat_obj,genes = makers)
saveplot(p3,filenames = glue('{sample}_dotplot'),width = 6,height = 3)

p4 <- FeaturePlot2(seurat_obj,gene = makers,ncol = 3,raster = FALSE)
saveplot(p4,filenames = glue('{sample}_feature_plot'),width = 6,height = 8)

###beta cell-----
beta_cell <- subset(seurat_obj,idents = c(0,3))
beta_cell <- NormalizeData(beta_cell) %>%  FindVariableFeatures()
beta_cell <- RunFastMNN(object.list = SplitObject(beta_cell, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

beta_cell <- beta_cell %>% FindClusters(resolution = 0.3)

p1 <- show_umap_plot(beta_cell,raster=FALSE)

beta_cell$cell_type <- 'Beta'
Idents(beta_cell) <- beta_cell$cell_type 
saveRDS(beta_cell,file = 'Rds/cell_type/beta_cell.Rds')

###Alpha cell-----
alpha_cell <- subset(seurat_obj,idents = c(1,2,5,10))
alpha_cell <- NormalizeData(alpha_cell) %>%  FindVariableFeatures()
alpha_cell <- RunFastMNN(object.list = SplitObject(alpha_cell, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

show_umap_plot(alpha_cell,label = F,raster=FALSE)
alpha_cell <- alpha_cell %>% FindClusters(resolution = 0.3)
p1 <- show_umap_plot(alpha_cell,raster=FALSE)

FeaturePlot(alpha_cell,features = 'GCG')
my_DotPlot(alpha_cell,genes = makers)

alpha_cell <- subset(alpha_cell,idents = c(5),invert = T)

alpha_cell$cell_type <- 'Alpha'
Idents(alpha_cell) <- alpha_cell$cell_type 
saveRDS(alpha_cell,file = 'Rds/cell_type/alpha_cell.Rds')

###delta cell-----
delta_cell <- subset(seurat_obj,idents = c(6,13,17))
delta_cell <- NormalizeData(delta_cell) %>%  FindVariableFeatures()
delta_cell <- RunFastMNN(object.list = SplitObject(delta_cell, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

show_umap_plot(delta_cell,label = F,raster=FALSE)
delta_cell <- delta_cell %>% FindClusters(resolution = 0.3)
FeaturePlot(delta_cell,features = 'PPY')
my_DotPlot(delta_cell,genes = makers)

delta_cell <- subset(delta_cell,idents = c(0,1,5))

delta_cell$cell_type <- 'Delta'
Idents(delta_cell) <- delta_cell$cell_type 
saveRDS(delta_cell,file = 'Rds/cell_type/delta_cell.Rds')


###gamma cell-----
gamma_cell <- subset(seurat_obj,idents = c(9))
gamma_cell <- NormalizeData(gamma_cell) %>%  FindVariableFeatures()
gamma_cell <- RunFastMNN(object.list = SplitObject(gamma_cell, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

show_umap_plot(gamma_cell,label = F,raster=FALSE)
gamma_cell <- gamma_cell %>% FindClusters(resolution = 0.3)
FeaturePlot(gamma_cell,features ='PPY')
my_DotPlot(gamma_cell,genes = makers)

gamma_cell <- subset(gamma_cell,idents = c(0,1))

gamma_cell$cell_type <- 'Gamma'
Idents(gamma_cell) <- gamma_cell$cell_type 
saveRDS(gamma_cell,file = 'Rds/cell_type/gamma_cell.Rds')

###ductal cell-----
ductal_cell <- subset(seurat_obj,idents = c(12))
ductal_cell <- NormalizeData(ductal_cell) %>%  FindVariableFeatures()
ductal_cell <- RunFastMNN(object.list = SplitObject(ductal_cell, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

show_umap_plot(ductal_cell,label = F,raster=FALSE)
ductal_cell <- ductal_cell %>% FindClusters(resolution = 0.3)
FeaturePlot(ductal_cell,features =  c('CFTR','KRT19'))
my_DotPlot(ductal_cell,genes = makers)

ductal_cell$cell_type <- 'Ductal'
Idents(ductal_cell) <- ductal_cell$cell_type 
saveRDS(ductal_cell,file = 'Rds/cell_type/ductal_cell.Rds')


###acinar cell-----
acinar_cell <- subset(seurat_obj,idents = c(4))
acinar_cell <- NormalizeData(acinar_cell) %>%  FindVariableFeatures()
acinar_cell <- RunFastMNN(object.list = SplitObject(acinar_cell, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

show_umap_plot(acinar_cell,label = F,raster=FALSE)
acinar_cell <- acinar_cell %>% FindClusters(resolution = 0.3)
FeaturePlot(acinar_cell,features =  c('REG1A','PRSS1'))
my_DotPlot(acinar_cell,genes = makers)

acinar_cell$cell_type <- 'Acinar'
Idents(acinar_cell) <- acinar_cell$cell_type 
saveRDS(acinar_cell,file = 'Rds/cell_type/acinar_cell.Rds')



##run for find other cell type-----
main_celltype_obj <- map(list.files('Rds/cell_type/',full.names = T),readRDS) %>% reduce(merge)


other_cells_seurat_obj <- subset(seurat_obj,cells = colnames(main_celltype_obj), invert = T) 

other_cells_seurat_obj <- NormalizeData(other_cells_seurat_obj) %>%  FindVariableFeatures()
other_cells_seurat_obj <- RunFastMNN(object.list = SplitObject(other_cells_seurat_obj, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

p11 <- DimPlot(other_cells_seurat_obj,label = T)
p12 <- my_DotPlot(other_cells_seurat_obj,genes = makers)

###Endo cell-----

endo_cell <- subset(other_cells_seurat_obj,idents = c(5))
endo_cell <- NormalizeData(endo_cell) %>%  FindVariableFeatures()
endo_cell <- RunFastMNN(object.list = SplitObject(endo_cell, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

endo_cell$cell_type <- 'Endothelial'
show_umap_plot(other_cells_seurat_obj,label = F,raster = F)
FeaturePlot(other_cells_seurat_obj,features = 'PECAM1')
FeaturePlot(other_cells_seurat_obj,features = 'CLEC14A')
my_DotPlot(other_cells_seurat_obj,genes = c('PECAM1','CLEC14A'))
show_umap_plot(endo_cell,group_by = 'cell_type',label = F,raster = F)

Idents(endo_cell) <- endo_cell$cell_type 
saveRDS(endo_cell,file = 'Rds/cell_type/endo_cell.Rds')

###Stellate cell-----
Stellate <- subset(other_cells_seurat_obj,idents = c(11))
Stellate <- NormalizeData(Stellate) %>%  FindVariableFeatures()
Stellate <- RunFastMNN(object.list = SplitObject(Stellate, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

Stellate <- Stellate %>% FindClusters(resolution = 0.3)

Stellate <- subset(Stellate,idents = 2)
Stellate$cell_type <- 'Stellate'
Idents(Stellate) <- Stellate$cell_type 
saveRDS(Stellate,file = 'Rds/cell_type/Stellate.Rds')

###immune cell

immune <- subset(other_cells_seurat_obj,idents = c(16))
immune <- NormalizeData(immune) %>%  FindVariableFeatures()
immune <- RunFastMNN(object.list = SplitObject(immune, split.by = "orig.ident")) %>%
  RunUMAP(reduction = "mnn", dims = 1:30) %>%
  FindNeighbors(reduction = "mnn", dims = 1:30)

immune$cell_type <- 'immune'
Idents(immune) <- immune$cell_type 
show_umap_plot(immune)

my_DotPlot(immune,genes = c('PTPRC','CCL3','CD68','CD74'))
saveRDS(immune,file = 'Rds/cell_type/immune.Rds')


##Run for integrate analysis------
all_cell_type <- list.files('Rds/cell_type/',full.names = T)

islet_intergrate <- map(all_cell_type,readRDS) %>% reduce(merge)
saveRDS(islet_intergrate,file = 'Rds/islet_intergrate.Rds')
library(SeuratWrappers)
library(harmony)
islet_intergrate_harmory <- NormalizeData(islet_intergrate) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
islet_intergrate_harmory <- islet_intergrate_harmory %>% RunHarmony(group.by.vars = "orig.ident")  
islet_intergrate_harmory <- islet_intergrate_harmory %>% 
  RunUMAP(min.dist = 0.3,umap.method = 'umap-learn',metric = 'correlation',reduction = "harmony", dims = 1:30)


show_umap_plot(islet_intergrate_harmory,group_by = 'cell_type',label = F,raster = F)
islet_intergrate_harmory <- islet_intergrate_harmory %>% FindNeighbors()
Idents(islet_intergrate_harmory) <- islet_intergrate_harmory$cell_type

islet_intergrate_harmory <- FindSubCluster(islet_intergrate_harmory,graph.name = 'RNA_snn',cluster = 'Alpha',resolution = 1)

DimPlot(islet_intergrate_harmory,group.by = 'sub.cluster',raster=FALSE,label = T) 
islet_intergrate_harmory <- subset(islet_intergrate_harmory,sub.cluster %in% c('Alpha_11','Alpha_17','Alpha_18'),invert = T)

islet_intergrate_harmory <- FindSubCluster(islet_intergrate_harmory,graph.name = 'RNA_snn',cluster = 'Beta',resolution = 1)
DimPlot(islet_intergrate_harmory,group.by = 'sub.cluster',raster=FALSE,label = T) 
islet_intergrate_harmory <- subset(islet_intergrate_harmory,sub.cluster %in% c('Beta_11','Beta_15'),invert = T)

islet_intergrate_harmory <- FindSubCluster(islet_intergrate_harmory,graph.name = 'RNA_snn',cluster = 'Delta',resolution = 2)
DimPlot(islet_intergrate_harmory,group.by = 'sub.cluster',raster=FALSE,label = T) 
islet_intergrate_harmory <- subset(islet_intergrate_harmory,sub.cluster %in% c('Delta_10'),invert = T)

p3 <- show_umap_plot(islet_intergrate_harmory,group_by = 'cell_type',label = F,raster = F)
saveplot(p3,filenames = glue('{sample}_umap_intergation'),width = 4,height = 4)

p4 <- show_umap_plot(islet_intergrate_harmory,group_by = 'orig.ident',raster = F)
saveplot(p4,filenames = glue('{sample}_umap_intergation_donor'),width = 5,height = 4)


p5 <- my_DotPlot(islet_intergrate_harmory,genes = makers)
saveplot(p5,filenames = glue('{sample}_dotplot'),width = 6,height = 3)


p6 <- FeaturePlot2(islet_intergrate_harmory,gene = makers,ncol = 3,raster = FALSE)
saveplot(p6,filenames = glue('{sample}_feature_plot'),width = 6,height = 8)


islet_intergrate_harmory <- islet_intergrate_harmory %>% 
  RunUMAP(min.dist = 0.3,reduction = "harmony", dims = 1:30)
show_umap_plot(islet_intergrate_harmory,group_by = 'cell_type',raster = F)
saveRDS(islet_intergrate_harmory,file = 'Rds/islet_intergrate_harmory.Rds')


islet_biobank_mutliome_metadata <- islet_intergrate_harmory@meta.data
saveRDS(islet_biobank_mutliome_metadata,'Rds/islet_biobank_mutliome_metadata.Rds')

##final----
peak_matrix <- readRDS('Rds/islet_reproducible_peak_matrix.Rds')
all(colnames(peak_matrix) %in% colnames(islet_intergrate_harmory))
islet_intergrate_harmory_final <- subset(islet_intergrate_harmory,cells = colnames(peak_matrix))
saveRDS(islet_intergrate_harmory_final,file = 'Rds/islet_intergrate_harmory_final.Rds')
islet_intergrate_harmory_final <- readRDS('Rds/islet_intergrate_harmory_final.Rds')


islet_biobank_mutliome_metadata_final <- islet_intergrate_harmory_final@meta.data
islet_biobank_mutliome_metadata_final <- islet_biobank_mutliome_metadata_final %>% rownames_to_column(var = 'barcode')
donor_info <- readxl::read_excel('Paper/donor_infomation.xlsx')
islet_biobank_mutliome_metadata_final <- left_join(islet_biobank_mutliome_metadata_final,donor_info,by = c('orig.ident'='Sample ID'))

islet_biobank_mutliome_metadata_final$status <- case_when(islet_biobank_mutliome_metadata_final$`Diabetes status` == 'Non-diabetic' ~ 'ND',
          .default = as.character(islet_biobank_mutliome_metadata_final$`Diabetes status`))


islet_biobank_mutliome_metadata_final$donor <- paste0('Biobank_multiome_',
                                                      islet_biobank_mutliome_metadata_final$Donor,
                                                      '_',islet_biobank_mutliome_metadata_final$status)
islet_biobank_mutliome_metadata_final <- as.data.frame(islet_biobank_mutliome_metadata_final)
rownames(islet_biobank_mutliome_metadata_final) <- islet_biobank_mutliome_metadata_final$barcode
identical(colnames(islet_intergrate_harmory_final),rownames(islet_biobank_mutliome_metadata_final))
islet_intergrate_harmory_final@meta.data <- islet_biobank_mutliome_metadata_final
saveRDS(islet_biobank_mutliome_metadata_final,'Rds/islet_biobank_mutliome_metadata_final.Rds')

saveRDS(islet_intergrate_harmory_final,file = 'Rds/islet_intergrate_harmory_final.Rds')
dior::seurat_write_h5(islet_intergrate_harmory_final,file = 'islet_intergrate_harmory_final.h5',assay.name = 'RNA',save.graphs = F)
#islet_intergrate_harmory_final <- readRDS('Rds/islet_intergrate_harmory_final.Rds')
DefaultAssay(islet_intergrate_harmory_final)
levels(islet_intergrate_harmory_final)
islet_makers <- FindAllMarkers(islet_intergrate_harmory_final,logfc.threshold = 1,min.pct = 0.5)
saveRDS(islet_makers,file = 'Rds/islet_makers.Rds')
islet_makers_filter <- islet_makers %>% group_by(cluster) %>% slice_max(n = 10,order_by = avg_log2FC)



p3 <- show_umap_plot(islet_intergrate_harmory_final,group_by = 'cell_type',label = F,raster = F)
saveplot(p3,filenames = glue('{sample}_umap_intergation_final'),width = 4,height = 4)

p4 <- show_umap_plot(islet_intergrate_harmory_final,group_by = 'orig.ident',raster = F)
saveplot(p4,filenames = glue('{sample}_umap_intergation_donor_final'),width = 5,height = 4)

p5 <- my_DotPlot(islet_intergrate_harmory_final,genes = makers)
saveplot(p5,filenames = glue('{sample}_dotplot_final'),width = 6,height = 3)

p6 <- FeaturePlot2(islet_intergrate_harmory_final,gene = makers,ncol = 3,raster = FALSE)
saveplot(p6,filenames = glue('{sample}_feature_plot_final'),width = 6,height = 8)



multiome_islet <- subset(islet_intergrate_harmory_final, idents = c('Alpha', 'Beta'))
multiome_islet <- NormalizeData(multiome_islet) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
multiome_islet <- multiome_islet %>% RunHarmony(group.by.vars = "orig.ident") %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% FindNeighbors()

Idents(multiome_islet) <- multiome_islet$cell_type

p1 <- VlnPlot_qc(multiome_islet,group.by = 'sample',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)
saveplot(p1,filenames = glue('{sample}_final_qc'),width = 6,height = 4)

p2 <- VlnPlot_qc(multiome_islet,group.by = 'donor',features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)
saveplot(p2,filenames = glue('{sample}_final_qc_donor'),width = 16,height = 4)


p7 <- show_umap_plot(multiome_islet,label = F)
saveplot(p7,filenames = glue('{sample}_alpha_beta_umap_intergation'),width = 4,height = 4)

p8 <- show_umap_plot(multiome_islet,group_by = 'donor',label = T)
saveplot(p8,filenames = glue('{sample}_alpha_beta_umap_intergation_donor'),width = 10,height = 4)

p6 <- my_DotPlot(multiome_islet,group_by = 'cell_type',genes = all_makers,scale = F)
saveplot(p6,filenames = glue('{sample}_dotplot'),width = 6,height = 3)

p4 <- FeaturePlot2(multiome_islet,gene = islet_top_makers,ncol = 3)
saveplot(p4,filenames = glue('{sample}_feature_plot'),width = 6,height = 6)
table(multiome_islet$donor,multiome_islet$cell_type)

saveRDS(multiome_islet,glue('Rds/{sample}_alpha_beta_final.Rds'),compress = FALSE)
multiome_islet <- readRDS(glue('Rds/{sample}_alpha_beta_final.Rds'))



multiome_beta_pseudobulk <- generate_pseudobulk(multiome_islet,multiome_islet$donor,mincell = 50)
#multiome_beta_pseudobulk <- multiome_beta_pseudobulk[,-str_which(colnames(multiome_beta_pseudobulk),'Pre')]
#remove Pre-T2D
write.csv(multiome_beta_pseudobulk,file = glue('{sample}_beta_pseudobulk.csv'))
saveRDS(multiome_beta_pseudobulk,glue('Rds/{sample}_beta_pseudobulk.Rds'))
multiome_beta_pseudobulk <- readRDS('Rds/biobank_multiome_beta_pseudobulk.Rds')
multiome_beta_pseudobulk_tpm <- get_tpm(multiome_beta_pseudobulk)
multiome_beta_pseudobulk_tpm <- multiome_beta_pseudobulk_tpm %>% as.data.frame() %>% 
    mutate(TPM_mean = rowMeans(.),gene = rownames(multiome_beta_pseudobulk_tpm)) %>%
    select(gene,TPM_mean)
data.table::fwrite(multiome_beta_pseudobulk_tpm,file = 'multiome_beta_pseudobulk_tpm.txt',col.names = F,sep = '\t')

multiome_beta_pseudobulk_tpm <- multiome_beta_pseudobulk_tpm %>% filter(TPM_mean >= 1)
data.table::fwrite(multiome_beta_pseudobulk_tpm[1],file = 'multiome_beta_UbiquitouslyExpressedGenes.txt',col.names = F,sep = '\t')


multiome_Alpha_pseudobulk <- generate_pseudobulk(multiome_islet,multiome_islet$donor,celltype = 'Alpha',mincell = 50)
#multiome_Alpha_pseudobulk <- multiome_Alpha_pseudobulk[,-str_which(colnames(multiome_Alpha_pseudobulk),'Pre')]
write.csv(multiome_Alpha_pseudobulk,file = glue('{sample}_Alpha_pseudobulk.csv'))
saveRDS(multiome_Alpha_pseudobulk,glue('Rds/{sample}_Alpha_pseudobulk.Rds'))
multiome_Alpha_pseudobulk <- readRDS('Rds/biobank_multiome_Alpha_pseudobulk.Rds')
multiome_Alpha_pseudobulk_tpm <- get_tpm(multiome_Alpha_pseudobulk)
multiome_Alpha_pseudobulk_tpm <- multiome_Alpha_pseudobulk_tpm  %>% as.data.frame() %>% 
  mutate(TPM_mean = rowMeans(.),gene = rownames(multiome_Alpha_pseudobulk_tpm)) %>%
  select(gene,TPM_mean)
data.table::fwrite(multiome_Alpha_pseudobulk_tpm,file = 'multiome_alpha_pseudobulk_tpm.txt',col.names = F,sep = '\t')

multiome_alpha_pseudobulk_tpm <- multiome_Alpha_pseudobulk_tpm %>% filter(TPM_mean >= 1)
data.table::fwrite(multiome_Alpha_pseudobulk_tpm[1],file = 'multiome_alpha_UbiquitouslyExpressedGenes.txt',col.names = F,sep = '\t')








##hg38 to hg19
blacllsit <- rtracklayer::import.bed('../../Encode blacklist/hg19-blacklist.v2.bed.gz')
library(GenomicRanges)
hg19_peaks <- readxl::read_excel('Paper/Peaks_hg19.xlsx')

hg19_peaks <- hg19_peaks[,1:3]
hg19_peaks_gr <- makeGRangesFromDataFrame(hg19_peaks)
hg19_peaks_gr <- hg19_peaks_gr[hg19_peaks_gr@seqnames %in% levels(peak_matrix@rowRanges@seqnames)]
hg19_peaks_gr <- hg19_peaks_gr[-queryHits(findOverlaps(hg19_peaks_gr,blacllsit, type="any"))]
rtracklayer::export.bed(hg19_peaks_gr,'Paper/gaowei_peaks_hg19.bed')



hg38_peaks_gr <- rtracklayer::import.bed('Paper/Peaks_hg38.bed')
rtracklayer::export.bed(peak_matrix@rowRanges,'Paper/Biobank_multiomr_Peaks_hg38.bed')
biobank_hg19_peaks_gr <- rtracklayer::import.bed('Paper/Biobank_multiomr_Peaks_hg19.bed')
biobank_hg19_peaks_gr <- biobank_hg19_peaks_gr[biobank_hg19_peaks_gr@seqnames %in% levels(peak_matrix@rowRanges@seqnames)]
biobank_hg19_peaks_gr <- biobank_hg19_peaks_gr[-queryHits(findOverlaps(biobank_hg19_peaks_gr,blacllsit, type="any"))]



length(subsetByOverlaps(hg19_peaks_gr,biobank_hg19_peaks_gr)) / length(hg19_peaks_gr)
length(subsetByOverlaps(hg19_peaks_gr,biobank_hg19_peaks_gr)) / length(biobank_hg19_peaks_gr)

##Create multiome object--------
library(Signac)
islet_singac <- CreateSeuratObject(
  counts = islet_intergrate_harmory_final@assays$RNA@counts,
  assay = "RNA",
  meta.data = islet_intergrate_harmory_final@meta.data
)
sparse_mtx <- peak_matrix@assays@data$PeakMatrix
colnames(sparse_mtx) <- colnames(peak_matrix)
rownames(sparse_mtx) <- paste(as.data.frame(peak_matrix@rowRanges)[[1]],as.data.frame(peak_matrix@rowRanges)[[2]],as.data.frame(peak_matrix@rowRanges)[[3]],sep = '-')

islet_singac[["ATAC"]] <- CreateChromatinAssay(
  counts = sparse_mtx,
  sep = c("-", "-")
)
saveRDS(islet_singac,'Rds/islet_multiome_seurat_singac.Rds')
islet_singac <- readRDS('Rds/islet_multiome_seurat_singac.Rds')




###------
library(Signac)
library(Seurat)
islet_singac@assays$RNA@counts
#show RNA gene by cell count matrix
islet_singac@assays$ATAC@counts
#show ATAC peak by cell count matrix
islet_singac@meta.data
#show cell metadata

cell_names <- sample(colnames(islet_singac),2000)
#random select 2000 cells
islet_singac@assays$RNA@counts[,cell_names]
islet_singac@assays$ATAC@counts[,cell_names]
islet_singac@meta.data[cell_names,]



library(data.table)
tpm <- data.table::fread('multiome_RNA_TPM_4.csv') %>% as.data.frame()
rownames(tpm) <- tpm$X
tpm <- tpm[,-1]

donor_info <- readxl::read_excel('Paper/donor_infomation.xlsx')

df <- t(tpm['SLC30A8',]) %>% as.data.frame() %>% rownames_to_column('sample') %>%
    mutate(ID = str_extract(sample,'[a-zA-Z0-9]+')) %>% 
    left_join(donor_info,by = 'ID') %>% filter(Donor != 'MM89')
df$SLC30A8 <- log1p(df$SLC30A8)
df$status <- df$`Diabetes status`

df$beta <- if_else(str_detect(df$sample,'beta1'),'beta1','beta2')

library(ggpubr)

data2 <- df %>% group_by(status) %>% summarise(mean = mean(SLC30A8))

pdf(file = 'SLC30A8_status.pdf',width = 3,height = 3)
ggplot(df,aes(status,SLC30A8,group=status)) +
    geom_boxplot() +
    geom_point() +
    geom_line(data = data2,aes(factor(status), y = mean),group  = 'status') +
    theme_bw()
dev.off()




p1 <- ggscatter(
    df, x = 'SLC30A8', y = 'HbA1c', add = "reg.line", color = 'status',
    palette = clustcol[1:length(unique(df[['status']]))],
    shape = 'beta',
    add.params = list(color = "blue", fill = "lightgray"),
    conf.int = TRUE)  +
    geom_text_repel(aes(label = sample)) +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
             method = 'spearman') +
    theme_bw() +
    theme(legend.title = element_blank())


pdf(file = 'SLC30A8_HbA1c.pdf',width = 12,height = 10)
p1
dev.off()


df_beta1 <- df %>% filter(str_detect(sample,'beta1'))
p2 <- ggscatter(
    df_beta1, x = 'SLC30A8', y = 'HbA1c', add = "reg.line", color = 'status',
    palette = clustcol[1:length(unique(df[['status']]))],
    add.params = list(color = "blue", fill = "lightgray"),
    conf.int = TRUE)  +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
             method = 'spearman') +
    theme_bw() +
    theme(legend.title = element_blank())
pdf(file = 'SLC30A8_HbA1c_beta1.pdf',width = 5,height = 4)
p2
dev.off()

df_beta2 <- df %>% filter(str_detect(sample,'beta2'))
p3 <- ggscatter(
    df_beta2, x = 'SLC30A8', y = 'HbA1c', add = "reg.line", color = 'status',
    palette = clustcol[1:length(unique(df[['status']]))],
    add.params = list(color = "blue", fill = "lightgray"),
    conf.int = TRUE)  +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
             method = 'spearman') +
    theme_bw() +
    theme(legend.title = element_blank())
pdf(file = 'SLC30A8_HbA1c_beta2.pdf',width = 5,height = 4)
p3
dev.off()


snATAC_beta_subtype <- data.table::fread('gaowei_data/biobank_snATAC_beta_subtype.csv') %>% as.data.frame()
rownames(snATAC_beta_subtype) <- snATAC_beta_subtype$X

snATAC_beta_subtype_cpm <- edgeR::cpm(snATAC_beta_subtype[,-1])
snATAC_beta_subtype_cpm <- snATAC_beta_subtype_cpm %>% as.data.frame() %>% rownames_to_column('peak')



peak2gene_10kb <- data.table::fread('gaowei_data/peak2gene_10kb.csv') %>% filter(gene =='SLC30A8')


snATAC_beta_SLC30A8 <- snATAC_beta_subtype_cpm %>% filter(peak %in% peak2gene_10kb$peak) %>% 
    pivot_longer(cols = -peak,names_to = 'sample',values_to = 'SLC30A8') %>% mutate(ID = str_extract(sample,'[a-zA-Z0-9]+'))  %>%
    left_join(donor_info,by = c('ID'='Donor')) %>% 
    filter(ID != 'MM89')
snATAC_beta_SLC30A8$SLC30A8 <- log1p(snATAC_beta_SLC30A8$SLC30A8)
snATAC_beta_SLC30A8$status <- snATAC_beta_SLC30A8$`Diabetes status`
snATAC_beta_SLC30A8$beta <- if_else(str_detect(snATAC_beta_SLC30A8$sample,'beta1'),'beta1','beta2')


snATAC_beta_SLC30A8_list <- snATAC_beta_SLC30A8 %>% group_split(peak)


plot_cor <- function(df){
    p1 <- ggscatter(
        df, x = 'SLC30A8', y = 'HbA1c', add = "reg.line", color = 'status',
        palette = clustcol[1:length(unique(df[['status']]))],
        shape = 'beta',
        title = unique(df[[1]]),
        add.params = list(color = "blue", fill = "lightgray"),
        conf.int = TRUE)  +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                 method = 'spearman') +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_blank())
    return(p1)
}
atac_cor_list <- map(snATAC_beta_SLC30A8_list,plot_cor)
pdf(file = 'SLC30A8_HbA1c_atac.pdf',width = 16,height = 6)
cowplot::plot_grid(plotlist = atac_cor_list,ncol = 4)
dev.off()

