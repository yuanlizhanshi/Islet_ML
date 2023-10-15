##basic function---
saveplot <- function(object,filenames = filename,width = NULL,height = NULL,dpi = 600){
  if (is.null(width)| is.null(height)) {
    stop('Must input figure width and height')
  }
  
  if(!dir.exists("Fig")){
    dir.create("Fig")
  }
  map(c("TIFF","PNG","PDF"),function(x){
    if(!dir.exists(paste0('Fig/',x))){
      dir.create(paste0('Fig/',x))
    }
  })
  
  pdf_filenames = paste0('Fig/PDF/',filenames,'.pdf')
  pdf_width = width*2
  pdf_height = height*2
  pdf_res = dpi
  
  tiff_filenames = paste0('Fig/TIFF/',filenames,'.tiff')
  tiff_width = width*1000
  tiff_height = height*1000
  tiff_res = dpi
  
  png_filenames = paste0('Fig/PNG/',filenames,'.png')
  png_width = width*500
  png_height = height*500
  if (dpi > 600) {
    png_res = 300
  }else{
    png_res = dpi/2
  }
  
  
  cairo_pdf(file = pdf_filenames,width = pdf_width, height = pdf_height,fallback_resolution = pdf_res)
  print(object)
  dev.off()
  
  tiff(filename = tiff_filenames,width = tiff_width, height = tiff_height, units = "px", res = tiff_res, compression = "lzw")
  print(object)
  dev.off()
  
  png(filename = png_filenames,width = png_width, height = png_height, res = png_res,units = "px")
  print(object)
  dev.off()
  
}
generate_pseudobulk <- function(sce,group_by = NULL,celltype = 'Beta',mincell = 50,ncore = NULL){
  if (is.null(group_by)) {
    stop('Must supply donor info')
  }
  beta <- subset(sce,idents = celltype)
  mat <- beta@assays$RNA@counts
  
  cell_type <-  names(which(table(beta$donor) > mincell))
  if (length(names(which(table(beta$donor) < mincell)))) {
    cat(names(which(table(beta$donor) < mincell)),'discarded due to low numbers of cells\n')
  }
  
  row_name <- rownames(mat)
  if (is.null(ncore)) {
    result <- map_dfc(cell_type,function(x,mat){
      col_index <- which(beta$donor == x)  
      mat <- mat[,col_index]
      mat_rowsum <- as.data.frame(Matrix::rowSums(mat)) 
      colnames(mat_rowsum) <- x
      cat(x,'pseudobulk generated with',length(col_index),' cells\n')
      return(mat_rowsum)
    },mat = mat)
  }else{
    require(furrr)
    plan(multisession, workers = ncore)
    result <- furrr::future_map_dfc(cell_type,function(x,mat){
      col_index <- which(colnames(mat) == x)  
      mat <- mat[,col_index]
      mat_rowsum <- as.data.frame(Matrix::rowSums(mat)) 
      colnames(mat_rowsum) <- x
      return(mat_rowsum)
    },mat = mat)
  }
  return(result)
}


##Seurat wrapper function----
run_seurat <- function(sce,cluster = FALSE, resolution = 0.3){
  if (cluster){
    sce <- SCTransform(sce) %>% RunPCA(verbose = FALSE) %>% 
      RunUMAP(dims = 1:30,umap.method = 'umap-learn',metric = 'correlation',verbose = FALSE)   %>%
      FindNeighbors() %>% 
      FindClusters(resolution = resolution)
  }else{
    sce <- SCTransform(sce) %>% RunPCA(verbose = FALSE) %>% 
      RunUMAP(dims = 1:30,umap.method = 'umap-learn',metric = 'correlation',verbose = FALSE) 
  }
  
  DefaultAssay(sce) <- 'RNA'
  sce <-  sce %>% NormalizeData(normalization.method = 'RC',scale.factor = 20000) %>% ScaleData(scale.max = 20)
  return(sce)
}

