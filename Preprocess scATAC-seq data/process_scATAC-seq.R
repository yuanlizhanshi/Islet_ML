library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRThreads(threads = 60) 
addArchRGenome("hg38")

all_samples <- str_extract(list.files('/data/gaowei/biobank_fastq_multiome/biobank_multiome_out/library/'),'.*(?=\\.)')

fragments <- paste0('/data/gaowei/biobank_fastq_multiome/biobank_multiome_out/',all_samples,'/outs/atac_fragments.tsv.gz')
names(fragments) <- all_samples
map(seq_along(fragments),function(i){
    file.copy(fragments[i],paste0('fragment_file/',names(fragments)[i],'_fragments.tsv.gz'),overwrite =T)
})
inputFiles <- list.files('fragment_file',pattern = 'gz$',full.names = T)

##function------

#check best parameters for LSI 
Optimization <- function(iterations,sampleCells,varFeatures,archr_obj = NULL,groub_by = 'Sample',useMatrix = "TileMatrix"){
    archr_obj <- addIterativeLSI(
        ArchRProj = archr_obj,
        useMatrix = useMatrix, 
        name = "IterativeLSI", 
        iterations = iterations, 
        clusterParams = list( 
            resolution = c(0.1, 0.8,1.5, 2)[1:(iterations-1)], 
            sampleCells = sampleCells, 
            n.start = 10
        ), 
        sampleCellsPre = sampleCells,
        varFeatures = varFeatures, 
        saveIterations = F,
        dimsToUse = 1:30,
        force = TRUE 
    )
    
    archr_obj <- addHarmony(
        ArchRProj = archr_obj,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "Sample",
        corCutOff = 0.8,
        force = TRUE
    )
    
    archr_obj <- addUMAP(
        ArchRProj = archr_obj, 
        reducedDims = "Harmony", 
        name = "UMAPHarmony", 
        nNeighbors = 30, 
        minDist = 0.3, 
        metric = "cosine",
        force = TRUE
    )
    p1 <- plotEmbedding(ArchRProj = archr_obj, colorBy = "cellColData", name = groub_by, embedding = "UMAPHarmony")
    saveplot(p1, filenames  = paste0('iterations_',iterations,'_Cell_',sampleCells,'_Var_',varFeatures),height = 4,width = 4)
    return(archr_obj)
}

subset_cells <- function(
    ArchRProj,embedding = 'UMAPHarmony',plot_object,features =NULL,thershold = NULL
){
    embedding_df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
    if (features %in% names(plot_object)) {
        plot_data <- plot_object[[features]][["layers"]][[1]][["data"]] 
    }else{
        plot_data <- plot_object[["layers"]][[1]][["data"]]  
    }
    rownames(plot_data) <- rownames(embedding_df)
    return(rownames(plot_data)[which(plot_data$color >= thershold)])
}


subdivideGRanges <- function (x, subsize = 500,core = NULL) 
{
    subdivide.range <- function (start.pos, end.pos, subsize = 100) 
    {
        width <- end.pos - start.pos + 1
        if (width < 2 * subsize) {
            stop("Width is less than 2 times subsize")
        }
        nchunks <- floor(width/subsize)
        relative.starts <- round(0:(nchunks - 1) * width/nchunks)
        relative.ends <- round(1:nchunks * width/nchunks) - 1
        return(list(start.pos = start.pos + relative.starts, end.pos = start.pos + 
                        relative.ends, length = nchunks))
    }
    
    subdivideIRanges <- function (x, subsize = 100) 
    {
        if (length(x) == 0) {
            return(x)
        }
        start.pos <- start(x)
        end.pos <- end(x)
        widths <- width(x)
        nsubranges <- pmax(1, floor(widths/subsize))
        out.start.pos <- numeric(sum(nsubranges))
        out.end.pos <- numeric(sum(nsubranges))
        out.idx <- 1
        for (i in 1:length(x)) {
            if (widths[i] < 2 * subsize) {
                out.start.pos[out.idx] <- start.pos[i]
                out.end.pos[out.idx] <- end.pos[i]
                out.idx <- out.idx + 1
            }
            else {
                sr <- subdivide.range(start.pos[i], end.pos[i], 
                                      subsize)
                out.start.pos[out.idx:(out.idx + sr$length - 1)] <- sr$start.pos
                out.end.pos[out.idx:(out.idx + sr$length - 1)] <- sr$end.pos
                out.idx <- out.idx + sr$length
            }
        }
        IRanges(start = out.start.pos, end = out.end.pos)
    }
    
    
    if (length(x) == 0) {
        return(x)
    }
    if (length(subsize) > 1) {
        stop("The subsize argument should be a single number: the desired width of the subdivided ranges")
    }
    if (is.null(core)) {
        x <- map(seq_along(x),function(i){
            fold = floor(width(x[i])/500) 
            if (fold < 1 ) {
                fold <- 1
            }
            resize(x[i],width = 500 * fold ,fix = 'center')   
        },.progress = T) 
        x <- do.call(c,x)
    }else{
        require(furrr)
        plan(multisession,workers = core)
        options(future.globals.maxSize= 100*1024*1024^2) 
        x <- furrr::future_map(seq_along(x),function(i){
            fold = floor(width(x[i])/500) 
            if (fold < 1 ) {
                fold <- 1
            }
            resize(x[i],width = 500 * fold ,fix = 'center')   
        },.progress = T) 
        x <- do.call(c,x)
    } 
    x <- sort(GenomicRanges::reduce(x))
    gr_list <- lapply(levels(seqnames(x)), function(seqlvl) {
        if (!any(seqnames(x) == seqlvl)) {
            return(GRanges())
        }
        rg <- ranges(x[seqnames(x) == seqlvl])
        GRanges(seqnames = seqlvl, ranges = subdivideIRanges(rg,subsize), seqlengths = seqlengths(x))
    })
    do.call(c, gr_list) 
    
}



##Scan fragment------
ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = str_extract(inputFiles,'(?<=/).*(?=_f)'),
    minTSS = 3, 
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)

##remove doublets-----
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, 
    knnMethod = "UMAP",
    LSIMethod = 1
)


##Create ArchR object-----
biobank_multiome_atac <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = "biobank_multiome_atac",
    copyArrows = TRUE)
biobank_multiome_atac <- filterDoublets(biobank_multiome_atac)
saveArchRProject(ArchRProj = biobank_multiome_atac, outputDirectory = "biobank_multiome_atac", load = TRUE)


##load scRNA-seq object-----
rna_h5 <- paste0('/data/gaowei/biobank_fastq_multiome/biobank_multiome_out/',all_samples,'/outs/filtered_feature_bc_matrix.h5')

biobank_multiome_rna_h5 <- import10xFeatureMatrix(
    input = rna_h5,
    names = all_samples
)

#biobank_multiome_atac <- loadArchRProject('biobank_multiome_atac/')
biobank_multiome_rna <- readRDS('../scRNA-seq/islet_intergrate_harmory.Rds')
RNA_ATAC_intersect_cells <- intersect(colnames(biobank_multiome_rna),getCellNames(biobank_multiome_atac))
biobank_multiome_rna <- subset(biobank_multiome_rna,cells = RNA_ATAC_intersect_cells)
biobank_multiome_atac <- biobank_multiome_atac[RNA_ATAC_intersect_cells,]

identical(colnames(biobank_multiome_rna),getCellNames(biobank_multiome_atac))

biobank_multiome_atac <- addCellColData(biobank_multiome_atac,
                                        data = biobank_multiome_rna$cell_type,
                                        cells = colnames(biobank_multiome_rna),
                                        name = 'cell_type',
                                        force = T)

biobank_multiome_rna_h5 <- biobank_multiome_rna_h5[,RNA_ATAC_intersect_cells]

biobank_multiome_atac <- addGeneExpressionMatrix(input = biobank_multiome_atac, 
                                                 seRNA = biobank_multiome_rna_h5, 
                                                 force = TRUE)

##Using scATAC-seq for dimension reduction------

biobank_multiome_atac <- Optimization(iterations = 4,sampleCells = 10000,varFeatures = 250000,archr_obj = biobank_multiome_atac)
p1 <- plotEmbedding(ArchRProj = biobank_multiome_atac, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony")
plotPDF(p1, name = "Plot-UMAP2Harmony-cell_type_250k_features.pdf", ArchRProj = biobank_multiome_atac, addDOC = FALSE, width = 5, height = 5)


##Using scRNA-seq for dimension reduction------

biobank_multiome_atac <- addIterativeLSI(
    ArchRProj = biobank_multiome_atac, 
    clusterParams = list(
        resolution = 0.3, 
        sampleCells = 30000,
        n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI",
    sampleCellsPre = 60000,
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
)
biobank_multiome_atac <- addHarmony(
    ArchRProj = biobank_multiome_atac,
    reducedDims = "LSI_RNA",
    name = "Harmony",
    groupBy = "Sample",
    corCutOff = 0.8,
    force = TRUE
)
biobank_multiome_atac <- addUMAP(
    ArchRProj = biobank_multiome_atac, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony_RNA", 
    nNeighbors = 30, 
    minDist = 0.3, 
    metric = "cosine",
    force = TRUE
)
p1 <- plotEmbedding(ArchRProj = biobank_multiome_atac, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony_RNA")
saveplot(p1,filenames = 'bPlot-UMAP2Harmony-cell_type_RNA',width = 5,height = 5)
plotPDF(p1, name = "Plot-UMAP2Harmony-cell_type_RNA.pdf", ArchRProj = biobank_multiome_atac, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = biobank_multiome_atac, outputDirectory = "biobank_multiome_atac", load = TRUE)



##Process the left 14 scATAC-seq-----
###Scan fragment file-----
biobank_atac_sample <- list.files('/data/kyh/biobank_ATAC_sra/',pattern = 'SAMP')

biobank_atac_fragments <- paste0('/data/kyh/biobank_ATAC_sra/',biobank_atac_sample,'/outs/fragments.tsv.gz')
names(biobank_atac_fragments) <- biobank_atac_sample

map(seq_along(biobank_atac_fragments),function(i){
    file.copy(biobank_atac_fragments[i],paste0('fragment_file_scATAC/',names(biobank_atac_fragments)[i],'_fragments.tsv.gz'),overwrite =T)
})
##Create ArchR object------

ArrowFiles <- createArrowFiles(
    inputFiles = list.files('fragment_file_scATAC',full.names = T),
    sampleNames = str_extract(list.files('fragment_file_scATAC',full.names = T),'SAMP-\\d+'),
    minTSS = 3, 
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, 
    knnMethod = "UMAP",
    LSIMethod = 1
)

ArrowFiles_multiome <- list.files('biobank_multiome_atac/ArrowFiles',full.names = T)
biobank_atac <- ArchRProject(
    ArrowFiles = c(ArrowFiles), 
    outputDirectory = "biobank_atac_all",
    copyArrows = TRUE)
biobank_atac <- filterDoublets(biobank_atac)
saveArchRProject(biobank_atac,outputDirectory = "biobank_atac_all")
biobank_atac <- loadArchRProject(path = "biobank_atac_all",showLogo = F)


##subset delta gamma cell
delta_cell <-  getCellNames(biobank_atac)[which(biobank_atac$cell_type == 'delta')]
gamma_cell <-  getCellNames(biobank_atac)[which(biobank_atac$cell_type == 'gamma')]


##subset Ductal Acinar cell----
Ductal_Acinar <- biobank_atac[biobank_atac$cell_type %in% c('Ductal','Acinar'),]
Ductal_Acinar <- Optimization(iterations = 5,sampleCells = 3000,varFeatures = 50000,archr_obj = Ductal_Acinar,groub_by = 'cell_type')
Ductal_Acinar <- addClusters(
    input = Ductal_Acinar,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 3,
    maxClusters = 10,
    force = T
)
Ductal_Acinar <- addImputeWeights(Ductal_Acinar)
Ductal_Acinar_plot <- plotEmbedding(ArchRProj = Ductal_Acinar, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p <- plotEmbedding(
    ArchRProj = Ductal_Acinar, 
    colorBy = "GeneScoreMatrix", 
    name = c('CFTR','KRT19','REG1A','PRSS1'), 
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(Ductal_Acinar)
)
plotPDF(plotList = p, 
        name = "Ductal_Acinar_Imputation.pdf", 
        ArchRProj = biobank_atac, 
        addDOC = FALSE, width = 5, height = 5)
labelOld <- str_sort(names(table(Ductal_Acinar$Clusters)))
labelNew <- c(
    "C1" = "Ductal","C2" = "Ductal","C3" = "Ductal","C4" = "Ductal","C5" = "Acinar",
    "C6" = "Acinar","C7" = "Acinar","C8" = "Acinar","C9" = "Acinar","C10" = "Acinar"
)
Ductal_Acinar$cell_type <- mapvalues(Ductal_Acinar$Clusters, names(labelNew), labelNew)

Ductal_Acinar_plot <- plotEmbedding(ArchRProj = Ductal_Acinar, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony")
plotPDF(Ductal_Acinar_plot, 
        name = "Ductal_Acinar_plot.pdf", 
        ArchRProj = biobank_atac, 
        addDOC = FALSE, width = 5, height = 5)

Ductal_cell <-  getCellNames(Ductal_Acinar)[which(Ductal_Acinar$cell_type == 'Ductal')]
Acinar_cell <-  getCellNames(Ductal_Acinar)[which(Ductal_Acinar$cell_type == 'Acinar')]



### subset Endo, Stellate, Immune cell----

es_cell <- biobank_atac[biobank_atac$cell_type %in% c('Endo','Stellate','Immune'),]
es_cell <- Optimization(iterations = 5,sampleCells = 5000,varFeatures = 100000,archr_obj = es_cell,groub_by = 'cell_type')
es_cell <- addClusters(
    input = es_cell,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 3,
    maxClusters = 20,
    force = T
)
es_cell <- addImputeWeights(es_cell)
es_cell_plot <- plotEmbedding(ArchRProj = es_cell, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
#'PECAM1','COL1A1','ZEB2' #Endo, Stellate, Immune
p <- plotEmbedding(
    ArchRProj = es_cell, 
    colorBy = "GeneScoreMatrix", 
    name = c('PECAM1','COL1A1','ZEB2'), 
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(es_cell)
)
plotPDF(plotList = p, 
        name = "es_cell_Imputation.pdf", 
        ArchRProj = biobank_atac, 
        addDOC = FALSE, width = 5, height = 5)
labelOld <- str_sort(names(table(es_cell$Clusters)))
labelNew <- c(
    "C1" = "Endo","C2" = "Endo","C3" = "Endo","C4" = "Immune","C5" = "Immune",
    "C6" = "Immune","C7" = "Immune","C8" = "Endo","C9" = "Endo","C10" = "Stellate",
    "C11" = "Stellate","C12" = "Stellate","C13" = "Stellate","C14" = "Stellate","C15" = "Stellate",
    "C16" = "Stellate","C17" = "Stellate","C18" = "Stellate","C19" = "Stellate","C20" = "Stellate")
es_cell$cell_type <- mapvalues(es_cell$Clusters, names(labelNew), labelNew)

es_cell_plot <- plotEmbedding(ArchRProj = es_cell, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony")
plotPDF(es_cell_plot, 
        name = "esi_plot.pdf", 
        ArchRProj = biobank_atac, 
        addDOC = FALSE, width = 5, height = 5)

Endo_cell <-  getCellNames(es_cell)[which(es_cell$cell_type == 'Endo')]
Stellate_cell <-  getCellNames(es_cell)[which(es_cell$cell_type == 'Stellate')]
Immune_cell <-  getCellNames(es_cell)[which(es_cell$cell_type == 'Immune')]

celltype2 <- list(
    Alpha = alpha_cell,beta = beta_cell,
    delta = delta_cell,gamma = gamma_cell,
    Ductal = Ductal_cell,Acinar = Acinar_cell,
    Endo = Endo_cell,Stellate = Stellate_cell,Immune = Immune_cell
)

saveRDS(celltype2,'Rds/biobank_atac_celltype.Rds')
celltype2 <- readRDS('Rds/biobank_atac_celltype.Rds')
names(celltype2) <- c('Alpha','Beta','Delta','Gamma','Ductal','Acinar','Endothelial','Stellate','immune')


celltype_df <- map_dfr(1:length(celltype2),function(i){
    df <- tibble(
        barcode = celltype2[[i]],
        cell_type = names(celltype2)[i]
    )
    return(df)
})

all(celltype_df$barcode %in% getCellNames(biobank_atac))
    
biobank_atac <- biobank_atac[celltype_df$barcode,]

all(celltype_df$barcode %in% getCellNames(biobank_atac))
all(getCellNames(biobank_atac) %in% celltype_df$barcode)


biobank_atac <- addCellColData(biobank_atac,data = celltype_df$cell_type,cells = celltype_df$barcode,name = 'cell_type',force = T)
p4 <- plotEmbedding(ArchRProj = biobank_atac, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony")
saveArchRProject(biobank_atac,outputDirectory = "biobank_atac_all")
biobank_atac <- loadArchRProject('biobank_atac_all')
##

test_Peak <- list.files('signac_macs2_temp',pattern = 'bed',full.names = T)
test_peak_name <- stringr::str_extract(test_Peak,'(?<=p/)\\w+(?=_)')
singac_macs2_rep_peaks <- indentifyReproduciblePeaks(summitFiles = test_Peak,summitNames = test_peak_name,reproducibility = 0)
singac_macs2_rep_peaks <- singac_macs2_rep_peaks[singac_macs2_rep_peaks@seqnames %in% paste0('chr',c(1:22,'X')),]


biobank_atac <- loadArchRProject('biobank_atac_all')
celltype1 <-  biobank_atac$cell_type
names(celltype1) <- biobank_atac$cellNames

biobank_multiome <- loadArchRProject('biobank_multiome_atac')
celltype2 <-  biobank_multiome$cell_type
names(celltype2) <- biobank_multiome$cellNames

celltype <- c(celltype1,celltype2)

##Combine all cells------
arrow_files <- c(list.files('biobank_multiome_atac/ArrowFiles/',full.names = T),
                 list.files('biobank_atac_all/ArrowFiles/',full.names = T))
atac_all <- ArchRProject(
    ArrowFiles = arrow_files, 
    outputDirectory = "Atac_all",
    copyArrows = TRUE)
atac_all <- filterDoublets(atac_all)
atac_all <- addCellColData(atac_all,data = celltype,cells = names(celltype),name = 'cell_type',force = T)
atac_all <- Optimization(iterations = 5,sampleCells = 10000,varFeatures = 40 * 10000,archr_obj = atac_all,groub_by = 'cell_type') 
saveArchRProject(atac_all,'Atac_all',load = T)
atac_all <- loadArchRProject('Atac_all')
atac_all <- atac_all[atac_all]

p1 <- plotEmbedding(ArchRProj = atac_all, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony")

atac_all <- atac_all[getCellNames(atac_all)[which(!is.na(atac_all$cell_type))],]


##-------
maker_gene <- plotBrowserTrack(
    ArchRProj = biobank_atac, 
    groupBy = "cell_type", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)
map(
    seq_along(maker_gene),
    function(i){
        filename = paste0('biobank_atac_',names(maker_gene)[i])
        saveplot(maker_gene[[i]],filenames = filename,width = 4,height = 6)
    }
)
##
p2 <- plotGroups(
    ArchRProj = atac_all, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
saveplot(p2,filenames = 'biobank_atac_TSS_enrichment',width = 6,height = 3)

p3 <- plotGroups(
    ArchRProj = atac_all, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
saveplot(p3,filenames = 'biobank_atac_log10_Frags',width = 6,height = 3)

p4 <- plotGroups(
    ArchRProj = atac_all, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "PromoterRatio",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
saveplot(p4,filenames = 'biobank_atac_PromoterRatio',width = 6,height = 3)

p5 <- plotFragmentSizes(ArchRProj = atac_all)
p6 <- plotTSSEnrichment(ArchRProj = atac_all)

legend <- cowplot::get_legend(
    p5 + theme(legend.box.margin = margin(0, 0, 0, 12),
               legend.position = "bottom",
               legend.title = element_blank())+
        guides(color = guide_legend(nrow = 2)) 
)

p7 <- cowplot::plot_grid(p5+theme(legend.position="none"),
                         p6+theme(legend.position="none"),
                         labels = c("A", "B"),
                         legend,
                         ncol = 2,
                         rel_heights = c(1, .2))
saveplot(p7,filenames = 'biobank_atac_fragment_size_TSS',width = 6,height = 3)
non_overlap_gr <- hg38_peaks_gr[-queryHits(findOverlaps(hg38_peaks_gr, singac_macs2_rep_peaks, type="any"))]
test_plot <- plotBrowserTrack(atac_all,
                               region = resize(non_overlap_gr[1],width = 10000,fix = 'center'),
                               upstream = 50000,
                               downstream = 50000,
                               groupBy = "cell_type")

saveplot(p4,filenames = 'biobank_atac_PromoterRatio',width = 6,height = 3)
##-----
markerGenes <- c(
    'GCG','INS','SST','PPY', # alpha beta delta gamma
    'CFTR','KRT19','REG1A','PRSS1', # Ductal Acinar
    'PECAM1','COL1A1','ZEB2' #Endo, Stellate, Immune
)
atac_all <- addImputeWeights(atac_all)


p$GCG
alpha_cell <- subset_cells(atac_all,plot_object = p, features = 'GCG',thershold = 0.4)
p$INS
beta_cell <- subset_cells(atac_all,plot_object = p, features = 'INS',thershold = 1.2)
p$SST
delta_cell <- subset_cells(atac_all,plot_object = p, features = 'SST',thershold = 0.8)
p$PPY
gamma_cell <- subset_cells(atac_all,plot_object = p, features = 'PPY',thershold = 0.6)
p$CFTR
p$KRT19
Ductal_cell <-
    union(
        subset_cells(
            atac_all,
            plot_object = p,
            features = 'KRT19',
            thershold = 0.5
        ),
        subset_cells(
            atac_all,
            plot_object = p,
            features = 'CFTR',
            thershold = 0.7
        )
    )
p$REG1A
p$PRSS1
Acinar_cell <-
    union(
        subset_cells(
            atac_all,
            plot_object = p,
            features = 'REG1A',
            thershold = 0.5
        ),
        subset_cells(
            atac_all,
            plot_object = p,
            features = 'PRSS1',
            thershold = 0.9
        )
    )

p$PECAM1
Endo_cell <- subset_cells(atac_all,plot_object = p, features = 'PECAM1',thershold = 0.6)

p$COL1A1
Stellate_cell <- subset_cells(atac_all,plot_object = p, features = 'COL1A1',thershold = 0.7)

p$ZEB2
Immune_cells <- subset_cells(atac_all,plot_object = p, features = 'ZEB2',thershold = 0.8)
selected_cell <- purrr::reduce(list(alpha_cell,beta_cell,delta_cell,gamma_cell,
                                    Ductal_cell,Acinar_cell,
                                    Endo_cell,Stellate_cell,Immune_cells),union)
  
atac_all <- atac_all[selected_cell,]

##Alpha beta cell----
alpha_beta <- atac_all[getCellNames(atac_all)[which(atac_all$cell_type %in% c('Alpha', 'Beta'))],]
alpha_beta <- Optimization(iterations = 4,sampleCells = 40000,varFeatures = 10000,archr_obj = alpha_beta,groub_by = 'cell_type')
p1 <- plotEmbedding(ArchRProj = alpha_beta, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony")
alpha_beta <- addClusters(
    input = alpha_beta,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.2,
    maxClusters = 10,
    force = T
)
p2 <- plotEmbedding(ArchRProj = alpha_beta, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

plotEmbedding(
    ArchRProj = alpha_beta,
    colorBy = "cellColData",
    name = "cell_type",
    embedding = "UMAPHarmony",
    highlightCells = getCellNames(alpha_beta)[which(alpha_beta$Clusters != 'C1')]
)

alpha_beta <- alpha_beta[getCellNames(alpha_beta)[which(alpha_beta$Clusters != 'C1')],]
saveArchRProject(alpha_beta,outputDirectory = 'Alpha_beta')
p1 <- plotEmbedding(ArchRProj = alpha_beta, colorBy = "cellColData", name = "cell_type", embedding = "UMAPHarmony")
saveplot(p1,filenames = 'ATAC_all-UMAP2Harmony-alpha_beta',width = 4,height = 4)

alpha_beta <- addImputeWeights(alpha_beta)
p <- plotEmbedding(
    ArchRProj = alpha_beta, 
    colorBy = "GeneScoreMatrix", 
    name = c('GCG','INS'), 
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(alpha_beta)
)

gene_score_a_b <- p$GCG + p$INS
saveplot(gene_score_a_b,filenames = 'ATAC_all-UMAP2Harmony-alpha_beta_gene_score',width = 8,height = 4)
saveArchRProject(alpha_beta,outputDirectory = 'Alpha_beta')
##


signac_islet_all <- merge(seurat_atac_obj,seurat_multiome_atac_obj)
saveRDS(signac_islet_all,file = 'Rds/singac_raw_all_object.Rds')


##Call peaks with mcas2
#first load all rawdata to create singac object and assign the corrsponding cell type label, which assigned in ArchR workflow
#macs2 method in singac is more stable with previous study
signac_islet_all <- readRDS('Rds/singac_raw_all_object.Rds')
signac_islet_all <- subset(signac_islet_all,cells = colnames(signac_islet_all)[which(!is.na(signac_islet_all$cell_type))])
singac_macs2_peaks <- CallPeaks(signac_islet_all,
                          macs2.path = '/home/kyh/miniconda3/envs/deeptools/bin/macs2',
                          group.by = 'cell_type',
                          cleanup = F,
                          extsize = 150,
                          shift = -75,
                          combine.peaks = FALSE,
                          additional.args = '--keep-dup all',
                          outdir = 'signac_macs2_temp2')