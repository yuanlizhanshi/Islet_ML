##GRN construction------
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(motifmatchr)
library(TFBSTools)
library(JASPAR2020)
all_Peak <- list.files('signac_macs2_temp',pattern = 'bed',full.names = T)
all_peak_name <- stringr::str_extract(all_Peak,'(?<=p/)\\w+(?=_)')

###subdivideGRanges method-----
all_peak <- map(list.files('signac_macs2_temp/','narrow',full.names = T),import,format = 'narrowpeak') |> 
    purrr::reduce(c) |> GenomicRanges::reduce()
all_peak_gr <- subdivideGRanges(all_peak)
all_peak_gr$peak_name <- paste0('peak',1:length(all_peak_gr))
saveRDS(all_peak_gr,'signac_macs2_temp/all_peak_gr.Rds')
all_peak_gr <- readRDS('signac_macs2_temp/all_peak_gr.Rds')
all_peaks_df <- as.data.frame(all_peak_gr)

##This could be very slow....----
all_narrow_Peak <- list.files('signac_macs2_temp',pattern = '.narrow',full.names = T)    
all_narrow_Peak_gr <- map(all_narrow_Peak,rtracklayer::import,format = 'narrowPeak')  %>% 
    setNames(stringr::str_extract(all_narrow_Peak,'(?<=p/)\\w+(?=_)')) 
all_peaks_overlap <- map(1:length(all_narrow_Peak_gr),function(i){
    col <- overlapsAny(all_peak_gr,all_narrow_Peak_gr[[i]]) |> as.numeric()
    return(col)
},.progress = T) %>% setNames(names(all_narrow_Peak_gr)) %>% bind_cols()
#10/s
all_peaks_overlap <- all_peaks_overlap %>% mutate(peak = all_peak_gr$peak_name)
#259k
saveRDS(all_peaks_overlap,'Rds/all_peaks_overlap_temp.Rds')
##Start hear----
all_peaks_overlap <- readRDS('Rds/all_peaks_overlap_temp.Rds')
# all_peaks_overlap %>% select(-peak) %>% 
#     mutate(total = rowSums(.)) %>% select(total) %>% 
#     summarise(max = max(.),min = min(.))

beta_peaks <- all_peaks_overlap %>% filter(Beta == 1)
beta_peaks_gr <- all_peak_gr[all_peak_gr$peak_name %in% beta_peaks$peak]
# alpha_beta_peaks %>% mutate(total = Alpha + Beta) %>% 
#     select(total) %>% summarise(max = max(.),min = min(.))
#196k

#promoter
candidate_promoter <- subsetByOverlaps(
    beta_peaks_gr,
    ChIPseeker::getPromoters(
        TxDb.Hsapiens.UCSC.hg38.knownGene,
        upstream = 2000,
        downstream = 500,
        by = 'transcript'
    )
)
beta_peaks_gr$promoter <- as.numeric(overlapsAny(beta_peaks_gr,candidate_promoter)) 
beta_peaks_promoter_gr <- beta_peaks_gr[beta_peaks_gr$promoter == 1]

##Match motif----
#expressed_gene <- data.table::fread('/data/kyh/islet_bam/expression/multiome_beta_pseudobulk_tpm.txt') 
#     filter(V2 >= 1)
# motif_set_JASPAR2020 <- motif_set_JASPAR2020[which(motif_symbol_JASPAR2020 %in% expressed_gene$V1)]
motif_set_Transfac <- readRDS('Rds/motif_set/Transfac_PWMatrixList.rds')
motif_set_annotation <- readRDS('Rds/motif_set/Transfac_anno.Rds')
motif_set_Transfac <- motif_set_Transfac[motif_set_annotation$Accession]

#http://genexplain.com/transfac/

####Transfac
beta_promoter_promoter_Transfac <- matchMotifs(motif_set_Transfac, beta_peaks_promoter_gr, genome = "hg38",out = "positions") 

beta_promoter_promoter_Transfac <- map(1:length(beta_promoter_promoter_Transfac),function(i){
    beta_promoter_promoter_Transfac[[i]]$motif <- names(beta_promoter_promoter_Transfac)[i]
    return(beta_promoter_promoter_Transfac[[i]])
},.progress = T) %>% purrr::reduce(c)
##Accelerate (4h~)----
TF_set <- motif_set_annotation$TFs 
names(TF_set) <- motif_set_annotation$Accession
beta_promoter_promoter_Transfac$TF <- TF_set[beta_promoter_promoter_Transfac$motif]
saveRDS(beta_promoter_promoter_Transfac,'beta_promoter_promoter_Transfac_tmp.Rds')
##----
beta_TF_target_promoter_Transfac <- peak_anno(beta_promoter_promoter_Transfac) %>% as.data.frame() 
saveRDS(beta_TF_target_promoter_Transfac,file = 'Rds/beta_TF_target_promoter_Transfac.Rds')
beta_TF_target_promoter_Transfac <- readRDS('Rds/beta_TF_target_promoter_Transfac.Rds')
###enhancer analysis----
enhancer_predictions <-
    data.table::fread(
        '/data/kyh/islet_abc/abc_out/beta_ABC_out/EnhancerPredictionsFull.txt'
    ) %>%
    filter(isSelfPromoter != TRUE,!is.na(TargetGeneExpression)) %>%
    mutate(distanceToTSS = abs((start + end)/2 - TargetGeneTSS))
enhancer_predictions_gr <- makeGRangesFromDataFrame(enhancer_predictions,keep.extra.columns = T)

####Transfac
enhancer_predictions_gr_Transfac <- matchMotifs(motif_set_Transfac, enhancer_predictions_gr, genome = "hg38",out = "positions") 
beta_TF_target_enhancer_Transfac <- map(1:length(enhancer_predictions_gr_Transfac),function(x){
    overlapped_gr <- subsetByOverlaps(enhancer_predictions_gr,enhancer_predictions_gr_Transfac[[x]]) 
    overlapped_gr$motif <-  names(enhancer_predictions_gr_Transfac)[x]
    return(overlapped_gr)
},.progress = T) %>%  purrr::reduce(c) %>% as_tibble()  

beta_TF_target_enhancer_Transfac$TF <- TF_set[beta_TF_target_enhancer_Transfac$motif]
saveRDS(beta_TF_target_enhancer_Transfac,'Rds/beta_TF_target_enhancer_Transfac.Rds')

beta_TF_target_promoter_Transfac <- beta_TF_target_promoter_Transfac[,c(8,10,12,11)]
colnames(beta_TF_target_promoter_Transfac) <- c('TF','TargetGene','distanceToTSS','Type')
beta_TF_target_promoter_Transfac$Type <- rep('promoter',nrow(beta_TF_target_promoter_Transfac))

beta_TF_target_enhancer_Transfac <- beta_TF_target_enhancer_Transfac[,c(29,9,27)]
colnames(beta_TF_target_enhancer_Transfac) <- c('TF','TargetGene','distanceToTSS')
beta_TF_target_enhancer_Transfac$Type <- rep('enhancer',nrow(beta_TF_target_enhancer_Transfac))

beta_TF_target_unrefined_GRN_Transfac <- bind_rows(beta_TF_target_promoter_Transfac,beta_TF_target_enhancer_Transfac) 
beta_TF_target_unrefined_GRN_Transfac <- unique(beta_TF_target_unrefined_GRN_Transfac)
beta_TF_target_unrefined_GRN_Transfac <- separate_rows(beta_TF_target_unrefined_GRN_Transfac,TF,sep = ';') %>% unique()
saveRDS(beta_TF_target_unrefined_GRN_Transfac,'Rds/beta_TF_target_unrefined_GRN_Transfac.Rds')

##ATF4 beta cell track------
beta_ArchR <- loadArchRProject('beta_ArchR')

beta_narrow_peak <- rtracklayer::import('signac_macs2_temp/Beta_peaks.narrowPeak', format = 'narrowPeak')
beta_narrow_peak <- beta_narrow_peak[beta_narrow_peak@seqnames %in% paste0('chr',c(1:22,'X'))] 
seqlevels(beta_narrow_peak) <- paste0('chr',c(1:22,'X'))
beta_ArchR <- addPeakSet(beta_ArchR,peakSet = beta_narrow_peak,force = T)
beta_ArchR <- addPeakMatrix(beta_ArchR)
get_SummarizedExperiment <- function(SummarizedExperiment,return_gr = TRUE) {
    gr <- rowData(SummarizedExperiment) |> as.data.frame()
    assay_df <- reduce(SummarizedExperiment@assays@data,bind_cols) %>% 
        setNames(names(SummarizedExperiment@assays@data))
    gr <- gr %>% bind_cols(assay_df)
    return(gr)
}
#diff_peak
beta_diff_peak <- getMarkerFeatures(
    ArchRProj = beta_ArchR, 
    useMatrix = "PeakMatrix",
    groupBy = "status",
    testMethod = "wilcoxon",
    maxCells = 40000,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "T2D",
    bgdGroups = "ND"
)
beta_diff_peak_df <- get_SummarizedExperiment(beta_diff_peak)
beta_diff_peak_sign <- beta_diff_peak_df %>% filter(Pval < 0.01)
beta_diff_peak_sign_gr <- beta_diff_peak_sign |> makeGRangesFromDataFrame()
#risk_snp
risk_snp <- import.bed('bed/finemap.SLC30A8.rs13266634_hg38.bed')
risk_snp <- risk_snp %>% resize(width = 100,fix = 'center') 

##hic_loop
hic_loop <- data.table::fread('bed/EnhancerPredictions.bedpe') %>% filter(str_detect(V7,'SLC30A8'))
beta_TF_ATF4_link <- tibble(
    chr =  hic_loop$V1,
    start = (hic_loop$V2  + hic_loop$V3)/2,
    end = hic_loop$V5,
    type = 'ATF4-SLC30A8',
    value =  hic_loop$V8
)
hic_loop_gr <- map_dfr(1:nrow(beta_TF_ATF4_link),function(i){
    
    if (beta_TF_ATF4_link$start[i] < beta_TF_ATF4_link$end[i]) {
        df <- tibble_row(chr = beta_TF_ATF4_link$chr[i],
                         start = beta_TF_ATF4_link$start[i],
                         end = beta_TF_ATF4_link$end[i],
                         value = beta_TF_ATF4_link$value[i]) 
    }else{
        df <- tibble_row(chr = beta_TF_ATF4_link$chr[i],
                         start = beta_TF_ATF4_link$end[i],
                         end = beta_TF_ATF4_link$start[i],
                         value = beta_TF_ATF4_link$value[i]) 
    }
    
}) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
#ATF4
##Enhancers
hic_loop_gr_enhancer <- unique(hic_loop[,1:3]) %>% setNames(c('chr','start','end')) %>% makeGRangesFromDataFrame()
hic_loop_Transfac <- matchMotifs(motif_set_Transfac, GenomicRanges::reduce(hic_loop_gr_enhancer), genome = "hg38",out = "positions") 
hic_loop_enhancer_Transfac <- map(1:length(hic_loop_Transfac),function(x){
    if (length(hic_loop_Transfac[[x]]) >0 ) {
        overlapped_gr <- subsetByOverlaps(hic_loop_gr_enhancer,hic_loop_Transfac[[x]]) 
        overlapped_gr$motif <-  names(hic_loop_Transfac)[x]
        return(overlapped_gr)
    }
},.progress = T) %>%  purrr::reduce(c) %>% as_tibble()
hic_loop_enhancer_Transfac$TF <- TF_set[hic_loop_enhancer_Transfac$motif]
hic_loop_enhancer_Transfac <- separate_rows(hic_loop_enhancer_Transfac,TF,sep = ';') %>% unique()
hic_loop_enhancer_Transfac_ATF_4 <- hic_loop_enhancer_Transfac %>% filter(TF== 'ATF4') %>% 
    makeGRangesFromDataFrame() %>% GenomicRanges::reduce()
#promoter
beta_TF_target_promoter_Transfac <- readRDS('Rds/beta_TF_target_promoter_Transfac.Rds') %>%
    filter(TF == 'ATF4',symbol == 'SLC30A8') %>% makeGRangesFromDataFrame()
hic_loop_enhancer_Transfac_ATF_4 <- c(hic_loop_enhancer_Transfac_ATF_4,beta_TF_target_promoter_Transfac) %>% GenomicRanges::reduce()
hic_loop_enhancer_Transfac_ATF_4 <- subsetByOverlaps(hic_loop_enhancer_Transfac_ATF_4,beta_narrow_peak)


p <- plotBrowserTrack(
    ArchRProj = beta_ArchR, 
    groupBy = "status",
    useGroups = c('ND','T2D'),
    geneSymbol = c('SLC30A8'), 
    pal = c("#88c4e8","#db6968"),
    features = list(Peaks = getPeakSet(beta_ArchR),
                    Differential_Peaks = beta_diff_peak_sign_gr,
                    ATF4_binding = hic_loop_enhancer_Transfac_ATF_4),
    scCellsMax = 40000,
    upstream = 50000,
    downstream = 350000,
    loops = hic_loop_gr
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-ND_T2D-with-ABC.pdf", 
        ArchRProj = beta_ArchR, 
        addDOC = FALSE, width = 8, height = 6)


##benchmark with Ground-truth-----
###benchmark1------
##PDX-1(PDX1),MAFA,E47(TCF3), NEUROD1 known TF motif in INS promoter
bench1 <- beta_TF_target_unrefined_GRN_JASPAR %>% filter(str_detect(TF, 'PDX1|MAFA|TCF3|NEUROD1'),TargetGene == 'INS')
bench1 <- beta_TF_target_unrefined_GRN_Transfac %>% filter(str_detect(TF, 'PDX1|MAFA|TCF3|NEUROD1'),TargetGene == 'INS')
###benchmark2------
bench2 <-
    beta_TF_target_unrefined_GRN_JASPAR %>% filter(
        str_detect(TF, 'HNF4G|HNF4A|HNF1A|TCF4|NEUROD1|NFIX'),
        str_detect(TargetGene, 'HNF4G|HNF4A|HNF1A|TCF4|NEUROD1|NFIX')
    ) 
bench2 <-
    beta_TF_target_unrefined_GRN_Transfac %>% filter(
        str_detect(TF, 'HNF4G|HNF4A|HNF1A|TCF4|NEUROD1|NFIX'),
        str_detect(TargetGene, 'HNF4G|HNF4A|HNF1A|TCF4|NEUROD1|NFIX')
    ) %>% separate_rows(TF,sep = ';') %>%
     filter(
        str_detect(TF, 'HNF4G|HNF4A|HNF1A|TCF4|NEUROD1|NFIX'),
        str_detect(TargetGene, 'HNF4G|HNF4A|HNF1A|TCF4|NEUROD1|NFIX')
    ) 

#2000,500
NFIX(NF1A)-NEUROD1
NFIX(NF1A)-TCF4


beta_narrow_peak <- import('signac_macs2_temp/Beta_peaks.narrowPeak',format = 'narrowPeak')
beta_narrow_peak <- beta_narrow_peak[beta_narrow_peak@seqnames %in% paste0('chr',c(1:22,'X'))]
temp_promoter <- ChIPseeker::getPromoters(
    TxDb.Hsapiens.UCSC.hg38.knownGene,
    upstream = 5000,
    downstream = 5000,
    by = 'transcript'
) |> peak_anno()
temp_promoter_df <- as.data.frame(temp_promoter)
temp_promoter <- temp_promoter[which(temp_promoter$symbol == 'HNF4G')] |> GenomicRanges::reduce()

tf_HNF4A_var2 <- matchMotifs(motif_set_JASPAR2020[str_which(motif_symbol_JASPAR2020,'HNF4A')], 
                             temp_promoter,
                             genome = "hg38",out = "positions")
tf_HNF4A_var2 <- purrr::reduce(tf_HNF4A_var2,c) |> GenomicRanges::reduce()
#c('HNF4G','HNF4A','NEUROD1','TCF4')
p <- plotBrowserTrack(
    ArchRProj = alpha_beta, 
    groupBy = "cell_type", 
    geneSymbol = c('HNF4G'), 
    features = list(promoter = temp_promoter,
                    narrowpeak = beta_narrow_peak,
                    beta_peak = beta_peaks_gr,
                    HNF4A = tf_HNF4A_var2),
    upstream = 10000,
    downstream = 10000
)
saveplot(p$HNF4G,filenames = 'HNF4G_alpha_beta_track_Plot',width = 5,height = 4)

###benchmark3------
#FOXA1,FOXA2,HNF4A,HNF6(ONECUT1),PAX6,MAFA,PTF1A
#FOXA2,HNF4A,HNF6,PAX6,MAFA PRF1A USF1 to PDX1
bench3 <- beta_TF_target_enhancer_JASPAR2 %>% 
    filter(str_detect(TF, 'FOXA1|FOXA2|HNF4A|ONECUT1|PAX6|MAFA|PTF1A'),
           TargetGene == 'PDX1')
bench3 <- beta_TF_target_enhancer_Transfac %>% 
    filter(str_detect(TF, 'FOXA1|FOXA2|HNF4A|ONECUT1|PAX6|MAFA|PTF1A'),
           TargetGene == 'PDX1')  %>% separate_rows(TF,sep = ';') %>%
    filter(str_detect(TF, 'FOXA1|FOXA2|HNF4A|ONECUT1|PAX6|MAFA|PTF1A'),
           TargetGene == 'PDX1') %>% distinct(.keep_all = T)
###benchmark4------
#TF:BMAL1 not in JASPER database

#human CLOCK target
CLOCK<- c("PARD6A","VEGFC", "PARD6G","ITGB1", "CDC42", "HGF", "SIPA1L2","DOCK4", "PRKD1",
 "ADCY2", "MAGI2", "MAPK1", "BCAR1", "PARD6B","RAPGEF1", "EFNA5", "AKT2","GNAI3",
"NGFR","F2RL3", "VEGFB", "KIT", "GNAI1", "DUSP7", "DUSP6", "CDC42", "MAP3K4" ,
 "MAP3K7","TAB2","MAX","MAPK8", "DUSP1", "MECOM", "AKT2","GNG12", "MAPK1",
"RPS6KA1", "JUN", "MAP4K2","DUSP5", "MKNK2", "TGFB2", "PPP3CA","CRY2","FBXL3",
 "BHLHE40", "BHLHE41","PER1","ABCC8", "ADCY2", "KCNJ11","CAMK2D","KCNMA1","ATP1B3" ,
 "MAPK1", "MAPK8", "ABCC8", "KCNJ11","DGKG","INPP4B","PTEN","GOSR2")

bench4 <- beta_TF_target_unrefined_GRN_JASPAR %>% 
    filter(str_detect(TF, 'CLOCK'))
length(which(CLOCK %in% bench4$TargetGene))/length(CLOCK)

bench4 <- beta_TF_target_unrefined_GRN_Transfac %>% 
    filter(str_detect(TF, 'CLOCK'))
length(which(CLOCK %in% bench4$TargetGene))/length(CLOCK)

###benchmark5------
bench5_ep <- data.table::fread('bed/Enhancer-promoter_hg38.bed')
colnames(bench5_ep) <- c('chr','start','end','gene','frame','starnd')
bench5_gr <- makeGRangesFromDataFrame(bench5_ep,keep.extra.columns = T)    
    
length(subsetByOverlaps(reduce(bench5_gr),reduce(enhancer_predictions_gr))) / length(reduce(bench5_gr))

length(subsetByOverlaps(reduce(bench5_gr),reduce(enhancer_predictions_gr))) / length(reduce(enhancer_predictions_gr))
