library(Seurat)
library(tidyverse)
library(xgboost)
library(Matrix)
library(ggsignif)
library(cowplot)

all_nd = 12+18+5+8+5+33+6+12+7+6
all_nd_cell_num = 6297+227+10620 + 228 + 11247 + 24296 + 150 + 117 +6640 +11356
t2d = 7+21+4+6+4+6+9
all_t2d_cell_num = 89 + 9902 + 130 + 278 + 5827 +13795
##function------
gene_length <- readRDS('Rds/human_ensembl2symbol_with_length.Rds')
get_tpm <- function(count_mtx){
    tpm <- function(counts,len) {
        x <- counts/len
        return(t(t(x)*1e6/colSums(x)))
    }
    count_mtx$gene_name <- rownames(count_mtx)
    count_mtx <- count_mtx %>% left_join(gene_length[3:4],by = c('gene_name')) %>% 
        filter(!is.na(length))
    count_tpm <- tpm(count_mtx[,str_which(colnames(count_mtx),'sample|multiome')],count_mtx$length)
    rownames(count_tpm) <- count_mtx$gene_name
    return(count_tpm)
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
#
train_classifier <- function(donor,seurat_obj,count_mtx = NULL,remove = c('donor','cell','rate'),downsample = FALSE,ncore = 60){
    label <- if_else(seurat_obj$status == 'ND',0,1)
    names(label) <-  colnames(seurat_obj)
    ##
    
    
    train_data <- count_mtx[,which(seurat_obj$donor != donor)]
    train_label <- label[which(seurat_obj$donor != donor)]
    
    if (downsample) {
        min_type <-  as_tibble(table(train_label)) %>% slice_min(n) %>% pull(train_label) %>% as.numeric()
        min_num <- as_tibble(table(train_label)) %>% slice_min(n) %>% pull(n) %>% as.numeric()
        index <-  c(which(train_label == min_type),sample(which(train_label != min_type),min_num)) 
        train_data <- train_data[,index]
        train_label <- train_label[index]
        cat('Down sample',min_num, 'cells from',if_else(min_type == 1,'ND','T2D'),'\n')
    }
    
    classifier <- xgboost(
        data = Matrix::t(train_data),
        label = train_label,
        max.depth = 6,
        eta = 0.2,
        nthread = ncore,
        nrounds = 100,
        verbose = 0,
        objective = "binary:logistic"
    )
    #
    test_data <- count_mtx[,which(seurat_obj$donor == donor)]
    test_label <- label[which(seurat_obj$donor == donor)]
    
    pred_label <- predict(classifier, Matrix::t(test_data))
    if (unique(test_label) == 1) {
        accuracy_cell <- test_label[which(pred_label >= 0.5)]
    }else{
        accuracy_cell <- test_label[which(pred_label <= 0.5)]
        
    }
    
    cat('Prediction accuracy rate:',length(accuracy_cell)/length(test_label),'\n')
    if (remove == 'donor') {
        if (length(accuracy_cell)/length(test_label) < 0.15) {
            return(donor)
        }
    }else if(remove == 'cell'){
        return(accuracy_cell)
    }else{
        df <- tibble(donor = donor,
                     accuracy_rate = length(accuracy_cell)/length(test_label))
        return(df)
    }
}
#
get_singlecell_gene_exp <- function(i,gene_vec, seuart_list){
    gene_exp <- FetchData(seuart_list[[i]],c(gene_vec,'donor','status'))
    gene_exp$dataset <- names(seuart_list)[i]
    return(gene_exp)
}
get_DEGs <- function(pseudobulk){
    donor_num <- length(str_which(colnames(pseudobulk),'ND')) + length(str_which(colnames(pseudobulk),'T2D'))
    if (donor_num != ncol(pseudobulk)) {
        stop('Detect inconsisted donor numbers')
    }
    
    coldata <- data.frame(row.names=colnames(pseudobulk),condition = factor(if_else(str_detect(colnames(pseudobulk),'ND'),'ND','T2D')))
    deseq2.obj <- DESeqDataSetFromMatrix(as.matrix(pseudobulk), colData = coldata,design=~condition)
    deseq2.obj  <- DESeq(deseq2.obj)
    deseq2.obj.res <- results(deseq2.obj)
    deseq2.obj.res.df <-  as.data.frame(deseq2.obj.res)  %>% rownames_to_column(var = 'gene')
    return(deseq2.obj.res.df)
}
##calculate cor------
calculate_cor <- function(x,y,method = c("pearson","spearman")){
    cor_res <- cor.test(log1p(x),log1p(y),exact = F,method = method) 
    return(list(cor = as.numeric(cor_res$estimate),pvalue = cor_res$p.value))
}
get_gene_pair_cor <- function(gene1,gene2,method = c("pearson","spearman"), pseudobulk,sample_ratio = NULL){
    if (is.null(pseudobulk)) {
        stop('Must with pseudobulk')
    }
    df <- tibble(gene1 = pseudobulk[gene1,],gene2 = pseudobulk[gene2,],sample = colnames(pseudobulk))
    df <- df %>% filter(gene1 > 0.1 & gene2 > 0.1)
    colnames(df)[1:2] <- c(gene1,gene2)
    
    if (nrow(df) < 30) {
        cor_res <- list(cor = 0,pvalue = 0)
    }else{
        cor_res <- calculate_cor(pseudobulk[gene1,],pseudobulk[gene2,],method = method)
    }
    cor_res <- calculate_cor(pseudobulk[gene1,],pseudobulk[gene2,],method = method)
    return(cor_res)
}
batch_calculate_cor <- function(rep_num, GRN,pseudobulk = NULL,sample_ratio = NULL,ncore = NULL) {
    ##pearson
    cor_pearson = paste0('cor_pearson_rep', rep_num)
    pvalue_pearson = paste0('pvalue_pearson_rep', rep_num)
    pvalue_pearson_adj = paste0('pvalue_pearson_adj_rep', rep_num)
    if (!is.null(sample_ratio)) {
        col_index <- sample(1:ncol(pseudobulk),ceiling(ncol(pseudobulk) * sample_ratio))
        pseudobulk <- pseudobulk[, col_index]
    }
    if (!is.null(ncore)) {
        require(furrr)
        options(future.globals.maxSize= 50*1024*1024^2)
        plan(multisession, workers = ncore) 
        suppressWarnings({TF_target_cor_pearson <-
            future_map2(
                GRN$TF,
                GRN$TargetGene,
                function(TF, target,pseudobulk,sample_ratio){
                    res <- get_gene_pair_cor(
                        TF,
                        target,
                        method = 'pearson',
                        pseudobulk = pseudobulk,
                        sample_ratio = sample_ratio
                    )
                    return(res)
                },pseudobulk = pseudobulk,sample_ratio = sample_ratio)})
    }else{
        suppressWarnings({TF_target_cor_pearson <-
            map2(
                GRN$TF,
                GRN$TargetGene,
                function(TF, target,pseudobulk,sample_ratio){
                    res <- get_gene_pair_cor(
                        TF,
                        target,
                        method = 'pearson',
                        pseudobulk = pseudobulk,
                        sample_ratio = sample_ratio
                    )
                    return(res)
                },pseudobulk = pseudobulk,sample_ratio = sample_ratio,.progress = T
            )})
    }
    
    GRN <- GRN %>% mutate(
        !!cor_pearson := map_dbl(TF_target_cor_pearson, function(x) {
            return(x$cor)
        }),
        !!pvalue_pearson := map_dbl(TF_target_cor_pearson, function(x) {
            return(x$pvalue)
        }),
        !!pvalue_pearson_adj := p.adjust(.data[[pvalue_pearson]], method = "BH")
    )
    
    ##spearman
    cor_spearman = paste0('cor_spearman_rep', rep_num)
    pvalue_spearman = paste0('pvalue_spearman_rep', rep_num)
    pvalue_spearman_adj = paste0('pvalue_spearman_adj_rep', rep_num)
    if (!is.null(ncore)) {
        suppressWarnings({TF_target_cor_spearman <-
            future_map2(GRN$TF,
                        GRN$TargetGene,
                        function(TF, target,pseudobulk,sample_ratio) {
                            res <- get_gene_pair_cor(
                                TF,
                                target,
                                method = 'spearman',
                                pseudobulk = pseudobulk,
                                sample_ratio = sample_ratio
                            )
                            return(res)
                        }, pseudobulk = pseudobulk,sample_ratio = sample_ratio)
        })
    }else{
        suppressWarnings({
            TF_target_cor_spearman <-
                map2(GRN$TF,GRN$TargetGene,
                     function(TF, target,pseudobulk,sample_ratio) {
                         res <- get_gene_pair_cor(
                             TF,
                             target,
                             method = 'spearman',
                             pseudobulk = pseudobulk,
                             sample_ratio = sample_ratio
                         )
                         return(res)
                     }, 
                     pseudobulk = pseudobulk,sample_ratio = sample_ratio,.progress = T)
        })
    }
    GRN <- GRN %>%
        mutate(
            !!cor_spearman := map_dbl(TF_target_cor_spearman, function(x) {
                return(x$cor)
            }),
            !!pvalue_spearman := map_dbl(TF_target_cor_spearman, function(x) {
                return(x$pvalue)
            }),
            !!pvalue_spearman_adj := p.adjust(.data[[pvalue_spearman]], method = "BH")
        )
    return(GRN)
}
##Plot cor------
plot_cor <- function(gene1,gene2,method = c("pearson","spearman"),pseudobulk = NULL,sample_ratio = NULL){
    require(ggpubr)
    if (!is.null(sample_ratio)) {
        col_index <- sample(1:ncol(pseudobulk),
                            ceiling(ncol(pseudobulk) * sample_ratio))
        pseudobulk <- pseudobulk[, col_index]
    }else{
        pseudobulk <- pseudobulk
    }
    get_exp <- function(gene,pseudobulk){
        as.numeric(pseudobulk[gene,])
    }
    df <- tibble(get_exp(gene1,pseudobulk),get_exp(gene2,pseudobulk))
    colnames(df) <- c(gene1,gene2)
    
    sp <- ggscatter(df, x = gene1, y = gene2,
                    add = "reg.line", 
                    add.params = list(color = "blue", fill = "lightgray"), 
                    conf.int = TRUE 
    ) + 
        stat_cor(
            aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
            label.x = 1,
            method = method
        ) +
        scale_x_continuous(trans = 'log1p') +
        scale_y_continuous(trans = 'log1p')
    theme_bw() 
    return(sp)
}
plot_cor2 <- function(df,pseudobulk){
    tf = df[1,1][[1]]
    target =  df[1,2][[1]]
    if (abs(df[1,4][[1]]) > abs(df[1,7][[1]])) {
        method = 'pearson'
    }else{
        method = 'spearman'
    }
    p1 <- plot_cor(tf,target,method,pseudobulk = pseudobulk)
    return(p1)
}
plot_cor_temp <- function(gene1,gene2,method = c("pearson","spearman"),pseudobulk = NULL,group_by = c('type','disease','dataset'),sample_ratio = NULL){
    require(ggpubr)
    if (!is.null(sample_ratio)) {
        col_index <- sample(1:ncol(pseudobulk),
                            ceiling(ncol(pseudobulk) * sample_ratio))
        pseudobulk <- pseudobulk[, col_index]
    }else{
        pseudobulk <- pseudobulk
    }
    
    
    get_exp <- function(gene,pseudobulk){
        as.numeric(pseudobulk[gene,])
    }
    
    
    df <- tibble(gene1 = get_exp(gene1,pseudobulk),gene2 = get_exp(gene2,pseudobulk),sample = colnames(pseudobulk))
    df$type <- if_else(df$sample %in% sample_10x,'10X','Non-10x')
    df$disease <- if_else(df$sample %in% sample_ND,'ND','T2D')
    df$dataset <- case_when(
        str_detect(df$sample,'Biobank_multiome')  ~ 'Biobank multiome',
        str_detect(df$sample,'sample1_Donor')  ~ 'sample1',
        str_detect(df$sample,'sample10')  ~ 'sample10',
        str_detect(df$sample,'sample12')  ~ 'sample12',
        str_detect(df$sample,'sample14')  ~ 'sample14',
        str_detect(df$sample,'sample16')  ~ 'sample16',
        str_detect(df$sample,'sample3')  ~ 'sample3',
        str_detect(df$sample,'sample4')  ~ 'sample4',
        str_detect(df$sample,'sample5')  ~ 'sample5',
        str_detect(df$sample,'sample6')  ~ 'sample6',
        str_detect(df$sample,'sample7')  ~ 'sample7'
    )
    df <- df %>% filter(gene1 > 0.5 & gene2 > 0.5)
    colnames(df)[1:2] <- c(gene1,gene2)
    sp <- ggscatter(df, x = gene1, y = gene2,
                    add = "reg.line", color = group_by, palette = clustcol[1:length(unique(df[[group_by]]))],
                    add.params = list(color = "blue", fill = "lightgray"), 
                    conf.int = TRUE 
    ) 
    
    if (gene1 == 'INS') {
        sp <- sp +  scale_x_log10() +
            scale_y_log10() +
            stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
                     method = method)
        theme_bw()
    } else{
        sp <- sp +  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                             method = method) +
            scale_x_continuous(trans = 'log1p') +
            scale_y_continuous(trans = 'log1p')
        theme_bw()
    }
    
    return(sp)
}
plot_singlecell_gene_exp <- function(df) {
    gene_name <-  colnames(df)[1]
    colnames(df)[1] <- 'gene'
    
    nd <- df %>% filter(status == 'ND') %>% pull(gene)
    t2d <- df %>% filter(status == 'T2D') %>% pull(gene)
    wilcoxtest <- wilcox.test(nd,t2d)
    pvalue <- signif(wilcoxtest$p.value, 2)
    if (pvalue < 1e-22) {
        pvalue <- 'p < 1e-22'
    }else{
        pvalue <- paste0('p = ',pvalue)
    }
    
    p1 <- ggplot(data= df, aes(x=status,y=gene))+
        geom_violin(aes(fill=status),colour= NA,width=0.5) +
        scale_fill_manual(values = c("#475eb1", "#bc192b")) +
        geom_boxplot(width=0.1, position = position_dodge(0.9))+
        ylim(0,max(df[1])*1.2) +
        theme_bw()+
        labs(x = NULL,y = 'Normalized expression level',
             #title = paste0(gene_name,' ',unique(df[['dataset2']]))
             title = unique(df[['dataset2']])
        ) +
        guides(fill = 'none') +
        geom_signif(
            comparisons = list(c("ND", "T2D")),
            tip_length = 0,
            margin_top = 0.05,
            map_signif_level = function(p,gene = gene_name,real_p = pvalue){
                star <- case_when(
                    p < 0.001 ~ '***',
                    p < 0.01 ~ '**',
                    p <= 0.05 ~ '*',
                    p > 0.05 ~ 'ns'
                )
                return(paste0(real_p,'\n',star))
            }
        ) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(color = 'black',size = rel(1.5)),
              panel.grid = element_blank())
    return(p1)
}


##Paper reproduce part-----
t2d_dataset <- readRDS('Rds/t2d_dataset.Rds')


###Start Machine Learning 


####sample3------
sample3 <- t2d_dataset$sample3
sample3$status <- if_else(str_detect(sample3$donor,'ND'),'ND','T2D')
sample3$cell_type <- 'Beta'

sample3_counts <- as(edgeR::cpm(sample3@assays$RNA@counts),'dgCMatrix')
sample3_counts <- sample3_counts[which(rowMeans(sample3_counts) > 1),]
sample3_bad_donor <- map(unique(sample3$donor),train_classifier,seurat_obj = sample3,downsample = T,
                         count_mtx = sample3_counts,remove = 'donor',ncore = 1,.progress = T)
sample3_no_bad_donor <- subset(sample3,donor %in% unlist(sample3_bad_donor),invert =T)
#3 round
while (!is.null(unlist(sample3_bad_donor))) {
    sample3_counts <- as(edgeR::cpm(sample3_no_bad_donor@assays$RNA@counts),'dgCMatrix')
    sample3_counts <- sample3_counts[which(rowMeans(sample3_counts) > 1),]
    sample3_bad_donor <- map(unique(sample3_no_bad_donor$donor),train_classifier,seurat_obj = sample3_no_bad_donor,downsample = T,
                             count_mtx = sample3_counts,remove = 'donor',ncore = 70,.progress = T)
    sample3_no_bad_donor <- subset(sample3_no_bad_donor,donor %in% unlist(sample3_bad_donor),invert =T)
}
###sample14------
sample14 <- t2d_dataset$sample14
sample14$status <- if_else(str_detect(sample14$donor,'ND'),'ND','T2D')

sample14_counts <- as(edgeR::cpm(sample14@assays$RNA@counts),'dgCMatrix')
sample14_counts <- sample14_counts[which(rowMeans(sample14_counts) > 1),]
sample14_bad_donor <- map(unique(sample14$donor),train_classifier,seurat_obj = sample14,downsample = T,
                          count_mtx = sample14_counts,remove = 'donor',ncore = 70,.progress = T)
sample14_no_bad_donor <- subset(sample14,donor %in% unlist(sample14_bad_donor),invert =T)
#3 round
while (!is.null(unlist(sample14_bad_donor))) {
    sample14_counts <- as(edgeR::cpm(sample14_no_bad_donor@assays$RNA@counts),'dgCMatrix')
    sample14_counts <- sample14_counts[which(rowMeans(sample14_counts) > 1),]
    sample14_bad_donor <- map(unique(sample14_no_bad_donor$donor),train_classifier,seurat_obj = sample14_no_bad_donor,downsample = T,
                              count_mtx = sample14_counts,remove = 'donor',ncore = 70,.progress = T)
    sample14_no_bad_donor <- subset(sample14_no_bad_donor,donor %in% unlist(sample14_bad_donor),invert =T)
}
###sample16-------
sample16 <- t2d_dataset$sample16
sample16$status <- if_else(str_detect(sample16$donor,'ND'),'ND','T2D')

sample16_counts <- as(edgeR::cpm(sample16@assays$RNA@counts),'dgCMatrix')
sample16_counts <- sample16_counts[which(rowMeans(sample16_counts) > 1),]
sample16_bad_donor <- map(unique(sample16$donor),train_classifier,seurat_obj = sample16,downsample = T,
                          count_mtx = sample16_counts,remove = 'donor',ncore = 70,.progress = T)
sample16_no_bad_donor <- subset(sample16,donor %in% unlist(sample16_bad_donor),invert =T)
#3 round
while (!is.null(unlist(sample16_bad_donor))) {
    sample16_counts <- as(edgeR::cpm(sample16_no_bad_donor@assays$RNA@counts),'dgCMatrix')
    sample16_counts <- sample16_counts[which(rowMeans(sample16_counts) > 1),]
    sample16_bad_donor <- map(unique(sample16_no_bad_donor$donor),train_classifier,seurat_obj = sample16_no_bad_donor,downsample = T,
                              count_mtx = sample16_counts,remove = 'donor',ncore = 70,.progress = T)
    sample16_no_bad_donor <- subset(sample16_no_bad_donor,donor %in% unlist(sample16_bad_donor),invert =T)
}
###sample17-------
sample17 <- t2d_dataset$sample17
sample17$donor <- paste0(sample17$sample,'_',sample17$donor,'_',sample17$status )
sample17_counts <- as(edgeR::cpm(sample17@assays$RNA@counts),'dgCMatrix')
sample17_counts <- sample17_counts[which(rowMeans(sample17_counts) > 1),]
sample17_bad_donor <- map(unique(sample17$donor),train_classifier,seurat_obj = sample17,
                          count_mtx = sample17_counts,remove = 'donor',ncore = 70,.progress = T)
sample17_no_bad_donor <- subset(sample17,donor %in% unlist(sample17_bad_donor),invert =T)
#0 round
while (!is.null(unlist(sample17_bad_donor))) {
    sample17_counts <- as(edgeR::cpm(sample17_no_bad_donor@assays$RNA@counts),'dgCMatrix')
    sample17_counts <- sample17_counts[which(rowMeans(sample17_counts) > 1),]
    sample17_bad_donor <- map(unique(sample17_no_bad_donor$donor),train_classifier,seurat_obj = sample17_no_bad_donor,downsample = T,
                              count_mtx = sample17_counts,remove = 'donor',ncore = 70,.progress = T)
    sample17_no_bad_donor <- subset(sample17_no_bad_donor,donor %in% unlist(sample17_bad_donor),invert =T)
}

##HPAP (sample12) data------
##start here 
sample12 <- t2d_dataset$sample12
sample12$status <- if_else(str_detect(sample12$donor,'ND'),'ND','T2D')

#
sample12_counts <- as(edgeR::cpm(sample12@assays$RNA@counts),'dgCMatrix')
sample12_counts <- sample12_counts[which(rowMeans(sample12_counts) > 1),]
sample12_bad_donor <- map(unique(sample12$donor),train_classifier,seurat_obj = sample12,count_mtx = sample12_counts,remove = 'donor',ncore = 70,.progress = T)
#

sample12_no_bad_donor <- subset(sample12,donor %in% unlist(sample12_bad_donor),invert =T)
#Two round
while (!is.null(unlist(sample12_bad_donor))) {
    sample12_counts <- as(edgeR::cpm(sample12_no_bad_donor@assays$RNA@counts),'dgCMatrix')
    sample12_counts <- sample12_counts[which(rowMeans(sample12_counts) > 1),]
    sample12_bad_donor <- map(unique(sample12_no_bad_donor$donor),train_classifier,seurat_obj = sample12_no_bad_donor,downsample = T,
                              count_mtx = sample12_counts,remove = 'donor',ncore = 70,.progress = T)
    sample12_no_bad_donor <- subset(sample12_no_bad_donor,donor %in% unlist(sample12_bad_donor),invert =T)
}

#####multiome beta1 beta2------
multiome_beta <- t2d_dataset$Biobank_multiome

#
multiome_counts <- as(edgeR::cpm(multiome_beta@assays$RNA@counts),'dgCMatrix')
multiome_counts <- multiome_counts[which(rowMeans(multiome_counts) > 1),]
multiome_bad_donor <- map(unique(multiome_beta$donor),train_classifier,seurat_obj = multiome_beta,count_mtx = multiome_counts,remove = 'donor',ncore = 70,.progress = T)
#
multiome_no_bad_donor <- subset(multiome_beta,donor %in% unlist(multiome_bad_donor),invert =T)
while (!is.null(unlist(multiome_bad_donor))) {
    multiome_counts <- as(edgeR::cpm(multiome_no_bad_donor@assays$RNA@counts),'dgCMatrix')
    multiome_counts <- multiome_counts[which(rowMeans(multiome_counts) > 1),]
    multiome_bad_donor <- map(unique(multiome_no_bad_donor$donor),train_classifier,seurat_obj = multiome_no_bad_donor,
                              count_mtx = multiome_counts,remove = 'donor',ncore = 70,.progress = T)
    multiome_no_bad_donor <- subset(multiome_no_bad_donor,donor %in% unlist(multiome_bad_donor),invert =T)
}

t2d_dataset_no_bad_donor <- list(
    Biobank_multiome = multiome_no_bad_donor,
    sample3 = sample3_no_bad_donor,
    sample16 = sample16_no_bad_donor,
    sample12 = sample12_no_bad_donor,
    sample14 = sample14_no_bad_donor,
    sample17 = sample17_no_bad_donor
)
saveRDS(t2d_dataset_no_bad_donor,'Rds/t2d_dataset_no_bad_donor.Rds')

##
#t2d_dataset_no_bad_donor <- readRDS('Rds/t2d_dataset_no_bad_donor.Rds')
all_dataset_beta_DEGs_MAST <- map(1:length(t2d_dataset_no_bad_donor),function(i,dataset){
    Idents(dataset[[i]]) <- dataset[[i]]$status
    beta_DEGs <- FindMarkers(dataset[[i]],test.use = 'MAST',logfc.threshold = 0.1,ident.1 = 'T2D',ident.2 = 'ND') %>%
        rownames_to_column('gene') %>% mutate(dataset = names(dataset)[i])
    return(beta_DEGs)
},dataset = t2d_dataset_no_bad_donor) %>% setNames(names(t2d_dataset_no_bad_donor))
saveRDS(all_dataset_beta_DEGs_MAST,'all_dataset_beta_DEGs_MAST.Rds')

donor_info <- map_dfr(t2d_dataset_no_bad_donor,function(x){as.data.frame(table(x$status,x$donor))}) %>% filter(Freq > 0)

#donor_info$Var1 |> table() ;ND 68 T2D 30 exclude  ND 14  T2D 18
map_int(t2d_dataset_no_bad_donor,function(x){
    table(x$status)[2]
}) |> sum()
#23250 T2D Beta cells
#30482 ND beta cells

ncol(all_dataset$sample1)  + ncol(all_dataset$sample4) + ncol(all_dataset$sample6) + ncol(all_dataset$sample10)
#28452


all_dataset_beta_DEGs_MAST <- readRDS('Rds/all_dataset_beta_DEGs_MAST.Rds')
all_dataset_beta_DEGs_up <- map(all_dataset_beta_DEGs_MAST,function(x){
    x %>% filter(avg_log2FC > 0,p_val_adj <= 0.05) 
})
all_dataset_beta_DEGs_down <- map(all_dataset_beta_DEGs_MAST,function(x){
    x %>% filter(avg_log2FC < 0,p_val_adj <= 0.05) 
})

get_freq <- function(dataset,gene){
    res <-  if_else(gene %in% dataset$gene,1,0)
    df <- tibble(res)
    colnames(df) <- unique(dataset$dataset)
    return(df)
}
up_gene_freq <- map_dfr(unique(bind_rows(all_dataset_beta_DEGs_up)$gene),function(x){
    df <- map_dfc(all_dataset_beta_DEGs_up,get_freq,gene = x) %>% mutate(gene = x,.before = 1)  
    df <- df %>% mutate(sum = pmap_dbl(dplyr::select(df,where(is.numeric)) , ~ sum(c(...))))
    return(df)
},.progress = T) %>% mutate(type = 'Up regulated')
down_gene_freq <- map_dfr(unique(bind_rows(all_dataset_beta_DEGs_down)$gene),function(x){
    df <- map_dfc(all_dataset_beta_DEGs_down,get_freq,gene = x) %>% mutate(gene = x,.before = 1)  
    df <- df %>% mutate(sum = pmap_dbl(dplyr::select(df,where(is.numeric)) , ~ sum(c(...))))
    return(df)
},.progress = T) %>% mutate(type = 'Down regulated')
beta_DEGs_freq <- bind_rows(up_gene_freq,down_gene_freq) %>% 
    group_by(type)  %>% arrange(desc(sum))  %>% mutate(n = 1:n()) 
writexl::write_xlsx(beta_DEGs_freq,'DEGs_MSAT.xlsx')

beta_DEGs_freq_num <- beta_DEGs_freq %>% group_by(type,sum) %>% mutate(n = 1:n()) %>% summarise(freq_num = max(n))


p1 <- ggplot(beta_DEGs_freq_num,aes(sum,y = freq_num,fill = type)) +
    geom_col(position = 'dodge') +
    scale_y_log10() +
    labs(x = 'Frequncey',y = 'Gene number',fill = '',title = 'Overlapped DEG numbers') +
    geom_text(aes(y = freq_num + 0.25 * (freq_num),label = freq_num),position = position_dodge(width = 1)) +
    scale_fill_manual(values = c("#aed4e9","#93cc82")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())         
saveplot(p1,filenames = 'DEG_freqency_MAST',width = 3.2,height = 2)
data.table::fwrite(beta_DEGs_freq[,1:9],file = 'beta_DEGs_freq_MAST.xls',sep = '\t')



