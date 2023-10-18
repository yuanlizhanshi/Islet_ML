library(tidyverse)
library(furrr)
##function------
calculate_cor <- function(x,y,method = c("pearson","spearman")){
    cor_res <- cor.test(log1p(x),log1p(y),exact = F,method = method) 
    return(list(cor = cor_res$estimate,pvalue = cor_res$p.value))
}
get_gene_pair_cor <- function(gene1,gene2,method = c("pearson","spearman"), pseudobulk,sample_ratio = NULL){
    if (is.null(pseudobulk)) {
        stop('Must with pseudobulk')
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
        col_index <- sample(1:ncol(beta_pseudobulk_tpm),ceiling(ncol(beta_pseudobulk_tpm) * sample_ratio))
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
                     }, pseudobulk = pseudobulk,sample_ratio = sample_ratio,.progress = T)
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

plot_cor <- function(gene1,gene2,method = c("pearson","spearman"),pseudobulk = NULL,sample_ratio = NULL){
    require(ggpubr)
    if (!is.null(sample_ratio)) {
        col_index <- sample(1:ncol(beta_pseudobulk_tpm),
                            ceiling(ncol(beta_pseudobulk_tpm) * sample_ratio))
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
##TPM----
beta_pseudobulk_tpm <- readRDS('GRN/beta_pseudobulk_tpm.Rds')
beta_pseudobulk_tpm_mean <- beta_pseudobulk_tpm  %>% as.data.frame() %>% 
    mutate(TPM_mean = rowMeans(.),gene = rownames(beta_pseudobulk_tpm)) %>%
    dplyr::select(gene,TPM_mean) %>% 
    filter(TPM_mean >= 1)


###10 iteration-----
for (i in 1:10) {
    cat('Calculate rep', i, '\n')
    beta_TF_target_unrefined_GRN_Transfac_unique <-
        batch_calculate_cor(
            i,
            beta_TF_target_unrefined_GRN_Transfac_unique,
            beta_pseudobulk_tpm,
            sample_ratio = 0.4,
            ncore = 60
        )
}
#60 core in 5 hours
#Using python will be more fast i core in 4 h
beta_TF_target_unrefined_GRN_Transfac_unique <- beta_TF_target_unrefined_GRN_Transfac_unique %>% rowwise() %>% mutate(
    mean_peason_cor = mean(c_across(contains('cor_pearson'))),
    mean_peason_cor_padj = mean(c_across(contains('pvalue_pearson_adj'))),
    mean_spearman_cor = mean(c_across(contains('cor_spearman'))),
    mean_spearman_cor_padj = mean(c_across(contains('pvalue_spearman_adj'))), 
) 
##Elegant but slow (11h)

beta_TF_target_unrefined_GRN_Transfac_unique <- 
    beta_TF_target_unrefined_GRN_Transfac_unique %>%
    mutate(mean_peason_cor = pmap_dbl(select(beta_TF_target_unrefined_GRN_Transfac_unique,contains('cor_pearson')) , ~ mean(c(...)),.progress = T),
           mean_peason_cor_padj = pmap_dbl(select(beta_TF_target_unrefined_GRN_Transfac_unique,contains('pvalue_pearson_adj')) , ~ mean(c(...)),.progress = T),
           mean_spearman_cor = pmap_dbl(select(beta_TF_target_unrefined_GRN_Transfac_unique,contains('cor_spearman')) , ~ mean(c(...)),.progress = T),
           mean_spearman_cor_padj = pmap_dbl(select(beta_TF_target_unrefined_GRN_Transfac_unique,contains('pvalue_spearman_adj')) , ~ mean(c(...)),.progress = T))
##More faster (8 min)
saveRDS(beta_TF_target_unrefined_GRN_Transfac_unique,file = 'GRN/beta_TF_target_unrefined_GRN_Transfac_unique_10_iteration.Rds')  

beta_TF_target_unrefined_GRN_Transfac_unique_cor7 <- beta_TF_target_unrefined_GRN_Transfac_unique %>%    
    dplyr::select(TF,TargetGene,Type,mean_peason_cor,mean_peason_cor_padj,mean_spearman_cor,mean_spearman_cor_padj) %>% filter(
    abs(mean_peason_cor) >= 0.7 | abs(mean_spearman_cor) >= 0.7,
    mean_peason_cor_padj <= 0.05 |mean_spearman_cor_padj <= 0.05
)
saveRDS(beta_TF_target_unrefined_GRN_Transfac_unique_cor7,file = 'GRN/beta_TF_target_unrefined_GRN_Transfac_cor7.Rds')  


