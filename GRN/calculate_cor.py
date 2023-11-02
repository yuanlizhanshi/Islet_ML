import pandas as pd
import scipy.stats as stats
import numpy as np
import math
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import warnings
from functools import partial
import sys

warnings.filterwarnings("ignore")

def fast_calculate_cor(x,y, method = "pearson|spearman"):
    x = np.asarray(x)
    y = np.asarray(y)
    if method == 'pearson':
        n = len(x)
        sum_x = np.sum(x)
        sum_y = np.sum(y)
        sum_xy = np.sum(x * y)
        sum_xsq = np.sum(x**2)
        sum_ysq = np.sum(y**2)
        
        numerator = (n * sum_xy) - (sum_x * sum_y)
        denominator = np.sqrt(((n * sum_xsq) - (sum_x**2)) * ((n * sum_ysq) - (sum_y**2)))
        
        correlation_coefficient = numerator / denominator
        
        # Calculate the degrees of freedom
        df = n - 2
        
        # Calculate the t-statistic
        t_statistic = correlation_coefficient * np.sqrt(df / (1 - correlation_coefficient**2))
        
        # Calculate the p-value
        p_value = 2 * (1 - stats.t.cdf(np.abs(t_statistic), df))

        return correlation_coefficient, p_value
    elif method == 'spearman':
        ranks_x = np.argsort(x).argsort()
        ranks_y = np.argsort(y).argsort()
        
        n = len(x)
        sum_ranks_x = np.sum(ranks_x)
        sum_ranks_y = np.sum(ranks_y)
        sum_ranks_xy = np.sum(ranks_x * ranks_y)
        sum_ranks_xsq = np.sum(ranks_x**2)
        sum_ranks_ysq = np.sum(ranks_y**2)
        
        numerator = (n * sum_ranks_xy) - (sum_ranks_x * sum_ranks_y)
        denominator = np.sqrt(((n * sum_ranks_xsq) - (sum_ranks_x**2)) * ((n * sum_ranks_ysq) - (sum_ranks_y**2)))
        
        correlation_coefficient = numerator / denominator
        
        # Calculate the degrees of freedom
        df = n - 2
        
        # Calculate the t-statistic
        t_statistic = correlation_coefficient * np.sqrt(df / (1 - correlation_coefficient**2))
        
        # Calculate the p-value
        p_value = 2 * (1 - stats.t.cdf(np.abs(t_statistic), df))

        return correlation_coefficient, p_value
    else:
        raise ValueError
    
def get_gene_pair_cor(gene1,gene2,method = "pearson|spearman", pseudobulk = None):
    if pseudobulk is None:
        raise ValueError
    try:
        df = pd.DataFrame({'gene1' : pseudobulk.loc[[gene1]].values.flatten().tolist(),'gene2' : pseudobulk.loc[[gene2]].values.flatten().tolist()})
    except:
        print(gene1,gene2)
    df = df[(df['gene1'] > 0.1) & (df['gene2'] > 0.1)]
    if len(df.index) < 30:
        cor_res = (0,0)
    else:
        cor_res =  fast_calculate_cor(df['gene1'].tolist(),df['gene2'].tolist(),method = method)
    return(cor_res)
    
def batch_calculate_cor(rep_num, GRN,pseudobulk = None,sample_ratio = None):
    if pseudobulk is None:
        raise ValueError
    #sample by 40% columns
    if sample_ratio is None:
        temp_tmp = pseudobulk
    else:
        temp_tmp = pseudobulk.sample(frac=sample_ratio,axis='columns')
    
    #pearson
    cor_pearson = 'cor_pearson_rep' + str(rep_num)
    pvalue_pearson = 'pvalue_pearson_rep' + str(rep_num)
    pvalue_pearson_adj = 'pvalue_pearson_adj_rep' + str(rep_num)

    pearson_res = []
    for i in range(0,len(GRN)):
        pearson_res.append(get_gene_pair_cor(GRN['TF'].tolist()[i],GRN['TargetGene'].tolist()[i],method = 'pearson',pseudobulk = temp_tmp))

    pearson_res = list(zip(*pearson_res))
    GRN.loc[:,cor_pearson] = pearson_res[0]
    GRN.loc[:,pvalue_pearson] = pearson_res[1]
    GRN.loc[:,pvalue_pearson_adj] = multipletests(pearson_res[1],method = 'fdr_bh')[1].tolist()

    #pearman
    cor_spearman = 'cor_spearman_rep' + str(rep_num)
    pvalue_spearman = 'pvalue_spearman_rep' + str(rep_num)
    pvalue_spearman_adj = 'pvalue_spearman_adj_rep' +str(rep_num)


    spearman_res = []
    for i in range(0,len(GRN)):
        spearman_res.append(get_gene_pair_cor(GRN['TF'].tolist()[i],GRN['TargetGene'].tolist()[i],method = 'spearman',pseudobulk = temp_tmp))

    spearman_res = list(zip(*spearman_res))
    GRN.loc[:,cor_spearman] = spearman_res[0]
    GRN.loc[:,pvalue_spearman] = spearman_res[1]
    GRN.loc[:,pvalue_spearman_adj] = multipletests(spearman_res[1],method = 'fdr_bh')[1].tolist()


    return(GRN)
## vectorization, more faster

def batch_calculate_cor2(rep_num, GRN,pseudobulk = None,sample_ratio = None):
    if pseudobulk is None:
        raise ValueError
    #sample by 40% columns
    if sample_ratio is None:
        temp_tmp = pseudobulk
    else:
        temp_tmp = pseudobulk.sample(frac=sample_ratio,axis='columns')
    
    #pearson
    cor_pearson = 'cor_pearson_rep' + str(rep_num)
    pvalue_pearson = 'pvalue_pearson_rep' + str(rep_num)
    pvalue_pearson_adj = 'pvalue_pearson_adj_rep' + str(rep_num)

    # pearson_res = []
    # for i in range(0,len(GRN)):
    #     pearson_res.append(get_gene_pair_cor(GRN['TF'].tolist()[i],GRN['TargetGene'].tolist()[i],method = 'pearson',pseudobulk = temp_tmp))
    pearson_res = list(map(partial(get_gene_pair_cor,method = 'pearson',pseudobulk = temp_tmp), GRN['TF'].tolist(), GRN['TargetGene'].tolist()))
    pearson_res = list(zip(*pearson_res))
    GRN.loc[:,cor_pearson] = pearson_res[0]
    GRN.loc[:,pvalue_pearson] = pearson_res[1]
    GRN.loc[:,pvalue_pearson_adj] = multipletests(pearson_res[1],method = 'fdr_bh')[1].tolist()

    #pearman
    cor_spearman = 'cor_spearman_rep' + str(rep_num)
    pvalue_spearman = 'pvalue_spearman_rep' + str(rep_num)
    pvalue_spearman_adj = 'pvalue_spearman_adj_rep' +str(rep_num)

    # spearman_res = []
    # for i in range(0,len(GRN)):
    #     spearman_res.append(get_gene_pair_cor(GRN['TF'].tolist()[i],GRN['TargetGene'].tolist()[i],method = 'spearman',pseudobulk = temp_tmp))
    spearman_res = list(map(partial(get_gene_pair_cor,method = 'spearman',pseudobulk = temp_tmp), GRN['TF'].tolist(), GRN['TargetGene'].tolist()))
    spearman_res = list(zip(*spearman_res))
    GRN.loc[:,cor_spearman] = spearman_res[0]
    GRN.loc[:,pvalue_spearman] = spearman_res[1]
    GRN.loc[:,pvalue_spearman_adj] = multipletests(spearman_res[1],method = 'fdr_bh')[1].tolist()
    return(GRN)


##load data
def usage():
    print('Usage: python script.py [input_tpm] [input_GRN] [out_GRN_cor]')




if __name__ == '__main__':
    try:      
        tpm = pd.read_csv(sys.argv[1],sep= '\t')
        tpm.index = tpm['gene']
        tpm = tpm.iloc[:,3:]
        grn = pd.read_csv(sys.argv[2],sep='\t')
        
        for i in tqdm(range(1,11)):
            print(f'Calculate rep{i}')
            grn = batch_calculate_cor2(i,grn,pseudobulk = tpm,sample_ratio= 0.4)
        grn.to_csv(sys.argv[3],sep='\t')
    except:
        usage()
