# 此文件用于分析所有基因表达数据的基因表达量分布
import os,re
import pickle
from copy import copy,deepcopy
import pandas as pd
import numpy as np
import statistics
import time
from tqdm import trange, tqdm
from progressbar import ProgressBar, Percentage, Bar, Timer, ETA, FileTransferSpeed
widgets = ['Progress: ', Percentage(), ' ', Bar('#'), ' ', Timer(), ' ', ETA(), ' ', FileTransferSpeed()]
progress = ProgressBar(widgets=widgets)
import scanpy as sc
def data_pretreat():
    sc.set_figure_params(facecolor="white", figsize=(8, 8), dpi=100, color_map='viridis_r')
    sc.settings.verbosity = 3  # 设置日志等级: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    # print(os.getcwd())  # 查看当前路径
    # os.chdir('./filtered_gene_bc_matrices/scanpy') #修改路径
    outdir = './outdir'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # 导入 stereo-Seq 数据
    adata = sc.read_h5ad('./dataSets/spatial_unified_domain_5.h5ad')
    adata.var_names_make_unique()  # 索引去重，若上一步中使用 `var_names='gene_ids'` 则这一步非必须进行
    # 基础过滤：去除表达基因100以下的细胞；去除在3个细胞以下表达的基因。
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    mat = pd.DataFrame(data=adata.X.todense(), index=adata.obs_names.values, columns=adata.var_names.values)
    spatial = adata.obsm['spatial'].join(pd.DataFrame(adata.obs[['orig.ident', 'timepoint', 'unified_domain_5']]))
    cell_names = adata.obs_names.values
    datalists = mat.values
    data_dict = {}
    col = mat.columns
    print(cell_names)
    gene_data_all = {}  # 用于记录所有时间段基因的表达分布
    for i in tqdm(range(len(cell_names))):
        item = cell_names[i]
        timepoint = item[0:3]
        data_dict[item] = []
        for j in range(len(col)):
            if (datalists[i][j] != 0):
                data_dict[item].append(col[j] + ':' + str(datalists[i][j]))
                gene_id = col[j]
                gene_score = datalists[i][j]
                if (gene_id not in gene_data_all.keys()):
                    gene_data_all[gene_id] = {}
                    gene_data_all[gene_id][timepoint] = [gene_score]
                elif (timepoint not in gene_data_all[gene_id].keys()):
                    gene_data_all[gene_id][timepoint] = [gene_score]
                else:
                    gene_data_all[gene_id][timepoint].append(gene_score)
    with open('dataSets/gene_data_allIndents.pkl', 'wb') as f:
        pickle.dump(gene_data_all, f)
def getAverage(data_list):
    return sum(data_list)/len(data_list)
def getSS(data_list):
    mean = np.mean(data_list)
    # 计算离差平方和
    diff_sum = sum((x - mean) ** 2 for x in data_list)
    return diff_sum
def getSS_average():
    # 分析每个基因在每个再生阶段的SS值（离差平方和）
    with open('dataSets/gene_data_allIndents.pkl', 'rb') as f:
        adata = pickle.load(f)
    # gene_names = list(adata.keys())
    timepoint_dict = {'A1_': '0h', 'A2_': '0h', 'F1_': '10d', 'F2_': '10d', 'D1_': '3d', 'D2_': '3d', 'C1_': '1.5d',
                      'C2_': '1.5d', 'E1_': '5d', 'E2_': '5d', 'B1_': '12h', 'B2_': '12h', 'W_1': 'WT', 'W_2': 'WT'}
    gene_data_all = {}
    for gene in tqdm(adata):
        gene_data_all[gene] = {}
        for timepoint in adata[gene]:
            gene_data_all[gene][timepoint_dict[timepoint]] = getAverage(adata[gene][timepoint])
    print(gene_data_all)
    with open('dataSets/gene_data_allIndents_average.pkl', 'wb') as f:
        pickle.dump(gene_data_all, f)
    gene_act_all = []
    gene_act_all_Data = {}
    for item in gene_data_all:
        if (len(gene_data_all[item]) == 7):
            gene_act_all.append(item)
    for gene in gene_act_all:
        timepoint_all = []
        for timepoint in gene_data_all[gene]:
            timepoint_all.append(gene_data_all[gene][timepoint])
        gene_data_all[gene]['ss'] = getSS(timepoint_all)
        gene_act_all_Data[gene] = gene_data_all[gene]
    with open('dataSets/gene_data_allIndents_ss.pkl', 'wb') as f:
        pickle.dump(gene_act_all_Data, f)
    print(gene_act_all_Data)
    pass
# data_pretreat() #每个基因在每个再生阶段的表达量分布
# getSS_average() #每个基因在每个再生阶段的离差平方和分析