import os
# GPU_NUMBER = [1,5]
# os.environ["CUDA_VISIBLE_DEVICES"] = ",".join([str(s) for s in GPU_NUMBER])
import random
import openpyxl
import pandas as pd
import datasets
from collections import Counter
from tqdm import tqdm
from sklearn.cluster import AffinityPropagation
from sklearn import metrics

from geneformer import EmbExtractor
import pickle
import numpy as np
# 此文件用于将训练好的模型中的基因 向量表示取出来，对其完成分类任务
# initiate EmbExtractor
# 获取所有基因的向量表达
def getEbeddings(dataset_path,dataset_target_path,target_path):
    with open('./tk_dict/gene_tk.pkl', 'rb') as f:
        gene_tk = pickle.load(f)
    gene_num = len(gene_tk)
    gene_ebeddings = {}
    i = 0
    j = 0
    data_all1 = []
    label = 'label'
    data_all2 = datasets.Dataset.load_from_disk(
        dataset_path)
    print('data_reload..')
    for item in tqdm(data_all2):
        item['emb_label'] = label + str(j)
        data_all1.append(item)
        i += 1
        if (i >= 3000):
            j += 1
            i = 0
    data_all1 = datasets.Dataset.from_list(data_all1)
    data_all1.save_to_disk(
        dataset_target_path)
    filter_dict = {}
    print(j)
    for k in tqdm(range(j + 1)):
        filter_dict['emb_label'] = 'label' + str(k)
        embex = EmbExtractor(model_type="Pretrained",
                             filter_data=filter_dict,
                             emb_mode="gene",
                             # emb_layer=-1,
                             emb_layer=-1,
                             max_ncells=3500,
                             forward_batch_size=12,
                             nproc=16,
                             token_dictionary_file='./tk_dict/gene_tk.pkl')
        embs = embex.extract_embs(
            # "/home/wwd/codebox/Geneformer/new/planaria2/results/model3_aug_bet/240410_geneformer_CellClassifier_L2048_B16_LR1e-05_LSlinear_WU500_E38_Oadamw_F4/timepoint/checkpoint-187012",
            '/home/wwd/codebox/Geneformer/new/planaria2/results/model3_aug_bet/240410_geneformer_CellClassifier_L2048_B16_LR1e-05_LSlinear_WU500_E38_Oadamw_F4/timepoint/checkpoint-247123',
            dataset_target_path,
            "./cluster_model/",
            "output_prefix")

        with open("./cluster_model/output_prefix.csv", 'r') as f:
            data_emb = f.readlines()
            for i in range(1, len(data_emb)):
                data_list = data_emb[i]
                data_list = data_list.strip('\n')
                data_list = data_list.split(',')
                if (data_list[0] not in gene_ebeddings):
                    gene_ebeddings[data_list[0]] =[float(item) for item in data_list[1:]]
    with open(target_path, 'wb') as f:
        pickle.dump(gene_ebeddings, f)
    print(len(gene_ebeddings))
getEbeddings('/home/wwd/codebox/Geneformer/new/planaria2/results/timepoint_data_aug.datasets','/home/wwd/codebox/Geneformer/new/planaria2/results/timepoint_data_aug_labelled.datasets','/home/wwd/codebox/Geneformer/new/planaria2/cluster_model/gene_ebeddings_aug2.pkl')




