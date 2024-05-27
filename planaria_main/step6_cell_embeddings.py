import os

import umap
from matplotlib.colors import LinearSegmentedColormap
GPU_NUMBER = [1,2,3,4,5]
os.environ["CUDA_VISIBLE_DEVICES"] = ",".join([str(s) for s in GPU_NUMBER])
import matplotlib.pyplot as plt
import datasets
from collections import Counter
from tqdm import tqdm
from geneformer import EmbExtractor
import pickle
import numpy as np
def getEbeddings(dataset_path,dataset_target_path,target_path):
    #  获取所有标签细胞的embeddings，使用各位平均值与标签作图
    with open('./tk_dict/gene_tk.pkl', 'rb') as f:
        gene_tk = pickle.load(f)
    cell_ebeddings = {}
    labels=[]
    data_all1 = []
    label = 'label'
    data_all2 = datasets.Dataset.load_from_disk(
        dataset_path)
    print('data_reload..')
    for item in tqdm(data_all2):
        item['emb_label'] = label + str(item['label'])
        if(label + str(item['label']) not in labels):
            labels.append(label + str(item['label']))
        data_all1.append(item)
    print(Counter(data_all2['label']))
    data_all1 = datasets.Dataset.from_list(data_all1)
    data_all1.save_to_disk(
        dataset_target_path)
    filter_dict = {}
    for k in tqdm(labels):
        filter_dict['emb_label'] = k
        embex = EmbExtractor(model_type="Pretrained",
                             filter_data=filter_dict,
                             emb_mode="cell",
                             # emb_layer=-1,
                             emb_layer=0,
                             max_ncells=50000,
                             forward_batch_size=12,
                             nproc=16,
                             token_dictionary_file='./tk_dict/gene_tk.pkl')
        embs = embex.extract_embs(
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
                if(k not in cell_ebeddings):
                    cell_ebeddings[k]=[[float(item) for item in data_list[1:]]]
                else:
                    cell_ebeddings[k].append([float(item) for item in data_list[1:]])

    with open(target_path, 'wb') as f:
        pickle.dump(cell_ebeddings, f)
    print(len(cell_ebeddings))
    print(cell_ebeddings.keys())
# getEbeddings('/home/wwd/codebox/Geneformer/new/planaria2/results/timepoint_data_aug.datasets',
#              '/home/wwd/codebox/Geneformer/new/planaria2/results/timepoint_data_aug_labelled.datasets',
#              '/home/wwd/codebox/Geneformer/new/planaria2/cluster_model/cell_ebeddings.pkl')
