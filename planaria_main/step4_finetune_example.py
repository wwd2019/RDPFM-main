# 此文件用于涡虫标签数据的微调与测试
import pickle
import os,re
import pickle
from collections import Counter
GPU_NUMBER = [3]
os.environ["CUDA_VISIBLE_DEVICES"] = ",".join([str(s) for s in GPU_NUMBER])
import seaborn as sns; sns.set()
from datasets import load_from_disk
from sklearn.metrics import accuracy_score, f1_score
from transformers import BertForSequenceClassification
from transformers import Trainer
from transformers.training_args import TrainingArguments
from geneformer import DataCollatorForCellClassification
from tqdm import tqdm
import torch
import pandas as pd
from progressbar import ProgressBar, Percentage, Bar, Timer, ETA, FileTransferSpeed
import scanpy as sc
import loompy as lp
import datasets
import datetime
import random
import subprocess
import numpy as np
def compute_metrics(pred):
    labels = pred.label_ids
    preds = pred.predictions.argmax(-1)
    # calculate accuracy and macro f1 using sklearn's function
    acc = accuracy_score(labels, preds)
    # 评价
    macro_f1 = f1_score(labels, preds, average='macro')
    return {
      'accuracy': acc,
      'macro_f1': macro_f1
    }
def makedata_labelled(input_path,output_path='./data_treated_sets/labelled_p_cell_with_label.pkl'):
    data_all=[]
    widgets = ['Progress: ', Percentage(), ' ', Bar('#'), ' ', Timer(), ' ', ETA(), ' ', FileTransferSpeed()]
    progress = ProgressBar(widgets=widgets)
    sc.set_figure_params(facecolor="white", figsize=(8, 8), dpi=100, color_map='viridis_r')
    sc.settings.verbosity = 3  # 设置日志等级: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    print(os.getcwd())  # 查看当前路径
    # os.chdir('./filtered_gene_bc_matrices/scanpy') #修改路径
    datdir = input_path
    # 导入 stereo-Seq 数据
    adata = sc.read_h5ad(datdir)
    adata.var_names_make_unique()  # 索引去重，若上一步中使用 `var_names='gene_ids'` 则这一步非必须进行
    # 用于存储分析结果文件的路径
    results_file = output_path
    # 基础过滤：去除表达基因200以下的细胞；去除在3个细胞以下表达的基因。
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    mat = pd.DataFrame(data=adata.X.todense(), index=adata.obs_names.values, columns=adata.var_names.values)
    # mat.to_csv(os.path.join(outdir, 'mat.csv'), index=True)
    spatial = adata.obsm['spatial'].join(pd.DataFrame(adata.obs[['orig.ident', 'timepoint', 'unified_domain_5']]))
    rep = ['1', '3', '5']
    reps = '|'.join(rep)
    spatial_rep1 = spatial[spatial['orig.ident'].str.contains(reps)]
    mat_rep1 = mat.loc[spatial['orig.ident'].str.contains(reps)]
    for slice_ in np.unique(spatial_rep1['timepoint']):
        print(slice_)
        spatial_rep1_slice = spatial_rep1.loc[spatial_rep1['timepoint'] == slice_, :].sort_index()
        mat_rep1_slice = mat_rep1.loc[spatial_rep1_slice.index, :].sort_index()
        print(mat_rep1_slice.shape)
    cell_names = adata.obs_names.values
    print(spatial.shape)
    print(len(cell_names))
    datalists = mat.values
    data_dict = {}
    col = mat.columns
    for i in tqdm(range(len(cell_names))):
        item = cell_names[i]
        data_dict[item] = []
        for j in range(len(col)):
            if (datalists[i][j] != 0):
                data_dict[item].append(col[j] + ':' + str(datalists[i][j]))
    i=0
    print(len(data_dict))
    j=0
    for cell in data_dict.keys():
        if(spatial['orig.ident'][i] and spatial['timepoint'][i] and spatial['unified_domain_5'][i]):
            data_all.append({'cell_genes': data_dict[cell],'length': len(data_dict[cell]), 'indent':spatial['orig.ident'][i],'timepoint':spatial['timepoint'][i],'unified_domain':spatial['unified_domain_5'][i]})
        else:
            data_all.append(
                {'cell_genes': data_dict[cell], 'length': len(data_dict[cell]), 'indent': 'unknown','timepoint':'unknown', 'unified_domain': 'unknown'})
            j+=1
        i+=1
    print(j)
    with open(output_path, 'wb') as f1:
        pickle.dump(data_all,f1)
def p_cell_finetune(dataset_path):
    # labels=['timepoint', 'unified_domain','timepoint_domain' ]
    labels = ['timepoint']
    labels2=['timepoint', 'unified_domain']
    gene2id={}
    with open('./tk_dict/gene_tk.pkl', 'rb') as f:
        gene2id = pickle.load(f)
    for label in labels:
        # 使用有标签的数据进行微调
        with open(dataset_path, 'rb') as f1:
            data_all = pickle.load(f1)
        data_all2=[]
        for cell in tqdm(data_all):
            if cell[label] != '' and cell[label] != 'unknown':
                data_all2.append(cell)
        print(len(data_all2))
        random.shuffle(data_all2)
        p_dataset = datasets.Dataset.from_list(data_all2)
        p_dataset = p_dataset.rename_column(label,'label')
        p_dataset = p_dataset.remove_columns(['indent']) # 这个无意义
        p_dataset = p_dataset.remove_columns([item for item in labels2 if item !=label])
        p_label_dict = {}
        i=0
        for id in Counter(p_dataset['label']):
            p_label_dict[id]=i
            i+=1
        print(p_label_dict)
        print(Counter(p_dataset['label']))
        def map_gene2id(example):
            example['input_ids'] = []
            for item in example['cell_genes']:
                item = item.split(':')
                item = item[0]
                if (item in gene2id):
                    example['input_ids'].append(int(gene2id[item]))
                example['input_ids'] = example['input_ids'][:2024]
            example['label'] = p_label_dict[example['label']]
            return example
        p_dataset = p_dataset.map(map_gene2id,num_proc=16)

        p_label_dict2=Counter(p_dataset['label'])
        with open('./tk_dict/'+label+'_dict', 'wb') as f2:
            pickle.dump(p_label_dict,f2)
        with open('./tk_dict/'+label+'_num_dict', 'wb') as f2:
            pickle.dump(p_label_dict2,f2)

        # # 使用有标签的数据进行微调
        # with open(dataset_path, 'rb') as f1:
        #      data_all = pickle.load(f1)
        # data_all2 = []
        # for cell in tqdm(data_all):
        #     if len(cell['cell_genes']) >= 2024:
        #         data_all_item=[]
        #         for i in range(15*(len(cell['cell_genes'])-2024)):
        #             random.shuffle(cell['cell_genes'])
        #             # 为防止数据泄露，去重
        #             data_all_item.append(cell['cell_genes'][:2024])
        #         data_all_item= [list(t) for t in set(tuple(sublist) for sublist in data_all_item)]
        #         for item in data_all_item:
        #             data_all2.append({'cell_genes':item,'timepoint':cell['timepoint'],'length':2024})
        #     else:
        #         data_all2.append(cell)
        # print(len(data_all2))
        # random.shuffle(data_all2)
        # p_dataset = datasets.Dataset.from_list(data_all2)
        # p_dataset = p_dataset.rename_column(label, 'label')
        # p_dataset = p_dataset.remove_columns(['indent'])  # 这个无意义
        # p_dataset = p_dataset.remove_columns([item for item in labels2 if item != label])
        # p_label_dict = {}
        # i = 0
        # for id in Counter(p_dataset['label']):
        #     p_label_dict[id] = i
        #     i += 1
        # print(p_label_dict)
        # print(Counter(p_dataset['label']))
        # def map_gene2id(example):
        #     example['input_ids'] = []
        #     for item in example['cell_genes']:
        #         item = item.split(':')
        #         item = item[0]
        #         if (item in gene2id):
        #             example['input_ids'].append(int(gene2id[item]))
        #         random.shuffle(example['input_ids'])
        #         example['input_ids'] = example['input_ids'][:2024]
        #     example['label'] = p_label_dict[example['label']]
        #     return example
        # p_dataset = p_dataset.map(map_gene2id, num_proc=16)
        p_dataset = datasets.Dataset.load_from_disk('./results/' + label + '_data_aug_bet.datasets')
        p_trainset = p_dataset.select([i for i in range(0, round(len(p_dataset) * 0.8))])
        p_evalset = p_dataset.select(
            [i for i in range(round(len(p_dataset) * 0.8), len(p_dataset))])
        print('fine_tune start')
        # set model parameters
        # max input size
        max_input_size = 2 ** 11  # 2048
        # set training hyperparameters
        # max learning rate
        max_lr = 1e-5
        # how many pretrained layers to freeze
        freeze_layers = 4
        # number gpus
        num_gpus = 2
        # number cpu cores
        num_proc = 6
        # batch size for training and eval
        geneformer_batch_size = 16
        # learning schedule
        lr_schedule_fn = "linear"
        # warmup steps
        warmup_steps = 500
        # number of epochs
        epochs = 40
        # optimizer
        optimizer = "adamw"
        # set logging steps
        logging_steps = round(len(p_trainset) / geneformer_batch_size / 10)
        print(p_label_dict.keys())
        # os.environ["CUDA_VISIBLE_DEVICES"] = '8'
        # reload pretrained model
        model = BertForSequenceClassification.from_pretrained("./model_train/model1_ep100/train/checkpoint-346320",
                                                              num_labels=len(p_label_dict.keys()),
                                                              output_attentions=False,
                                                              output_hidden_states=False).to('cuda')
        # model.load_state_dict(torch.load('./model_train/model2_ep30/model.safetensors'))
        # define output directory path
        current_date = datetime.datetime.now()
        datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"
        output_dir = f"./results/model3_low/{datestamp}_geneformer_CellClassifier_L{max_input_size}_B{geneformer_batch_size}_LR{max_lr}_LS{lr_schedule_fn}_WU{warmup_steps}_E{epochs}_O{optimizer}_F{freeze_layers}/"+label+'/'
        os.makedirs(output_dir, exist_ok=True)
        # ensure not overwriting previously saved model
        saved_model_test = os.path.join(output_dir, f"pytorch_model.bin")
        if os.path.isfile(saved_model_test) == True:
            raise Exception("Model already saved to this directory.")
        # make output directory
        subprocess.call(f'mkdir {output_dir}', shell=True)
        # set training arguments
        training_args = {
            "learning_rate": max_lr,
            "do_train": True,
            "do_eval": True,
            "evaluation_strategy": "epoch",
            "save_strategy": "epoch",
            "logging_steps": logging_steps,
            "group_by_length": True,
            "length_column_name": "length",
            "disable_tqdm": False,
            "lr_scheduler_type": lr_schedule_fn,
            "warmup_steps": warmup_steps,
            "weight_decay": 0.001,
            "per_device_train_batch_size": geneformer_batch_size,
            "per_device_eval_batch_size": geneformer_batch_size,
            "num_train_epochs": epochs,
            "load_best_model_at_end": True,
            "output_dir": output_dir,

        }
        training_args_init = TrainingArguments(**training_args)
        # create the trainer
        trainer = Trainer(
            model=model,
            args=training_args_init,
            data_collator=DataCollatorForCellClassification(),
            train_dataset=p_trainset,
            eval_dataset=p_evalset,
            compute_metrics=compute_metrics
        )
        # train the cell type classifier
        trainer.train()
        trainer.save_model()
    pass
# makedata_labelled('./data_treated_sets/spatial_unified_domain_5.h5ad')
p_cell_finetune('./data_treated_sets/labelled_p_cell_with_label.pkl')