# 此部分用于实现无监督学习
import os
os.environ["NCCL_DEBUG"] = "INFO"
os.environ["OMPI_MCA_opal_cuda_support"] = "true"
os.environ["CONDA_OVERRIDE_GLIBC"] = "2.56"
GPU_NUMBER = [6,7,8]
os.environ["CUDA_VISIBLE_DEVICES"] = ",".join([str(s) for s in GPU_NUMBER])
import pickle
import subprocess
import pytz
import torch
from datasets import load_from_disk
from transformers import BertConfig, BertForMaskedLM, TrainingArguments
from geneformer import GeneformerPretrainer
import os
import sys
from loompy import loompy
from tqdm import tqdm
import math
import pandas as pd
from geneformer import TranscriptomeTokenizer
from progressbar import ProgressBar, Percentage, Bar, Timer, ETA, FileTransferSpeed
import scanpy as sc
import datasets
import random
import numpy as np
# 首先根据geneformer要求，将基因名称化为秩值编码
def view_bar(message, num, total):
    # 进度条
    rate = num / total
    rate_num = int(rate * 30)
    rate_nums = math.ceil(rate * 100)
    r = '\r%s:[%s%s]%d%%\t%d/%d' % (message, "|" * rate_num, " " * (30 - rate_num), rate_nums, num, total,)
    sys.stdout.write(r)
    sys.stdout.flush()
def getDict():
    with open('./dataSets/smesg.transcripts.fa','r')as f1:
        datas=f1.readlines()
    genenames=[]
    genenames2=[]
    genedict={}
    for data in datas:
        if('SMEST'in data):
            data=data.split(' ')
            smest_id=data[0]
            smesg_id=data[1]
            smest_id=smest_id[1:]
            smest_id=smest_id.split('.')
            smest_id=smest_id[0]
            smesg_id=smesg_id.split('=')
            smesg_id=smesg_id[1]
            smesg_id=smesg_id.split('.')
            smesg_id=smesg_id[0]
            genenames.append(smesg_id)
            genenames2.append(smest_id)
            genedict[smest_id]=smesg_id
    return genenames,genenames2,genedict
def  data_pretreat_tokenizer():
    # 根据基因表达的频率排序基因id
    data_all = datasets.Dataset.load_from_disk("./data_treated_sets/other_valid_p_cell/p_all_unlabell_cell.dataset")
    gene_names, gene_names2, geneDict = getDict()

    gene_nums=dict(zip(gene_names,[0 for i in range(len(gene_names))]))
    for cell in tqdm(data_all):
        cell_genes = cell['cell_genes']
        for gene in gene_names:
            if (gene in cell_genes):
                gene_nums[gene] += 1
    for item in gene_nums:
        gene_nums[item]=float(gene_nums[item])/len(data_all)
    with open('./tk_dict/gene_dictionary.pkl','wb')as f:
        pickle.dump(gene_nums,f)
    gene_nums=sorted(gene_nums.items(), reverse=True)
    gene_tk= {}
    i=1
    for item in gene_nums:
        gene_tk[item[0]]=i
        i+=1
    with open('./tk_dict/gene_tk.pkl','wb') as f:
        pickle.dump( gene_tk,f)

    tk = TranscriptomeTokenizer({"cell_type": "cell_type"}, nproc=16,gene_median_file='./tk_dict/gene_dictionary.pkl',
        token_dictionary_file='./tk_dict/gene_tk.pkl')
    tk.tokenize_data("./loom_data",
                     "./tk_dict",
                     "p_cell_all")
    pass
def data_pretreat(data_path):
    geneDict = {}
    # 读取fa文件
    gene_names,gene_names2,geneDict=getDict()
    # # data_path 为经过处理的所有涡虫无标签数据路径
    data_all = []
    for item in data_path:
        with open(item, 'rb') as f1:
            data_all_item = pickle.load(f1)
        data_all += data_all_item
    data_all_dataset=datasets.Dataset.from_list(data_all)
    def map_id(example):
        example['cell_genes']=[geneDict[item] for item in example['cell_genes'] if 'SMEST' in item]
        example['length']=len(example['cell_genes'])
        return example
        pass
    data_all=data_all_dataset.map(map_id,num_proc=16)
    data_all.save_to_disk("./data_treated_sets/other_valid_p_cell/p_all_unlabell_cell.dataset")
    data_all=datasets.Dataset.load_from_disk("./data_treated_sets/other_valid_p_cell/p_all_unlabell_cell.dataset")
    print(len(data_all))
    # 去重
    gene_names=list(set(gene_names))
    # 获取表达矩阵
    p_cell_matrix=[]
    p_cell_nums=[]
    p_cell_type=[]
    with open("./data_treated_sets/gene_matrix.txt",'r')as f:
        data_lists=f.readlines()
        for data_list in tqdm(data_lists):
            data_list=data_list.strip('\n')
            data_list=data_list.split(' ')
            p_cell_matrix_item=[int(item) for item in data_list if item=='1' or item=='0']
            p_cell_nums.append(len(p_cell_matrix_item))
            p_cell_type.append('example_cell')
            p_cell_matrix.append(p_cell_matrix_item)

    # for cell in tqdm(data_all):
    #     cell_genes= cell['cell_genes']
    #     p_cell_nums.append(len(cell_genes))
    #     p_cell_type.append('unlabelled')
    #     p_cell_matrix_item=[]
    #     for gene in gene_names:
    #         if(gene in cell_genes):
    #             p_cell_matrix_item.append(1)
    #             gene_nums[gene_names.index(gene)]+=1
    #         else:
    #             p_cell_matrix_item.append(0)
    #     p_cell_matrix.append(p_cell_matrix_item)
    # with open("./data_treated_sets/gene_matrix.txt",'w')as f:
    #     for item in p_cell_matrix:
    #         for item0 in item:
    #             f.write(str(item0)+' ')
    #         f.write('\n')
    print('step success')
    # p_cell_matrix=np.array(p_cell_matrix)
    # p_cell_nums=np.array(p_cell_nums)
    # p_cell_type=np.array(p_cell_type)
    row_attrs = {
        "ensembl_id": gene_names,  # 基因的id
        # "gene_nums": gene_nums
        # "Gene": np.array(adata.var_names) ,
    }
    col_attrs = {
        "n_counts": p_cell_nums,
        "cell_type":p_cell_type
    }
    print(len(p_cell_matrix))

    p_cell_matrix=np.array(p_cell_matrix)
    print('step2 success')
    p_cell_matrix=p_cell_matrix.transpose()
    with open('./tk_dict/p_num.pkl','wb')as f1:
        pickle.dump(p_cell_nums,f1)
        print('nums saved')
    # 创建一个新的.loom文件
    loom_file = "./loom_data/p_cell.loom"
    loompy.create(loom_file,p_cell_matrix,row_attrs=row_attrs,col_attrs=col_attrs)
    print(loom_file+' saved successfully')
    # lp.create("./loom_data/p_cell.loom", p_cell_matrix.transpose(), row_attrs=row_attrs, col_attrs=col_attrs)
    data_pretreat_tokenizer()
    print('transfered success')
    pass
def data_get_us(input_path,output_path='./data_treated_sets/labelled_p_cell.pkl'):
    # 将已有的数据也加入涡虫总数据中
    outdir='./outdir'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
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
    # 基础过滤：去除表达基因10以下的细胞；去除在3个细胞以下表达的基因。
    sc.pp.filter_cells(adata, min_genes=20)
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
        spatial_rep1_slice.to_csv(os.path.join(outdir, 'spatial_rep1-slice_%s.csv' % slice_), index=True)
        mat_rep1_slice = mat_rep1.loc[spatial_rep1_slice.index, :].sort_index()
        mat_rep1_slice.to_csv(os.path.join(outdir, 'mat_rep1-slice_%s.csv' % slice_), index=True)
        print(mat_rep1_slice.shape)
    cell_names = adata.obs_names.values
    datalists = mat.values
    data_dict = {}
    col = mat.columns
    for i in tqdm(range(len(cell_names))):
        item = cell_names[i]
        data_dict[item] = []
        for j in range(len(col)):
            if (datalists[i][j] != 0):
                data_dict[item].append(col[j] + ':' + str(datalists[i][j]))
    for cell in data_dict.keys():
        data_all.append({'cell_genes': data_dict[cell],'length': len(data_dict[cell])})
    with open(output_path, 'wb') as f1:
        pickle.dump(data_all,f1)
    print(len(data_all))
    pass
# data_get_us('./data_treated_sets/spatial_unified_domain_5.h5ad')
# data_list=['./data_treated_sets/other_valid_p_cell/p_cell_unannotated.pkl','./data_treated_sets/labelled_p_cell.pkl']
# data_pretreat(data_list)
# getDict()
# data_all = datasets.Dataset.load_from_disk("./data_treated_sets/other_valid_p_cell/p_all_unlabell_cell.dataset")
data_all = datasets.Dataset.load_from_disk("./dataSets/p_all_unlabell_cell.dataset")
id2tk={}
with open('./tk_dict/gene_tk.pkl', 'rb') as f:
    id2tk=pickle.load(f)
data_all2=[]
for example in tqdm(data_all):
    # 上限2048
    input_ids=[id2tk[item] for item in example['cell_genes'] if id2tk[item]]
    random.shuffle(input_ids)
    example['input_ids'] = input_ids[:2048]
    example['length']=len(example['cell_genes'])
    if(example['length']>=20):
        data_all2.append(example)
random.shuffle(data_all2)
print(len(data_all2))
data_all2=datasets.Dataset.from_list(data_all2)
data_all2.save_to_disk('./tk_dict/data_all_tk2.datasets')

seed_num = 0
random.seed(seed_num)
np.random.seed(seed_num)
seed_val = 42
torch.manual_seed(seed_val)
torch.cuda.manual_seed_all(seed_val)

# set local time/directories
timezone = pytz.timezone("US/Eastern")
rootdir = "/parent_ouput_directory"

# set model parameters
# model type
model_type = "bert"
# max input size
max_input_size = 2**11  # 2048
# number of layers
num_layers = 6
# number of attention heads
num_attn_heads = 4
# number of embedding dimensions
num_embed_dim = 256
# intermediate size
intermed_size = num_embed_dim * 2
# activation function
activ_fn = "relu"
# initializer range, layer norm, dropout
initializer_range = 0.02
layer_norm_eps = 1e-12
attention_probs_dropout_prob = 0.02
hidden_dropout_prob = 0.018


# set training parameters
# total number of examples in Genecorpus-30M after QC filtering:
num_examples = 615702
# number gpus
num_gpus = 3
# batch size for training and eval
geneformer_batch_size = 8
# max learning rate
max_lr = 1e-5
# learning schedule
lr_schedule_fn = "linear"
# warmup steps
warmup_steps = 10_000
# number of epochs
# epochs = 50
epochs = 30
# optimizer
optimizer = "adamw"
# weight_decay
weight_decay = 0.001


# output directories
model_output_dir='./model_train/model1_ep100'
training_output_dir='./model_train/model1_ep100/train'
logging_dir='./model_train/model1_ep100/logs'

# ensure not overwriting previously saved model
model_output_file = os.path.join(model_output_dir, "pytorch_model.bin")
if os.path.isfile(model_output_file) is True:
    raise Exception("Model already saved to this directory.")


# make training and model output directories
subprocess.call(f"mkdir {training_output_dir}", shell=True)
subprocess.call(f"mkdir {model_output_dir}", shell=True)
subprocess.call(f"mkdir {logging_dir}", shell=True)
# 这里需要更改
# load gene_ensembl_id:token dictionary (e.g. https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/blob/main/token_dictionary.pkl)
# with open("/home/wwd/codebox/Geneformer/geneformer/token_dictionary.pkl", "rb") as fp:
#     token_dictionary = pickle.load(fp)
with open("./tk_dict/gene_tk.pkl", "rb") as fp:
    token_dictionary = pickle.load(fp)
# model configuration
config = {
    "hidden_size": num_embed_dim,
    "num_hidden_layers": num_layers,
    "initializer_range": initializer_range,
    "layer_norm_eps": layer_norm_eps,
    "attention_probs_dropout_prob": attention_probs_dropout_prob,
    "hidden_dropout_prob": hidden_dropout_prob,
    "intermediate_size": intermed_size,
    "hidden_act": activ_fn,
    "max_position_embeddings": max_input_size,
    "model_type": model_type,
    "num_attention_heads": num_attn_heads,
    "pad_token_id": token_dictionary.get("<pad>"),
    "vocab_size": len(token_dictionary),  # genes+2 for <mask> and <pad> tokens
}
# 设置基础的bert
config = BertConfig(**config)
model = BertForMaskedLM(config)
# 使用预训练好的模型
model=model.from_pretrained("./model_train/model1_ep100/train/checkpoint-923520")
model = model.train()

#
'''
分词器返回一个包含三个重要项目的字典：

input_ids是句子中每个标记对应的索引。
Attention_mask指示是否应注意令牌。
token_type_ids标识当存在多个序列时令牌属于哪一个序列。
通过解码返回您的输入input_ids：
'''
# define the training arguments
training_args = {
    "learning_rate": max_lr,
    "do_train": True,
    "do_eval": False,
    "group_by_length": True,
    "length_column_name": "length",
    "disable_tqdm": False,
    "lr_scheduler_type": lr_schedule_fn,
    "warmup_steps": warmup_steps,
    "weight_decay": weight_decay,
    "per_device_train_batch_size": geneformer_batch_size,
    "num_train_epochs": epochs,
    "save_strategy": "steps",
    "save_steps": np.floor(num_examples / geneformer_batch_size / 8),  # 8 saves per epoch
    "logging_steps": 1000,
    "output_dir": training_output_dir,
    "logging_dir": logging_dir,
}
training_args = TrainingArguments(**training_args)
#您将看到有关未使用某些预训练权重以及某些权重被随机初始化的警告。别担心，这是完全正常的！ BERT 模型的预训练头被丢弃，并替换为随机初始化的分类头。您将在序列分类任务中微调这个新模型头，将预训练模型的知识传递给它
print("Starting training.")
# 定义trainer,加入真实细胞数据训练，主要更改此文件
# define the new model
with open("./tk_dict/gene_tk.pkl", "rb") as fp:
    token_dictionary = pickle.load(fp)
trainer = GeneformerPretrainer(
    model=model,
    args=training_args,
    # pretraining corpus (e.g. https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main/genecorpus_30M_2048.dataset)
    train_dataset=load_from_disk('./tk_dict/data_all_tk.datasets'),
    # file of lengths of each example cell (e.g. https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/blob/main/genecorpus_30M_2048_lengths.pkl)
    example_lengths_file="./tk_dict/p_num.pkl",
    token_dictionary=token_dictionary,
)
# train
trainer.train()
# save model
trainer.save_model(model_output_dir)