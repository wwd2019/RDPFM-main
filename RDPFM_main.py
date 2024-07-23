# 此文件完成RDPFM的模型训练部分流程。数据集预处理部分请观看planaria_main中前两步
import argparse
import os
from sklearn.metrics import accuracy_score, f1_score

from planaria_main.makePng import getMinSS_SMESG

os.environ["NCCL_DEBUG"] = "INFO"
os.environ["OMPI_MCA_opal_cuda_support"] = "true"
os.environ["CONDA_OVERRIDE_GLIBC"] = "2.56"
GPU_NUMBER = [6,7,8]
os.environ["CUDA_VISIBLE_DEVICES"] = ",".join([str(s) for s in GPU_NUMBER])
import pytz
import json
import pickle
from collections import Counter

import joblib
from sklearn import svm
from sklearn.model_selection import GridSearchCV, LeaveOneOut
import torch
import seaborn as sns; sns.set()
from transformers import BertForSequenceClassification
from transformers import Trainer
from geneformer import DataCollatorForCellClassification, EmbExtractor
import datetime
import subprocess
from datasets import load_from_disk
from transformers import BertConfig, BertForMaskedLM, TrainingArguments
from geneformer import GeneformerPretrainer
from tqdm import tqdm
import datasets
import random
import numpy as np
# step unsupervised_learning
def unsupervised_learning(unlabelled_data_path="./planaria_main/tk_dict/data_all_tk.datasets",tk_path="./planaria_main/tk_dict/gene_tk.pkl",p_num_path="./planaria_main/tk_dict/p_num.pkl",pretrained_model_path="/home/wwd/codebox/Geneformer/",saved_model_path='./planaria_main/model_train/model1_ep100'):
    # unlabelled_data_path是无标签数据集读取路径
    # tk_path 是基因的token_id字典的保存路径
    # pretrained_model_path是已预训练好的模型读取路径
    # saved_model_path是再训练后需要保存的模型保存路径
    # set local time/directories
    timezone = pytz.timezone("US/Eastern")
    rootdir = "/parent_ouput_directory"
    # set model parameters
    # model type
    model_type = "bert"
    # max input size
    max_input_size = 2 ** 11  # 2048
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
    num_examples = 615702 #应该改成对应数据集的细胞样本数
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
    epochs = 50
    # optimizer
    optimizer = "adamw"
    # weight_decay
    weight_decay = 0.001
    # output directories
    model_output_dir = saved_model_path
    training_output_dir = saved_model_path+'/train'
    logging_dir = saved_model_path+'/logs'
    # ensure not overwriting previously saved model
    model_output_file = os.path.join(model_output_dir, "pytorch_model.bin")
    if os.path.isfile(model_output_file) is True:
        raise Exception("Model already saved to this directory.")

    # make training and model output directories
    subprocess.call(f"mkdir {training_output_dir}", shell=True)
    subprocess.call(f"mkdir {model_output_dir}", shell=True)
    subprocess.call(f"mkdir {logging_dir}", shell=True)
    with open(tk_path, "rb") as fp:
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
    model = model.from_pretrained(pretrained_model_path)
    model = model.train()

    #
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
    # 您将看到有关未使用某些预训练权重以及某些权重被随机初始化的警告。别担心，这是完全正常的！ BERT 模型的预训练头被丢弃，并替换为随机初始化的分类头。您将在序列分类任务中微调这个新模型头，将预训练模型的知识传递给它
    print("Starting training.")
    # 定义trainer,加入真实细胞数据训练，主要更改此文件
    # define the new model
    with open(tk_path, "rb") as fp:
        token_dictionary = pickle.load(fp)
    trainer = GeneformerPretrainer(
        model=model,
        args=training_args,
        # pretraining corpus (e.g. https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main/genecorpus_30M_2048.dataset)
        train_dataset=load_from_disk(unlabelled_data_path),
        # file of lengths of each example cell (e.g. https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/blob/main/genecorpus_30M_2048_lengths.pkl)
        example_lengths_file="./planaria_main/tk_dict/p_num.pkl",
        token_dictionary=token_dictionary,
    )
    # train
    trainer.train()
    # save model
    trainer.save_model(model_output_dir)
    pass
# unsupervised_learning()
# step supervise_finetune
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
def supervise_finetune():
    # labels=['timepoint', 'unified_domain','timepoint_domain' ]
    labels = ['timepoint']
    labels2 = ['timepoint', 'unified_domain']
    gene2id = {}
    with open('./planaria_main/tk_dict/gene_tk.pkl', 'rb') as f:
        gene2id = pickle.load(f)
    for label in labels:
        # 使用有标签的数据进行微调
        with open('./planaria_main/data_related/labelled_p_cell_with_label.pkl', 'rb') as f1:
            data_all = pickle.load(f1)
        data_all2 = []
        for cell in tqdm(data_all):
            if cell[label] != '' and cell[label] != 'unknown':
                data_all2.append(cell)
        print(len(data_all2))
        random.shuffle(data_all2)
        p_dataset = datasets.Dataset.from_list(data_all2)
        p_dataset = p_dataset.rename_column(label, 'label')
        p_dataset = p_dataset.remove_columns(['indent'])  # 这个无意义
        p_dataset = p_dataset.remove_columns([item for item in labels2 if item != label])
        p_label_dict = {}
        i = 0
        for id in Counter(p_dataset['label']):
            p_label_dict[id] = i
            i += 1
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

        p_dataset = p_dataset.map(map_gene2id, num_proc=16)
        p_label_dict2 = Counter(p_dataset['label'])
        with open('./planaria_main/tk_dict/' + label + '_dict', 'wb') as f2:
            pickle.dump(p_label_dict, f2)
        with open('./planaria_main/tk_dict/' + label + '_num_dict', 'wb') as f2:
            pickle.dump(p_label_dict2, f2)
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
        # p_dataset.save_to_disk('./planaria_main/dataSets/' + label + '_data_aug_bet.datasets')
        p_dataset = datasets.Dataset.load_from_disk('./planaria_main/dataSets/' + label + '_data_aug_bet.datasets')
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
        output_dir = f"./planaria_main/results/model3_low/{datestamp}_geneformer_CellClassifier_L{max_input_size}_B{geneformer_batch_size}_LR{max_lr}_LS{lr_schedule_fn}_WU{warmup_steps}_E{epochs}_O{optimizer}_F{freeze_layers}/" + label + '/'
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
            # output_dir=f"/home/wwd/codebox/Geneformer/new/frog/model",
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
# supervise_finetune()
# step getEmbedding
def getEmbedding(dataset_path='./planaria_main/results/timepoint_data_aug.datasets',model_path='./planaria_main/results/timepoint_data_aug_labelled.datasets',save_path='./planaria_main/model_results/gene_ebeddings_aug2.pkl'):
    with open('./planaria_main/tk_dict/gene_tk.pkl', 'rb') as f:
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
        dataset_path+'2')
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
                             token_dictionary_file='./planaria_main/tk_dict/gene_tk.pkl')
        embs = embex.extract_embs(
            # "/home/wwd/codebox/Geneformer/new/planaria2/results/model3_aug_bet/240410_geneformer_CellClassifier_L2048_B16_LR1e-05_LSlinear_WU500_E38_Oadamw_F4/timepoint/checkpoint-187012",
            model_path,
            dataset_path+'2',
            "save_path",
            "output_prefix")

        with open(save_path+"/output_prefix.csv", 'r') as f:
            data_emb = f.readlines()
            for i in range(1, len(data_emb)):
                data_list = data_emb[i]
                data_list = data_list.strip('\n')
                data_list = data_list.split(',')
                if (data_list[0] not in gene_ebeddings):
                    gene_ebeddings[data_list[0]] = [float(item) for item in data_list[1:]]
    with open(save_path+'_embeddings.pkl', 'wb') as f:
        pickle.dump(gene_ebeddings, f)
    print(len(gene_ebeddings))
    pass
# get cell states
def getState(dataset_path='./planaria_main/dataSets/timepoint_data_aug_bet.datasets'):
    label_dict={0:'WT',1:'5d',3:'12h',5:'1.5d',6:'3d',4:'10d',2:'0h'}
    # 读取path路径下文件，判断每个细胞的再生阶段
    model=BertForSequenceClassification.from_pretrained('./planaria_main/model_train')
    # print(model)
    cell_Dataset=datasets.load_from_disk(dataset_path)
    with open('./planaria_main/tk_dict/gene_tk.pkl', 'rb') as f:
        gene2id = pickle.load(f)
    # print(gene2id)
    def map_gene2id(example):
        example['input_ids'] = [int(gene2id[item]) for item in example['input_ids'] if item in gene2id.keys()]
        random.shuffle(example['input_ids'])
        example['input_ids'] = example['input_ids'][:2024]
        return example
    cell_Dataset=cell_Dataset.map(map_gene2id,num_proc=16)
    # set training hyperparameters
    # max learning rate
    max_lr = 1e-5
    # batch size for training and eval
    geneformer_batch_size = 16
    # optimizer
    optimizer = "adamw"
    training_args = {
        "learning_rate": max_lr,
        "do_train": True,
        "do_eval": True,
        "evaluation_strategy": "epoch",
        "save_strategy": "epoch",
        "group_by_length": True,
        "length_column_name": "length",
        "disable_tqdm": False,
        "weight_decay": 0.001,
        "per_device_train_batch_size": geneformer_batch_size,
        "per_device_eval_batch_size": geneformer_batch_size,
        "load_best_model_at_end": True,
        'output_dir':'results'
    }
    training_args_init = TrainingArguments(**training_args)
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
    # create the trainer
    trainer = Trainer(
        model=model,
        args=training_args_init,
        data_collator=DataCollatorForCellClassification(),
        compute_metrics=compute_metrics
    )
    predict=trainer.predict(cell_Dataset)
    preds = predict.predictions.argmax(-1)
    labels=[label_dict[item] for item in preds]
    print(labels)
    return labels
    pass
# step gene_classification
# 已知的正负样本
negative=['SMESG000010270', 'SMESG000042630', 'SMESG000066574', 'SMESG000025574', 'SMESG000033607', 'SMESG000022746', 'SMESG000022746', 'SMESG000006969', 'SMESG000037117', 'SMESG000017894', 'SMESG000032691', 'SMESG000050058', 'SMESG000050058', 'SMESG000050059', 'SMESG000065245', 'SMESG000065247', 'SMESG000065253', 'SMESG000065254', 'SMESG000065254', 'SMESG000065269']
positive=['SMESG000005389', 'SMESG000046293', 'SMESG000063577', 'SMESG000020104', 'SMESG000060238', 'SMESG000029273', 'SMESG000053725', 'SMESG000070594', 'SMESG000070103', 'SMESG000049472', 'SMESG000016858', 'SMESG000060978', 'SMESG000037081', 'SMESG000016751', 'SMESG000056159', 'SMESG000069029', 'SMESG000045786', 'SMESG000070150', 'SMESG000064333', 'SMESG000025543', 'SMESG000048533', 'SMESG000060148', 'SMESG000079942', 'SMESG000020739', 'SMESG000042345', 'SMESG000057772', 'SMESG000033647', 'SMESG000048749', 'SMESG000036986', 'SMESG000060356', 'SMESG000030852', 'SMESG000044282', 'SMESG000025390', 'SMESG000061064', 'SMESG000074748', 'SMESG000078501', 'SMESG000042131', 'SMESG000051493', 'SMESG000049972', 'SMESG000074761', 'SMESG000072824']
def get_data():
    with open('./planaria_main/model_results/gene_ebeddings_aug2.pkl', 'rb') as f:
        gene_ebedding = pickle.load(f)
    # with open('./cluster_model/gene_marker_dict.pkl', 'rb') as f:
    #     gene_marker = pickle.load(f)
    positive_data=[] # 确定为正样本
    negative_data=[] # 确定为负样本
    unknown_data=[]  # 不确定
    for gene in negative:
        if(gene in gene_ebedding):
            negative_data.append({'gene_id':gene,'label':'negative','tensor':[float(item) for item in gene_ebedding[gene]]})
    for gene in positive:
        if (gene in gene_ebedding):
            positive_data.append({'gene_id': gene, 'label': 'positive', 'tensor': [float(item) for item in gene_ebedding[gene]]})
    for gene in gene_ebedding:
        if(gene not in positive and gene not in negative):
            unknown_data.append({'gene_id': gene, 'label': 'unknown', 'tensor': [float(item) for item in gene_ebedding[gene]]})
    # print(len(positive_data)) 34
    # print(len(negative_data)) 13
    # print(len(unknown_data)) 15117
    return positive_data, negative_data, unknown_data
def getMinSS_SMESG(ss,data_dict):
    negative_items=[]
    for data in data_dict:
        if(data_dict[data]['ss']<=ss):
            data=data.split('-')
            data=data[0]
            negative_items.append(data)
    return negative_items
def SVM_analyse(ss_max,id,embedding_path='./planaria_main/model_results/gene_ebeddings_aug2.pkl',ss_path='./planaria_main/dataSets/gene_data_allIndents_ss.pkl'):
    # 创建SVM分类器对象
    clf = svm.SVC()
    parameters = {'kernel': ('linear', 'rbf','poly'), 'C': [0.001, 0.01, 0.1, 1, 10, 100, 1000]}
    clf = GridSearchCV(clf, parameters)
    positive_data, negative_data, unknown_data=get_data()
    ebeddings=[]
    labels=[]
    loo = LeaveOneOut()
    for data in positive_data:
        ebeddings.append(data['tensor'])
        labels.append(1)
    for data in negative_data:
        ebeddings.append(data['tensor'])
        labels.append(0)
    with open(embedding_path, 'rb') as f:
        gene_ebedding = pickle.load(f)
    with open('./planaria_main/dataSets/gene_data_allIndents_ss.pkl', 'rb') as f:
        data_dict = pickle.load(f)
    # print(len(data_dict))
    negative_items=getMinSS_SMESG(ss_max,data_dict)
    negative_items_len=len(negative_items)
    print('new_negative',len(negative_items))
    for data in negative_items:
        if(data in gene_ebedding):
            ebeddings.append(gene_ebedding[data])
        labels.append(0)
    ebeddings=torch.tensor(ebeddings)
    labels=torch.tensor(labels)
    # 训练模型
    for train_index, test_index in loo.split(ebeddings):
        ebeddings_train, ebeddings_test = ebeddings[train_index], ebeddings[test_index]
        labels_train,labels_test = labels[train_index],labels[test_index]
        clf.fit(ebeddings_train, labels_train)
        # 测试模型
        true_num=0
        predictions = clf.predict(ebeddings_test)
        for i in range(len(predictions)):
            if predictions[i] == labels_test[i]:
                true_num += 1
    # print('true_num',true_num,'all',len(predictions))
    predictions = clf.predict(ebeddings)
    true_num=0
    for i in range(len(predictions)):
        if predictions[i]==labels[i]:
            true_num+=1
    print('true_num',true_num,'all',len(predictions),'accuracy',float(true_num)/len(predictions))
    joblib.dump(clf, './planaria_main/model_results/SVM_models/'+str(id)+'_0_svm_model.pkl')
    all_ebeddings=[]
    gene_ids=[]
    gene_result={}
    for data in unknown_data:
        all_ebeddings.append(data['tensor'])
        gene_ids.append(data['gene_id'])
    results=clf.predict(all_ebeddings)
    for i in tqdm(range(len(results))):
        result=results[i]
        if(result==1):
            gene_result[gene_ids[i]]='positive'
        else:
            gene_result[gene_ids[i]] = 'negative'
    print(Counter(results))
    with open('./planaria_main/model_results/SVM_models/'+str(id)+'_result_gene_classification.pkl', 'wb')as f:
        pickle.dump(gene_result,f)
    return Counter(results),true_num,len(predictions), negative_items_len
def SVM_20(start,lenth,target_path,embedding_path='./planaria_main/model_results/gene_ebeddings_aug2.pkl',ss_path='./planaria_main/dataSets/gene_data_allIndents_ss.pkl'):
    SVM_result={}
    for i in range(21):
        print('id',i)
        results,true_num,predictions,add_negtive=SVM_analyse(i*lenth+start,i,embedding_path,ss_path)
        results=dict(results)
        if(1 in results): # 存在正例
            SVM_result[i]={'positive':int(results[1]),'negative':int(results[0]),'add_negative_nums':add_negtive,'accuracy':float(true_num)/predictions,'ss':i*lenth+start,'id':str(i)}
            print('positive',int(results[1]))
        else:
            SVM_result[i] = {'positive': 0, 'negative': int(results[0]),'add_negative_nums':add_negtive,
                             'accuracy': float(true_num) / predictions,'ss':i*lenth+start, 'id': str(i)}
    with open(target_path, "w") as file:
        json.dump(SVM_result, file)
def gene_classification(low,high,model_savepath,embedding_path='./planaria_main/model_results/gene_ebeddings_aug2.pkl',ss_path='./planaria_main/dataSets/gene_data_allIndents_ss.pkl'):
    SVM_20(low, high, model_savepath)
    pass
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--act',
                        default='unsupervised_learning')
    args = parser.parse_args()
    act_name = args.act
    if(act_name=='unsupervised_learning'):
        unsupervised_learning()
    if(act_name=='supervise_finetune'):
        supervise_finetune()
    if (act_name == 'getStates'):
        getState()
    if(act_name=='getEmbedding'):
        getEmbedding()
    if(act_name=='gene_classification'):
        SVM_20(0.02, 0.0001, './planaria_main/results/svm_result.json')

# # unsupervised_learning
#     unsupervised_learning()
# # supervise_finetune()
#     supervise_finetune()
# # step getEmbedding
#     getEmbedding()
# # step gene_classification
#     SVM_20(0.02,0.0001,'./results/svm_result.json')
