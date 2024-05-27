import json
import pickle
import random
from collections import Counter

import joblib
import numpy as np
import openpyxl
from sklearn import svm
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from torch import nn
from d2l import torch as d2l
import torch
from torch.nn import functional as F
from tqdm import tqdm
# negative=['SMED30002212','SMED30006871','SMED30007365','SMED30024071','SMED30020620','SMED30001259', 'SMED30034753','SMED30004716','SMED30022815','SMED30004877','SMED30019948']
# 已知的正负样本
negative=['SMESG000010270', 'SMESG000042630', 'SMESG000066574', 'SMESG000025574', 'SMESG000033607', 'SMESG000022746', 'SMESG000022746', 'SMESG000006969', 'SMESG000037117', 'SMESG000017894', 'SMESG000032691', 'SMESG000050058', 'SMESG000050058', 'SMESG000050059', 'SMESG000065245', 'SMESG000065247', 'SMESG000065253', 'SMESG000065254', 'SMESG000065254', 'SMESG000065269']
positive=['SMESG000005389', 'SMESG000046293', 'SMESG000063577', 'SMESG000020104', 'SMESG000060238', 'SMESG000029273', 'SMESG000053725', 'SMESG000070594', 'SMESG000070103', 'SMESG000049472', 'SMESG000016858', 'SMESG000060978', 'SMESG000037081', 'SMESG000016751', 'SMESG000056159', 'SMESG000069029', 'SMESG000045786', 'SMESG000070150', 'SMESG000064333', 'SMESG000025543', 'SMESG000048533', 'SMESG000060148', 'SMESG000079942', 'SMESG000020739', 'SMESG000042345', 'SMESG000057772', 'SMESG000033647', 'SMESG000048749', 'SMESG000036986', 'SMESG000060356', 'SMESG000030852', 'SMESG000044282', 'SMESG000025390', 'SMESG000061064', 'SMESG000074748', 'SMESG000078501', 'SMESG000042131', 'SMESG000051493', 'SMESG000049972', 'SMESG000074761', 'SMESG000072824']
def get_data():
    with open('./model_results/gene_ebeddings_aug2.pkl', 'rb') as f:
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
def getSS_SMESG(data_dict,neg_lenth=10,pos_lenth=10):
    negative_items=[]
    positive_items=[]
    data_dict_list=[]
    for data in data_dict:
        data_dict_list.append({'ss':data_dict[data]['ss'],'gene_id':data})
    data_list=sorted(data_dict_list,key=lambda x:x['ss'])
    for data in data_list[:neg_lenth]:
        data=data['gene_id']
        data = data.split('-')
        data=data[0]
        negative_items.append(data)
    for data in data_list[-1-pos_lenth:-1]:
        data = data['gene_id']
        data = data.split('-')
        data = data[0]
        positive_items.append(data)
    # for data in data_dict:
    #     if(data_dict[data]['ss']<=ss):
    #         data=data.split('-')
    #         data=data[0]
    #         negative_items.append(data)
    return negative_items,positive_items

def SVM_analyse(ss_max,id,embedding_path='/model_results/gene_ebeddings_aug2.pkl',ss_path='dataSets/gene_data_allIndents_ss.pkl'):
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
    with open('/model_results/gene_ebeddings_aug2.pkl', 'rb') as f:
        gene_ebedding = pickle.load(f)
    with open('dataSets/gene_data_allIndents_ss.pkl', 'rb') as f:
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
    joblib.dump(clf, './model_results/SVM_models/'+str(id)+'_0_svm_model.pkl')
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
    with open('SVM_models/'+str(id)+'_result_gene_classification.pkl', 'wb')as f:
        pickle.dump(gene_result,f)
    return Counter(results),true_num,len(predictions), negative_items_len
def SVM_20(start,lenth,target_path,embedding_path='/model_results/gene_ebeddings_aug2.pkl',ss_path='./dataSets/gene_data_allIndents_ss.pkl'):
    SVM_result={}
    for i in range(21):
        print('id',i)
        results,true_num,predictions,add_negtive=SVM_analyse(i*lenth+start,i,embedding_path,ss_path)
        results=dict(results)
        if(1 in results):
            SVM_result[i]={'positive':int(results[1]),'negative':int(results[0]),'add_negative_nums':add_negtive,'accuracy':float(true_num)/predictions,'ss':i*lenth+start,'id':str(i)}
            print('positive',int(results[1]))
        else:
            SVM_result[i] = {'positive': 0, 'negative': int(results[0]),'add_negative_nums':add_negtive,
                             'accuracy': float(true_num) / predictions,'ss':i*lenth+start, 'id': str(i)}
    with open(target_path, "w") as file:
        json.dump(SVM_result, file)
with open('dataSets/gene_data_allIndents_ss.pkl', 'rb') as f:
    data_dict = pickle.load(f)
getSS_SMESG(data_dict)
# 不同离差平方和阈值下的SVM预测再生基因
# SVM_20(start=0.01,lenth=0.001,target_path="SVM_models/result0.json")
# SVM_20(start=0.02,lenth=0.0001,target_path="SVM_models/result2.json")
# SVM_20(start=0.02,lenth=0.0001,target_path="SVM_models/result0_1.json")
SVM_20(start=0.02,lenth=0.0001,target_path="SVM_models/result0_0.json")
# # 提取新再生基因
new_grow_gene={}
new_grow_gene_scores=[]
for id in range(21):
    with open('SVM_models/'+str(id)+'_result_gene_classification.pkl', 'rb')as f:
        gene_result=pickle.load(f)
        new_grow_gene[id]=[]
        for gene in gene_result:
            if(gene_result[gene]=='positive'):
                new_grow_gene[id].append(gene)
                new_grow_gene_scores.append(gene)
for id in new_grow_gene:
    print('id:',id,'genes:',new_grow_gene[id])
print(Counter(new_grow_gene_scores))
with open('SVM_models/SVM_gene_classification2.json', 'w')as f:
    json.dump(new_grow_gene,f)
new_grow_gene_score=dict(Counter(new_grow_gene_scores))
print(new_grow_gene_score)


