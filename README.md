# RDPFM-main
![Loading DLow Overview](overview.png "DLow Overview") # 图画好后要补上
! 网址 #网站建好后要把网址放在这里
---
This repo contains the official implementation of our paper:
RDPFM: Regeneration Discovery via Planarian Foundation Model
# Installation 
### Datasets

##### * You can download raw data from NCBI (or preprocessed data and models from [BaiduYun](链接：https://pan.baidu.com/s/1jEreRoycgqEKj70bvbvIcQ?pwd=h5hd)）
### Environment
* **Tested OS:** Linux
* **Main Packages:**
    * Python >= 3.8
    * geneformer == 0.0.1
* **Note**: You can quickly start on RDPFM_min.py or gradually run the codes in planaria_main to learn more details.

### Pretrained Models
* Download our pretrained models at [BaiduYun](https://pan.baidu.com/s/1Ye6bHXcX6lNVMLaXJyzyWg)( password: y9ph) ,then copy the decompressed file to this path and overwrite the original file
* Download raw data from GEO or at [BaiduYun](https://pan.baidu.com/s/1zAonrS7bGnn22pAU3Ggb-g?pwd=i8qs)(password ：i8qs),then copy the decompressed file.
# Train
### unsupervised learning
```
python RDPFM_main.py --act unsupervised_learning
```

### supervised learning
```
python RDPFM_main.py --act supervise_finetune
```

# classify generative states
```
python RDPFM_main.py --act getStates
```
# get gene ebedding

```
python RDPFM_main.py --act getEmbedding
```
# get gene_classification

```
python RDPFM_main.py --act gene_classification
```




