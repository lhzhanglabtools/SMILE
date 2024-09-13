# SMILE

#![image-20221221161146872](SMILE_Overview.png)

## Overview

SMILE is designed for alignment and integration of spatially resolved transcriptomics data.

**a**. SMILE first normalizes the expression proﬁles for all spots and constructs a spatial neighbor network using the spatial coordinates. SMILE further employs a graph constrastive auto-encoder neural network to extract spatially aware embedding, and align the embeddings by minimazing the distribution of anchors, which are identified by MNN. 
**b**. SMILE can be applied to integrate SRT datasets and scRNA-seq to achieve alignment and simultaneous deconvolution of spots from different biological samples.



## Installation
The SMILE package is developed based on the Python libraries [Scanpy](https://scanpy.readthedocs.io/en/stable/), [PyTorch](https://pytorch.org/) and [PyG](https://github.com/pyg-team/pytorch_geometric) (*PyTorch Geometric*) framework, and can be run on CPU or GPU.



First clone the repository. 

```
git clone https://github.com/zhanglabtools/SMILE.git
cd SMILE-main
```

It's recommended to create a separate conda environment for running SMILE:

```
#create an environment called env_SMILE
conda create -n env_SMILE python=3.11

#activate your environment
conda activate env_SMILE
```

Install all the required packages.

```
pip install -r requiements.txt
```
The use of the mclust algorithm requires the rpy2 package (Python) and the mclust package (R). See https://pypi.org/project/rpy2/ and https://cran.r-project.org/web/packages/mclust/index.html for detail.

The torch-geometric library is also required, please see the installation steps in https://github.com/pyg-team/pytorch_geometric#installation

Install SMILE.

```
python setup.py build
python setup.py install
```



## Tutorials

Three step-by-step tutorials are included in the `Tutorial` folder and https://smile.readthedocs.io/en/latest/ to show how to use SMILE. 

- Tutorial 1: Integrating simulation data
- Tutorial 2: Integrating DLPFC slices 
- Tutorial 3: Integrating SRT of healthy and psoriatic skin 

## Support

If you have any questions, please feel free to contact us [zhanglh@whu.edu.cn](mailto:zhanglh@whu.edu.cn). 


