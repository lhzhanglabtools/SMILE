# SMILE

![SMILE_Overview](https://github.com/lhzhanglabtools/SMILE/blob/main/SMILE_overview.png)

## Overview

SMILE is designed for alignment and integration of spatially resolved transcriptomics data.


## Installation
The SMILE package is developed based on the Python libraries [Scanpy](https://scanpy.readthedocs.io/en/stable/), [PyTorch](https://pytorch.org/) and [PyG](https://github.com/pyg-team/pytorch_geometric) (*PyTorch Geometric*) framework.

First clone the repository. 

```
git clone https://github.com/lhzhanglabtools/SMILE.git
cd SMILE-main
```

It's recommended to create a separate conda environment for running SMILE:

```
#create an environment called env_SMILE
conda create -n env_SMILE python=3.11

#activate your environment
conda activate env_SMILE
```

Install all the required packages. The torch-geometric library is required, please see the installation steps in https://github.com/pyg-team/pytorch_geometric#installation
```
conda install pyg
conda install conda-forge::pytorch_scatter
conda install conda-forge::pytorch_cluster
conda install conda-forge::pytorch_sparse
```

The use of the mclust algorithm requires the rpy2 package (Python) and the mclust package (R). See https://pypi.org/project/rpy2/ and https://cran.r-project.org/web/packages/mclust/index.html for detail.

```
pip install -r requirements.txt
```

Install SMILE. We provide two optional strategies to install SMILE.
```
pip install stSMILE
```
Or

```
python setup.py build
python setup.py install
```



## Tutorials

Three step-by-step tutorials are included in the `Tutorial` folder to show how to use  SMILE. 

- [Tutorial 1: Running SIMLE on simulation data]()
- [Tutorial 2: Running SMILE on DLPFC slices]()
- [Tutorial 2 plus: Running semi-SMILE on DLPFC slices]()
- [Tutorial 3: Running SMILE on anterior and posterior sections of mouse brain]()
- [Tutorial 4: Running SMILE on normal skin and psoriasis diseased skin]()
- [Tutorial 5: Running SMILE on SRT data from Stereo-seq and Slide-seqV2 platforms]()

## Support

If you have any questions, please feel free to contact us [zhanglh@whu.edu.cn](mailto:zhanglh@whu.edu.cn). 


