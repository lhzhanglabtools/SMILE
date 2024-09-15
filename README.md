# SMILE

## Overview

SMILE is designed for multiscale dissection of spatial heterogeneity by integrating multi-slice spatial and single-cell transcriptomics.

## Installation
### 1. Prepare environment
To install SpaceFlow, we recommend using the [Anaconda Python Distribution](https://anaconda.org/) and creating an isolated environment, so that the SpaceFlow and dependencies don't conflict or interfere with other packages or applications. To create the environment, run the following script in command line:

```bash
conda create -n stsmile_env python=3.7
```

After create the environment, you can activate the `stsmile_env` environment by:
```bash
conda activate stsmile_env
```

### 2. Install SMILE
After successfully installed `Pytorch` with the version that `>=1.9.0`, install the SMILE package using `pip` by:
```bash                                          
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps stSMILE
```

## Tutorial

- [Tutorial: Integrating simulation data](https://github.com/lhzhanglabtools/SMILE/blob/main/tutorial/run_SMILE_on_simulation_data.ipynb)

## Support

If you have any questions, please feel free to contact us [zhanglh@whu.edu.cn](mailto:zhanglh@whu.edu.cn). 


