# scCube v1.0.0

## Simulating multiple variability in spatially resolved transcriptomics

[![python >=3.8](https://img.shields.io/badge/python-%3E%3D3.8-brightgreen)](https://www.python.org/) 

scCube is a python package coupled with an interactive website for reproducible, paired, and platform-diverse simulation of spatially-resolved transcriptomic data 

![avatar](images/workflow.jpg)

## Requirements and Installation
[![anndata 0.8.0](https://img.shields.io/badge/anndata-0.8.0-success)](https://pypi.org/project/anndata/) [![numpy 1.23.5](https://img.shields.io/badge/numpy-1.23.5-important)](https://pypi.org/project/numpy/) [![pandas 1.5.3](https://img.shields.io/badge/pandas-1.5.3-critical)](https://pypi.org/project/pandas/) [![scanpy 1.9.1](https://img.shields.io/badge/scanpy-1.9.1-informational)](https://github.com/scverse/scanpy) [![pot 0.8.2](https://img.shields.io/badge/pot-0.8.2-blueviolet)](https://pypi.org/project/POT/) [![matplotlib 3.6.3](https://img.shields.io/badge/matplotlib-3.6.3-ff69b4)](https://pypi.org/project/matplotlib/) [![seaborn 0.12.2](https://img.shields.io/badge/seaborn-0.12.2-9cf)](https://pypi.org/project/seaborn/) [![tqdm 4.64.1](https://img.shields.io/badge/tqdm-4.64.1-lightgrey)](https://pypi.org/project/tqdm/)

### Create and activate Python environment
For scCube, the python version need is over 3.8. If you have installed Python3.6 or Python3.7, consider installing Anaconda, and then you can create a new environment.
```
conda create -n sccube python=3.8
conda activate sccube
```
### Install pytorch
The version of pytorch should be suitable to the CUDA version of your machine. You can find the appropriate version on the [PyTorch website](https://pytorch.org/get-started/locally/).
Here is an example with CUDA11.6:
```
pip install torch --extra-index-url https://download.pytorch.org/whl/cu116
```
### Install other requirements
```
cd scCube-release
pip install -r requirements.txt
```
### Install scCube
```
python setup.py build
python setup.py install
```

## Quick Start
To use scCube we require two formatted `.csv` files as input (i.e. read in by pandas). We have included two test datasets 
in the [tutorial/demo_data folder](tutorial/demo_data) of this repository as examples to show how to use scCube. 

If you generate spot-based ST data, please refer to:
* [Demonstration of scCube on generating spot-based ST data](tutorial/demo_spot.ipynb)

If you generate image-based ST data, please refer to:
* [Demonstration of scCube on generating image-based ST data](tutorial/demo_image.ipynb)

For more details about the format of input and the description of parameters, see the [Tutorial Handbook](tutorial/handbook.md).


## Tutorials
TODO:


## About
scCube was developed by Jie Liao and Jingyang Qian. Should you have any questions, please contact Jingyang Qian at qianjingyang@zju.edu.cn.

