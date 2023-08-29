## Steps to create computing environment

`SCENIC` is a single-cell specific tool that atempts to use gene expression and known transcription factor netowrks to make inferences about active regulons at the single cell level. The initial implementation was written in R, but is not scalable past ~ 500 cells without using immense amount of computing power. Thus pySCENIC was develop, a faster implmetation.  

The instuctions below are here to provide the spte to generate a functional conda environemnt to complete SCENIC analysis.

```sh
conda create -n scenic_protocol
conda activate scenic_protocol

conda install pandas matplotlib seaborn
conda install -c anaconda cytoolz numpy

conda install scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph louvain
conda install -c conda-forge multicore-tsne
pip install scanpy

pip install pyscenic

pip install --user ipykernel
python -m ipykernel install --user --name=scenic_protocol

#fix a versioning error
pip install numpy==1.21

conda install -c conda-forge gxx_linux-64==3.4.29

```

^^ above ran into conflicts and did not work :(

```sh
conda create -n scPy
conda activate scPy

conda install -c anaconda python=3.6

conda install -c anaconda pandas matplotlib seaborn
conda install -c anaconda cytoolz numpy
conda install -c anaconda scikit-learn statsmodels numba pytables

```
