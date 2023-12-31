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

For this env to solve it is best to set you .condarc file is set to a flexible solver. You can do this with:
```sh
conda config --set channel_priority true
```

Now you should be able to install the required software with the following code block.
```sh
conda create -n scPy
conda activate scPy

conda install -c anaconda python=3.6

conda install -c anaconda pandas matplotlib seaborn
conda install -c anaconda cytoolz numpy
conda install -c anaconda scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph louvain
conda install -c conda-forge multicore-tsne
pip install scanpy

pip install pyscenic

pip install --user ipykernel
python -m ipykernel install --user --name=scPy

#test install
pyscenic -h

#it installed and works, but only with the old tf datasets... scPy2 attempts to get the most recent version running

```



still a no go due to old pyarrow, nee to upgrade eveything...
```sh
conda create -n scPy2
conda activate scPy2

conda install -c anaconda python=3.10

conda install -c anaconda pandas matplotlib seaborn
conda install -c anaconda cytoolz numpy
conda install -c anaconda scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph louvain
conda install -c conda-forge multicore-tsne
pip install scanpy

pip install pyscenic

pip install --user ipykernel
python -m ipykernel install --user --name=scPy2

#test install
pyscenic -h

conda list

#running into this error again: ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.29'
```


Software install above is continuing to fail, so try Docker:

```sh
singularity build aertslab-pyscenic-0.12.1.sif docker://aertslab/pyscenic:0.12.1
singularity build aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif docker://aertslab/pyscenic_scanpy:0.12.1_1.9.1
```
