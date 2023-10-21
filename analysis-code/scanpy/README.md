## Create a Conda env then convert to a kernel

Create the env
```sh
conda create -n scanpy
conda activate scanpy
```

Install analysis software
```sh
conda install -c conda-forge scanpy
conda install -c conda-forge python-igraph leidenalg
conda install -c conda-forge python-annoy
```

Install ipykernel to get conda env into a kernel
```sh
conda install -c anaconda ipykernel
```

Export the env as a kernel and enjoy
```sh
python -m ipykernel install --user --name=scanpy
```
