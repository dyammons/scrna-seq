## Monocle3 installation instructions

If your conda config is set to `strict`, it is recomended to change this setting to `flexible` which will increase likelihood of a successful solve.  
Do that by running the following command:
```sh
conda config --set channel_priority flexible
```

```sh
conda create -n r_monocle3
conda install --conda-forge r
conda install --channel=conda-forge r
conda install -c bioconda r-monocle3
```

When done you may want to set your conda config file back to `strict`.  
(Note: `strict` often results in faster solve times, but can be detrimental if installing software from muliple channels; as done above).
```sh
conda config --set channel_priority strict
```
