### These instructions are designed to help you create a conda envrionemnt to run the R based analysis code on the SUMMIT server.
#### First, create a .condarc file if there is not already one in your home directory
```sh
nano .condarc
```
#### If not already in the .condarc file, add the following contents:
```sh
pkgs_dirs:
  - /projects/$USER/.conda_pkgs
envs_dirs:
  - /projects/$USER/software/anaconda/envs
channels:
  - conda-forge
  - bioconda
  - r
  - defaults
channel_priority: strict
```
The .condarc code will specify where to save you Conda environment contents (in your projects space) and will set the priority order for Conda. The strict order is necessary to get the most recent version of R through conda-forge.

#### Next, make sure the following is in your ".bashrc" file:
```sh
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# User specific aliases and functions

export PS1='[\h \w]\$ '

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/curc/sw/anaconda3/2019.03/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/curc/sw/anaconda3/2019.03/etc/profile.d/conda.sh" ]; then
        . "/curc/sw/anaconda3/2019.03/etc/profile.d/conda.sh"
    else
        export PATH="/curc/sw/anaconda3/2019.03/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```
This will initiate your "base" Conda environment when launching a server.


While you're here you can also add some alias' to make navigating the server easier:
```sh
alias que='squeue -u $USER'
alias scratch='cd /scratch/summit/$USER'
alias project='cd /projects/$USER/'
```

Side note... if you ever want to activate a conda environment through a SLURM job you can activate conda with the following code:
```sh
source /curc/sw/anaconda3/latest
```

#### Now that Conda is installed, let's create a custom environment.
```sh
conda create -n r_env r-base r-essentials
```
This will create an R conda environemnt called "r_env" with some basic packages. It might take a little while to solve the environment, but once it is solved, accept the package plan and all the packages will be installed.

#### Activate the newly created environemnt with:
```sh
conda activate r_env
```

#### Whenever in this Conda environemnt you can launch an R session by simply typing "R" and you will be in an interactive R session.
```sh
R
```
Now you are in R and can begin installing packages as you would in an Rgui or Rstudio. One you will definitely need is Seurat.

#### Use this command to install Seurat:
```r
install.packges("Seurat")
```
Occasionally, R will have trouble installing dependencies, so if you run into trouble look at the error message and try to install any failed dependencies on their own. Then retry to main package you wish to install (in this case Seurat).


#### To run the first step of the Seurat analysis pipeline you will also need DoubletFinder:
```r
install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
```


Note: Sometimes packages will be unable to be installed through R, so exit R, and try to do it through conda install. An example (no need to run this now!):
```sh
conda install -c conda-forge r-hdf5r r-httpuv
```

### You should now have a Conda environment that can run most of the analysis code. Enjoy!
