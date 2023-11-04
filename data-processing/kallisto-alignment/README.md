## kallisto pipeline installation and useage on Alpine for scRNA-seq

### Softare installation
Check the [kallisto manual](https://pachterlab.github.io/kallisto/manual) for most recent installation instructions.  
Additionally, I would encourage you to install the software at `/projects/$USER/software/`.

```sh
#load in the Alpine module to be able to compile the sofware
module load cmake

#make/navigate to the dir to install
mkdir -p /projects/$USER/software/
cd /projects/$USER/software/

#download the file structure
git clone https://github.com/pachterlab/kallisto

#build the toolkit
cd kallisto
mkdir build
cd build
cmake ..
make

#this portion will not work on the server because we cannot install software in /usr/local/bin
# make install
```

### Index generation
I would encourage you to build the index at `/projects/$USER/references/canine/`.

The index for kallisto is unique fomr traditional aligners (STAR, HISAT2, etc.) in that it uses the cDNA fasta file (transcriptome) instead of the whole genome fasta file. Therefore, you will want to find the cooresponding `.cdna.all.fa.gz` file for the reference of interest. The curl command below will retreive the file required to index the 104 release of CanFam3.1 from ensembl.

```sh
curl -O ftp://ftp.ensembl.org/pub/release-104/fasta/canis_lupus_familiaris/cdna/Canis_lupus_familiaris.CanFam3.1.cdna.all.fa.gz
```

Once the file is retrieved simple run the following command.  
NOTE: there is no need to submit a job if you are on a compile or compute node. This is a quick and relatively low resource command.

```sh
kallisto index -i canFam3.1.fa.idx Canis_lupus_familiaris.CanFam3.1.cdna.all.fa
```

### Pseudoalign your samples

code written; need to bring over

