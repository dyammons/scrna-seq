## kallisto pipeline installation and useage on Alpine for scRNA-seq

<br>

### Softare installation
Rather than install kallisto manually, as done for [bulk analysis]() we will be using a wrapper package called `kb-python` which contains `kallisto` and `bustools`. These are the two tools that are needed to go from raw data to count matricies. 

\[Recommended\] You _should_ be able to create a new `env` pretty quickly with the following commands. I have not had success building the `env` using this approach and resorted to using `pip`. Let's call the new environment `kb-tools`.
```sh
conda create -n kb-tools
conda activate kb-tools
conda install -c bioconda kb-python
```

\[Less ideal\] You can install using a `pip` command inside a relavent `conda env` (check out the [source](https://pachterlab.github.io/kallistobustools/) if you have problems with install or usage).
```sh
conda activate <xxx>
pip install kb-python
```
<br>

### Index generation
I would encourage you to build the index at `/projects/$USER/references/canine/` (assuming you are woking with the canine genome).

The index for `kallisto` is unique from traditional aligners (STAR, HISAT2, etc.) in that it uses the cDNA fasta file (transcriptome) instead of the whole genome fasta file. Therefore, you will want to find the cooresponding `.cdna.all.fa.gz` file for the reference of interest. The curl command below will retreive the file required to index the 104 release of CanFam3.1 from ensembl.

```sh
curl -O ftp://ftp.ensembl.org/pub/release-104/fasta/canis_lupus_familiaris/cdna/Canis_lupus_familiaris.CanFam3.1.cdna.all.fa.gz
```

Once the file is retrieved simply run the following command to generate an indexed reference.
NOTE: there is no need to submit a job if you are on a compile or compute node. This is a quick and relatively low resource command.

```sh
kb ref -i canFam3.1.kb-python.idx -g transcripts_to_genes.txt -f1 Canis_lupus_familiaris.CanFam3.1.cdna.all.fa Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa.gz CanFam3.1.104.filtered.gtf

#[2023-11-05 20:18:35,549]    INFO [ref] Preparing Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa.gz, CanFam3.1.104.filtered.gtf
#[2023-11-05 20:19:01,204]    INFO [ref] Splitting genome Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa.gz into cDNA at /projects/dyammons@colostate.edu/references/canine/tmp/tmpw2m8o51u
#[2023-11-05 20:20:19,782]    INFO [ref] Concatenating 1 cDNAs to Canis_lupus_familiaris.CanFam3.1.cdna.all.fa
#[2023-11-05 20:20:20,695]    INFO [ref] Creating transcript-to-gene mapping at transcripts_to_genes.txt
#[2023-11-05 20:20:21,457]    INFO [ref] Indexing Canis_lupus_familiaris.CanFam3.1.cdna.all.fa to canFam3.1.kb-python.idx
```

<br>

### Pseudoalign your samples

Looks like the below block of code will work, convert to job script and run with 8 cores
```sh
kb count -i /projects/$USER/references/canine/canFam3.1.kb-python.idx -g /projects/dyammons@colostate.edu/references/canine/transcripts_to_genes.txt -x 10xv3 \
--h5ad -t 1 --tmp /scratch/alpine/dyammons@colostate.edu/tmp/tmp/ \
-o ../03_output/ \
/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/01_input/duod_norm_1/duod_norm_1_CKDL220032323-1A_HMLG2DSX5_S6_L002_R1_001.fastq.gz /scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/01_input/duod_norm_1/duod_norm_1_CKDL220032323-1A_HMLG2DSX5_S6_L002_R2_001.fastq.gz

#[2023-11-06 18:37:20,869]    INFO [count] Using index /projects/dyammons@colostate.edu/references/canine/canFam3.1.kb-python.idx to generate BUS file to ../03_output/2023-11-06_kb-count_output/duod_norm_1/ from
#[2023-11-06 18:37:20,870]    INFO [count]         /scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/01_input/duod_norm_1/duod_norm_1_CKDL220032323-1A_HMLG2DSX5_S6_L002_R1_001.fastq.gz
#[2023-11-06 18:37:20,870]    INFO [count]         /scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/01_input/duod_norm_1/duod_norm_1_CKDL220032323-1A_HMLG2DSX5_S6_L002_R2_001.fastq.gz
#[2023-11-06 18:56:37,785]    INFO [count] Sorting BUS file ../03_output/2023-11-06_kb-count_output/duod_norm_1/output.bus to /scratch/alpine/dyammons@colostate.edu/tmp/tmp/output.s.bus
#[2023-11-06 18:57:28,182]    INFO [count] Whitelist not provided
#[2023-11-06 18:57:28,184]    INFO [count] Copying pre-packaged 10XV3 whitelist to ../03_output/2023-11-06_kb-count_output/duod_norm_1/
#[2023-11-06 18:57:29,587]    INFO [count] Inspecting BUS file /scratch/alpine/dyammons@colostate.edu/tmp/tmp/output.s.bus
#[2023-11-06 18:57:43,716]    INFO [count] Correcting BUS records in /scratch/alpine/dyammons@colostate.edu/tmp/tmp/output.s.bus to /scratch/alpine/dyammons@colostate.edu/tmp/tmp/output.s.c.bus with whitelist #../03_output/2023-11-06_kb-count_output/duod_norm_1/10x_version3_whitelist.txt
#[2023-11-06 18:58:03,262]    INFO [count] Sorting BUS file /scratch/alpine/dyammons@colostate.edu/tmp/tmp/output.s.c.bus to ../03_output/2023-11-06_kb-count_output/duod_norm_1/output.unfiltered.bus
#[2023-11-06 18:58:25,007]    INFO [count] Generating count matrix ../03_output/2023-11-06_kb-count_output/duod_norm_1/counts_unfiltered/cells_x_genes from BUS file ../03_output/2023-11-06_kb-#count_output/duod_norm_1/output.unfiltered.bus
#[2023-11-06 18:58:48,060]    INFO [count] Reading matrix ../03_output/2023-11-06_kb-count_output/duod_norm_1/counts_unfiltered/cells_x_genes.mtx
#[2023-11-06 18:59:02,647] WARNING [count] 4432 gene IDs do not have corresponding gene names. These genes will use their gene IDs instead.
#[2023-11-06 18:59:02,651]    INFO [count] Writing matrix to h5ad ../03_output/2023-11-06_kb-count_output/duod_norm_1/counts_unfiltered/adata.h5ad

```

### Additional notes that ended up failing to get the task done (this section is here so we don't try it again later in life)
You will need additional software to generate the transcript-to-gene file, do this:
```sh
module load cmake/3.27.7
cd /projects/$USER/software/

git clone https://github.com/BUStools/bustools.git

mkdir build
cd build
cmake ..
make
```

```sh
/projects/$USER/software/bustools/build/src/bustools version
#bustools, version 0.43.1
```

While it is nice to have bustools now... lets install that packages we actually wanted....
```r
remotes::install_github("lambdamoses/BUStoolsR")
library(BUSpaRse)
fasta <- system.file("Canis_lupus_familiaris.CanFam3.1.cdna.all.fa", package = "BUSpaRse")

tr2g_fa <- tr2g_fasta(file = fasta, write_tr2g = FALSE, save_filtered = FALSE)
head(tr2g_fa)
```

The above didn't work, so manually extrat data with:
```sh
grep ">" Canis_lupus_familiaris.CanFam3.1.cdna.all.fa | cut -d" " -f1,4,5,6,7 | sed 's/ *:/ /g' | cut -d" " -f1,3,7,9 | cut -c 2- | sed 's/ /\t/g' > tr2g_k9.tmp
echo -e "transcript\tgene\tgene_biotype\tgene_name" | cat - tr2g_k9.tmp > tr2g_k9.tsv
rm *.tmp
```
