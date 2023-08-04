Welcome. The instructions provided here are designed to help you create a reference then run raw (fastq) single cell sequencing data through a cellranger pipeline.

# Steps to create reference genome and run cellranger:
0. [Get to a compile node](#navigate-off-of-a-login-node)
0. [Install cellranger](#install-cell-ranger)
0. [Download and prepare a reference genome](#download-and-prepare-a-reference-genome)
0. [Download and prepare the GTF annotation file](#download-and-prepare-the-gtf-files)
0. [Convert the GTF and genome to a cellranger reference](#convert-gtf-file-and-genome-to-cell-ranger-reference-file)
0. [Get raw data in an accessible location](#get-raw-data-in-place)
0. [Run cellranger counts to align data](#run-cell-ranger-counts)


## Navigate off of a login node
Before getting everything set up, you should navigate off a login node. If you spawned a server using JupyterHub, then you can ssh to the node in which your server spawned.  

<br>

If unsure what node you are on, you can check with `squeue -u $USER`.
```sh
#output
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
2588422    amilan sys/dash dyammons  R      10:38      1 c3cpu-c15-u1-1
```

<br>

The nodelist value associated with the "sys/dash" name is the node that your JupyterHub session is running on and you can move there by running:
```sh
#note: change the node to match the output of your squeue -u $USER command
ssh c3cpu-c15-u1-1
```

<br>

If you ssh'd in using a terminal you will be on a login node and can move to a compile node by running:
```sh
ssh acompile
```
These nodes are designed for light computing and optimizing scripts. For any "real" computing we will be submitting our scripts to a job manager (SLURM).  

When you first launch a server you will be placed on a login node. The login node is designed to ba landing place to get you onto the server and you should move to a compile node if you plan to do any work. Thus, almost everytime you launch a server you should run the above command to get to a compile node.

<br>

## Install Cell Ranger
Questions? Check out 10x Genomics cell ranger [installation page](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).  
I recommend downloading cellranger in your projects space on the Alpine server. Navigate to your desired location then install cellranger.  
Something like this path should work well: `/projects/$USER/software/`.

<br>

#### If you need a hint to get started here is some code:
```sh
mkdir -p /projects/$USER/software/
cd /projects/$USER/software/
```

<br>

#### Download and unpack cellranger:
Here is a command to download cellranger. 10x frequently updates the link, so if it fails to run, check out the 10x [genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) as they may have updated the url.
```sh
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1639406708&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Mzk0MDY3MDh9fX1dfQ__&Signature=a2fdE-2x1h3umjXTQjwakASZyEpEgGhypuS2aL~0gXxUSfZQhG96g66p5faY-WPKQLqhY10n6HrcWgrxo~Oi6IAJfmgqvVTO8JyJvFc5A7M3Mn9~zafNk7OuWOX~gjj-Zqf77RYec1KpjpxBVFQATzCLIXMjn~OVb7Hr~Hwih-74JF9C5jteDwsIB0vkBpiOBOWlsHbb02DkTfpDVcT9d5X9cWYg3rkJRCHqifaJdjpb~wTnrVrwC2e0iS0~G4Qp8anTHB4Tc-RdniMPy8VSGdv4shcUdZGkXOgncXvg8ql1qitz-gJhm3bNrd9xZ6pAgmt~M4623JJE73CiKSDmGA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

<br>

Unpack the tar ball after it finishes downloading.
``` sh
tar -xzvf cellranger*tar.gz
```

<br>

#### Test that you have access to cellranger (be sure to check version number matches with version installed):
```sh
export PATH=/projects/$USER/software/cellranger-6.1.2:$PATH
cellranger
```
If you see a help dialog in your terminal then all should be good. Refer to the 10x website if you would like to do further testing.

To ensure you have access to cellranger when computing, there is an "export" command in each bash script. The path should be correct, but double check to make sure the path is correct.

<br>

## Download and prepare a reference genome:
Navigate to your references directory with `cd /projects/$USER/references/`. Then use the command below to pull down the canine genome. If you are interested in a different genome you can pull down any genome using a similar command, you just need to modify the path according to the ensembl ftp webpage.

Note: when navigating the ensembl ftp website the ftp url will likely lack the word “ensembl” – be sure to add it before “pub” (added to cmd below)
#### Get the reference files:
```sh
#don't forget the "." at the end of the command!
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/fasta/canis_lupus_familiaris/dna/*.dna.toplevel*.fa.gz .
```
Questions? Here is a link to the [ensembl ftp help page](http://uswest.ensembl.org/info/data/ftp/rsync.html).

<br>

#### Ensure files came down correctly
Whenever you retrieve data from an outside source it is always a good idea to check that the data was not altered during transfer.

The way to do this is to check the hash, a 128-bit value that is unique to each file. The value on the ensembl ftp site should be in a file called CHECKSUM, so we will retrevie this file then cross reference the hash with the value of the downloaded file. If a file was altered in any way the MD5 hash will change, making it so you can confirm that your files came through uncorrupted. The following code walks you through the process. NOTE: Ensembl uses unix `sum` command, not `md5sum` to calculate the hash, so you have to do the same to verify the file did not get corrupted in transit.
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/bos_taurus/dna/CHECKSUMS .

grep ".dna.toplevel" CHECKSUMS
# output: $ 00439 799942 Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz

sum *.dna.toplevel*.fa.gz
# output: $ 00439 799942
```
Since this is only one file you can visually inspect to make sure the numbers match from both outputs.  
If the do not match you should delete the file you initially pulled down and re-download the file, as something likely went wrong!

<br>

#### Unzip the genome file:
```sh		
gunzip *.dna.toplevel*.fa.gz
```

<br>

## Download and prepare the GTF files:
Explore the [ensembl ftp website](https://uswest.ensembl.org/info/data/ftp/index.html) to find the annotation (GTF) file you need.  

#### Pull the GTF from ensembl:
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/gtf/canis_lupus_familiaris/*  .
```

<br>

#### Check md5sums:
```sh
grep "Bos_taurus.ARS-UCD1.2.110.gtf.gz" CHECKSUMS
# output: $ 34025 14341 Bos_taurus.ARS-UCD1.2.110.gtf.gz

sum Bos_taurus.ARS-UCD1.2.110.gtf.gz
# output: 34025 14341
```

<br>

#### Prepare the GTF file:
```sh	
gunzip Bos_taurus.ARS-UCD1.2.110.gtf.gz
```

<br>

#### Filter the GTF file with cellranger mkgtf:
Create bash script called “mkgtf.sh” in your references directory:
```sh
touch mkgtf.sh
```
Then copy over the contents of the [mkgtf.sh script](./mkgtf.sh).

If using a Jupyterhub portal then you can use the file navigator panel to locate and open the file. Alternatively you can edit the file using the command `nano mkgtf.sh`. Note: the Jupyterhub text editor is more user friendly.

<br>

Once you have the mkgtf.sh bash script run it with the follwoing command:
```sh
bash mkgtf.sh > mkgtf.log 2>&1 &
```
For reference, the “&” on end of the command makes it so the script runs in the background; check progress with cmd: `jobs -l` (that’s a lowercase L)
		
The output will be a filtered gtf file: "*_filtered.gtf". The goal of this step is to remove unwanted annotations to make subsequent steps easier in terms of file size. The script provided will keep all protein coding annoations as well as other important annotations, such as immunoglobulin genes. The mininium recommended filter is to select all protein coding annotations, the inclusion of additional annotations is optional. At this point, if there are any additional annotations that are not included in the annotation file, you can `cat` them to include them in the alignment process.

<br>

## Convert the gtf file and genome to Cell Ranger reference file:

#### Create the bash and sbatch scripts in your references directory:
```sh
touch mkref.sh cute_cellrngr_mkref.sbatch
```
Copy the contents of [mkref.sh](./mkref.sh) and [cute_cellrngr_mkref.sbatch](./cute_cellrngr_mkref.sbatch) to their respective files then customize them to make sure all the paths/settings match the needs of your run.

<br>

#### Once the files are in place submit a SLURM job using the following cmd: 
```sh
sbatch cute_cellrngr_mkref.sbatch
```	
Should be completed in under 1 hour

<br>

#### You can check progress with this cmd: 
```sh
squeue -u $USER
```
A few notes on cellranger mkref: First, here is a link to the [10x mkref documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references) I recommend looking it over to ensure you understand the process. Second, gtf annotation files contain a fair amount of information in them, but the default settings in cellranger will only look for annotations associated with the feature type of "exon" and ignore all others. In the context of the 10x platform and short read sequencing it is imporant to note there is a strong 3' bias in read mapping, so you may find that you want to include reads that map the 3' untranslated regions (utrs). It is possible to modify the gtf file to extend 3' utrs (currently working on this), but another strategy is to use the output files from cellranger to extend annotations and alter the count matrix to include reads that map to 3' uts - that is where the End Sequencing Analysis Toolkit (ESAT) comes in.

To investigate whether or not you should extend 3' annotations, I recommend looking at metagene plots ([code provided](https://github.com/dyammons/K9-PBMC-scRNAseq/blob/main/analysisCode/metaGenePlot.md), but underdevelopment) to determine how many reads are affected by short 3' utr annotations. From there you can decide how to proceed.

#### You should have a reference when the job finishes!

<br>

## 5. Get raw data in place
From here on out you should be working in your scratch space. 
#### If you are not already in your scratch space you can navigate there with:
```sh
cd /scratch/summit/$USER
```
#### Let's make some directories for organization:
```sh
mkdir project_01
cd project_01
mkdir 01_input 02_scripts
cd 01_input
```
The process of getting your raw data onto the server will vary based on where your data is stored. Regardless you will want to put it in your 01_input directory in your scratch space. You can put it in a sub directory or leave it in the main 01_input directory, just note where you put it!

Useful command to move (pull or push) data:
```sh
rsync -avzP -e 'ssh -p 22' <source path> <user name with "\" before the "@">@login.rc.colorado.edu:/scratch/summit/<user name>/project_01/01_input/
```
The above command will send all the files in the directory you are located in on a local terminal to the server, so just navigate to the directory containing your fastq files then run the command. 

The file name(s) should looks something like this: \<sample name\>_S7_L004_R1_001.fastq.gz

You will need to know the sample name to run cellranger counts.

<br>

## Run Cell Ranger counts
Now that you have everything in place running the final step should be a breeze!

Complete the following step in your 02_scripts directory. If you are not already there use this command: `cd /scratch/summit/$USER/project_01/02_scripts/`.

<br>

#### Create the bash and sbatch scripts to run cellranger counts:
```sh
touch cellrngr_cnts.sh cute_cellrngr_cnts.sbatch
```
Copy the contents of [cellrngr_cnts.sh](./cellrngrCnts.sh) and [cute_cellrngr_cnts.sbatch](./cute_cellrngrCnts.sbatch) to their respective files then customize them to make sure all the paths/options match the needs of your run. Once everything looks good you can submit the SLURM job. 
```sh
sbatch cute_cellrngr_cnts.sbatch
```
The job should take 2-24 hours to run and will create several files for downstream use.

If you do not request enough time for the job you can easily resume. All you have to do is go into the cellranger counts output folder (likely named 'run_counts*') and delete the "_lock" file. Once the file is deleted, you can submit the job again using the same sbatch script.

## Congratulations! You have completed the alignment of your data to the genome and can now get a sense of what the data look like.
