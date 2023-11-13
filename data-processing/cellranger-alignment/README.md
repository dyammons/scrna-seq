Welcome! The instructions provided here are designed to help you create a reference then run raw `.fastq` single-cell RNA sequencing data through a Cell Ranger pipeline.

>Updated ~~August 4, 2023~~ ~~October 30, 2023~~ November 10, 2023 (Designed to be run on Alpine and modified for workshop)

# Steps to create an indexed reference and run cellranger
0. [Get a server launched](#get-a-server-launched)
0. [Install Cell Ranger](#install-cell-ranger-now-optional)
0. [Download and prepare a reference genome](#download-and-prepare-a-reference-genome)
0. [Get raw data in an accessible location](#get-raw-data-in-place)
0. [Run Cell Ranger counts to align data](#run-cell-ranger-counts)


## Get a server launched

To launch a server using JupyterLab, visit https://ondemand-rmacc.rc.colorado.edu/ and login using your CSU credentials.
In the OnDemand portal we will launch in interactive Jupyter session to gain access to Alpine (check [this](https://curc.readthedocs.io/en/latest/gateways/OnDemand.html) out for more details).

<br>

<details>
  <summary>Additional information that we will skip over</summary>  

<br>

When you first launch a server you will be placed on a login node. The login node is designed to be landing place to get you onto the server and you should move to a compute/compile node if you plan to do any work. Thus, almost everytime you launch a server you should run the following commands to get off the login node. For today we will not be doing any computational tasks on the node, so we skip over this.

If you spawned a server using JupyterHub, then you can ssh to the node in which your server spawned.  

<br>

If unsure what node you are on, you can check with `squeue -u $USER`.
```sh
#output
#JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
#3672920  acompile sys/dash dyammons  R      11:49      1 c3cpu-a7-u34-4
```

<br>

The nodelist value associated with the "sys/dash" name is the node that your JupyterHub session is running on and you can move there by running:
```sh
#note: change the node to match the output of your squeue -u $USER command
ssh c3cpu-a7-u34-4
```

</details>

<br>

These sessions are designed for light computing, compiling, and optimizing scripts. For any "real" computing we will be submitting our scripts to a job manager (SLURM).  

<br>

<details>
  <summary>Click me for an overview of Alpine file systems.</summary>  

<br>

```sh
#2GB permanent storage where we start out
/home/$USER/

#250GB permanent storage -- where we will put reference
/projects/$USER/

#10TB temporary storage (deleted every 90 days) -- where raw data goes
/scratch/alpine/$USER/
```

</details>

<br>

## Install Cell Ranger (Now optional)

`Cell Ranger` is now installed in a software module on Alpine, so there is no need to install the software manually. If you want to use an older version then you may have to complete a manual install following the instructions below. Otherwise, you can simply load the software by running the following:
```sh
module load cellranger
```

To ensure everything is functioning properly run:
```sh
cellranger --version
```

and/or

```sh
cellranger --help
```

<details>
  <summary>Click for manual install instructions</summary>

<br>
 
Questions? Check out 10x Genomics cell ranger [installation page](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).  
I recommend downloading cellranger in your projects space on the Alpine server. Navigate to your desired location then install cellranger.  
Something like this path should work well: `/projects/$USER/software/`.

If you need a hint to get started here is some code:
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
</details>

<br>

## Download and prepare a reference genome

For this workshop we will be using pre-indexed references, so we will essentially be skipping this step today. This is a key part esspecially when working in non-traditional animal models. The hidden code below walks through the process of how you would go about generating an indexed reference for canFam3.1. The specific steps taken to generate each of the references avalible for use today can be found in [:file\_folder: reference-indexing](/data-processing/reference-indexing).

<br>

To get the references in place we will create a [symbolic link](https://linuxize.com/post/how-to-create-symbolic-links-in-linux-using-the-ln-command/) in a `references` directory located in out `projects` space.

```sh
mkdir -p /projects/$USER/references/canine/
cd /projects/$USER/references/canine/
```

<details>
	<summary>Bovine</summary>

<br>
 
```sh
mkdir -p /projects/$USER/references/bovine/
cd /projects/$USER/references/bovine/
```
</details>

<br>

Three different commands are provided to link up a reference, so use the one for the genome you wish to use for alignment.

<details>
	<summary>canFam3.1</summary>
	
```sh
ln -sf /scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/canFam31/canine_ref_genome_cellranger-7.1.0_canFam3.1_base/ canine_ref_genome_cellranger-7.1.0_canFam3.1_base
reference=/projects/$USER/references/canine/canine_ref_genome_cellranger-7.1.0_canFam3.1_base
```

</details>


<details>
  <summary>canFam4</summary>
	
```sh
ln -sf /scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/gsd/canine_ref_genome_cellranger_7_1_0_gsd_UU_Cfam_GSD_1_0_110_base/ canine_ref_genome_cellranger_7_1_0_gsd_UU_Cfam_GSD_1_0_110_base
reference=/projects/$USER/references/canine/canine_ref_genome_cellranger_7_1_0_gsd_UU_Cfam_GSD_1_0_110_base
```

</details>

<details>
  <summary>ROS</summary>

```sh
ln -sf /scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/ros/canine_ref_genome_cellranger_7_1_0_ROS_Cfam_1_0_110_base/ canine_ref_genome_cellranger_7_1_0_ROS_Cfam_1_0_110_base
reference=/projects/$USER/references/canine/canine_ref_genome_cellranger_7_1_0_ROS_Cfam_1_0_110_base
```

</details>

<br>

<details>
  <summary>Bovine</summary>

```sh
ln -sf /projects/dyammons@colostate.edu/references/bovine/bovine_ref_genome/ bovine_ref_genome
reference=/projects/$USER/references/bovine/bovine_ref_genome
```

</details>

<br>

Once the files are properly linked we should now be able to see it, so let's check.
```sh
ls
#canine_ref_genome_cellranger-7.1.0_canFam3.1_base #should be a teal color
```

<br>

<details>
  <summary>If all else fails we will point directly to my scratch space</summary>

CanFam4
```sh
reference=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/gsd/canine_ref_genome_cellranger_7_1_0_gsd_UU_Cfam_GSD_1_0_110_base/
```

CanFam3.1
```sh
reference=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/canFam31/canine_ref_genome_cellranger-7.1.0_canFam3.1_base/
```

ROS
```sh
reference=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/ros/canine_ref_genome_cellranger_7_1_0_ROS_Cfam_1_0_110_base/
```

Bovine
```sh
reference=/projects/dyammons@colostate.edu/references/bovine/bovine_ref_genome/
```

</details>


<details>
  <summary>Bonus: Click for instructions to pull down required files and generate a genomic index</summary>

<br>

Navigate to your references directory with `cd /projects/$USER/references/canine/`. Then use the command below to pull down the canine genome. If you are interested in a different genome you can pull down any genome using a similar command, you just need to modify the path according to the ensembl ftp webpage.

Note: when navigating the ensembl ftp website the ftp url will likely lack the word “ensembl” – be sure to add it before “pub” (added to cmd below)  

<br>

Download the reference files from ensembl:
```sh
#don't forget the "." at the end of the command!
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/fasta/canis_lupus_familiaris/dna/*.dna.toplevel*.fa.gz .
```
Questions? Here is a link to the [ensembl ftp help page](http://ensembl.org/info/data/ftp/rsync.html).

<br>

Ensure files came down correctly:
Whenever you retrieve data from an outside source it is always a good idea to check that the data was not altered during transfer.

The way to do this is to check the hash, a 128-bit value that is unique to each file. The value on the ensembl ftp site should be in a file called CHECKSUM, so we will retrevie this file then cross reference the hash with the value of the downloaded file. If a file was altered in any way the hash will change, making it so you can confirm that your files came through uncorrupted. The following code walks you through the process. NOTE: Ensembl uses unix `sum` command, not `md5sum` to calculate the hash, so you have to do the same to verify the file did not get corrupted in transit.
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/fasta/canis_lupus_familiaris/dna/CHECKSUMS .

grep ".dna.toplevel" CHECKSUMS
# output: $ 32065 708687 Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa.gz

sum *.dna.toplevel*.fa.gz
# output: $ 32065 708687
```
Since this is only one file you can visually inspect to make sure the numbers match from both outputs.  
If they do not match you should delete the file you initially pulled down and re-download the file, as something likely went wrong!

<br>

Unzip the genome file:
```sh		
gunzip *.dna.toplevel*.fa.gz
```

<br>

#### Download and prepare the GTF files:
Explore the [ensembl ftp website](https://uswest.ensembl.org/info/data/ftp/index.html) to find the annotation (GTF) file you need.  

Download the GTF from ensembl:
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/gtf/canis_lupus_familiaris/*  .
```

<br>

Check sums to ensure they were not corrupted in transit:
```sh
grep "Canis_lupus_familiaris.CanFam3.1.104.gtf.gz" CHECKSUMS
# output: $ 61947 17598 Canis_lupus_familiaris.CanFam3.1.104.gtf.gz

sum Canis_lupus_familiaris.CanFam3.1.104.gtf.gz
# output: 61947 17598
```

<br>

Unzip the GTF file:
```sh
gunzip Canis_lupus_familiaris.CanFam3.1.104.gtf.gz
# rm *.gtf.gz #uncomment and run if you want to remove unnecessary files
```

<br>

#### Filter the GTF file with cellranger mkgtf:
Create a bash script called “mkgtf.sh” in your `/projects/references/canine/` directory:
```sh
touch mkgtf.sh
```
Then copy over the contents of the [mkgtf.sh script](./mkgtf.sh).

If using a Jupyterhub portal then you can use the file navigator panel to locate and open the file. Alternatively you can edit the file using the command `nano mkgtf.sh`. Note: the Jupyterhub text editor is more user friendly.

<br>

The goal of this step is to remove unwanted annotations to make subsequent steps easier in terms of file size. The script provided will keep all `protein coding` annoations as well as a few other important annotations, such as `immunoglobulin genes`. The mininium recommended filter is to select all `protein coding` annotations, the inclusion of additional annotations is optional. At this point, if there are any additional annotations that are not included in the annotation file, you can `cat` them to include them in the alignment process.  

To check what biotypes are present in the gtf file you can run:
```sh
grep -oP 'gene_biotype \K\S+' *.gtf | cut -d"\"" -f2 | sort -u

###output:
#IG_C_gene
#IG_V_gene
#TR_C_gene
#TR_J_gene
#TR_V_gene
#protein_coding
```

If it turns out all the biotypes are ones that you want included (as is the case above) then this step really isn't necessary, but no harm in running it either.

<br>

Once you have the mkgtf.sh bash script run it with the following command:
```sh
bash mkgtf.sh > mkgtf.log 2>&1 &
```
For reference, the `&` on end of the command makes it so the script runs in the background; check progress with cmd: `jobs -l` (that’s a lowercase L)
		
The output should be a filtered gtf file: `*_FILTERED.gtf`. 

<br>

### Convert the gtf file and genome to Cell Ranger reference file

Create the sbatch script in your /projects/references/canine/ directory:
```sh
touch cute_cellrngr_mkref.sbatch
```
Copy the contents of [cute_cellrngr_mkref.sbatch](./cute_cellrngr_mkref.sbatch) to the newly created file then customize it to ensure all the paths/options match the needs of your run.

<br>

Once the files are in place submit a SLURM job using the following command: 
```sh
sbatch cute_cellrngr_mkref.sbatch
```	
Should be completed in under 1 hour.

<br>

You can check progress with this cmd: 
```sh
squeue -u $USER
```
A few notes on cellranger mkref:  
First, here is a link to the [10x mkref documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references) I recommend looking it over to ensure you understand the process.  

Second, gtf annotation files contain a fair amount of information in them, but the default settings in cellranger will only look for annotations associated with the feature type of `exon` and ignore all others. In the context of the 10x platform and short read sequencing it is imporant to note there is a strong 3' bias in read mapping, so you may find that you want to include reads that map to `three_prime_utr` (3' untranslated regions). It is possible to modify the gtf file to convert all `three_prime_utr` data points to `exon`, however my brief testing has revealed this has a profoundly negative impact on the alignment quality.

Third, there are a few tool kits that will extend annoations in the 3' direction to increase alignment. The tool I have used is [End Sequencing Analysis Toolkit](https://genome.cshlp.org/content/26/10/1397.full.html) (ESAT), but I am not a huge fan of this tool and I do not think this step is necessary.  

If you're curious about how strong the 3' bias is, I recommend looking at metagene plots ([code provided](https://github.com/dyammons/K9-PBMC-scRNAseq/blob/main/analysisCode/metaGenePlot.md), but the code is underdeveloped/abandoned) to determine how many reads are affected by short 3' utr annotations. From there you can decide if you want to take additional steps to handle the 3' bias.

#### You should have a reference when the job finishes!

</details>

<br>

## Get raw data in place
From here on out you should be working in your `scratch` space.  
If you are not already in your scratch space you can navigate there with:
```sh
cd /scratch/alpine/$USER
```

<br>

We will be using a pretty standard directory structure which uses the main directories `project_scrna_01`, `01_input`, and `02_scripts`.

<details>
<summary>Directory tree</summary>

<br>

Personal structure
```sh
project_scrna_01/
├── 01_input -> /scratch/alpine/dyammons@colostate.edu/proj02_k9_pbmc/01_input/
└── 02_scripts
    └── cute_cellrngr_cnts.sbatch
```

Structure of symbolic link input directory
```sh
/scratch/alpine/dyammons@colostate.edu/proj02_k9_pbmc/01_input/
├── Healthy_PBMC_2
│   ├── Healthy_PBMC_2_S7_L004_R1_001.fastq.gz
│   └── Healthy_PBMC_2_S7_L004_R2_001.fastq.gz
├── ...
└── k9_PBMC_OSA_3
    ├── MD5.txt
    ├── k9_PBMC_OSA_3_CKDL210024467-1a-SI_TT_A2_HVWVGDSX2_S2_L002_I1_001.fastq.gz
    ├── k9_PBMC_OSA_3_CKDL210024467-1a-SI_TT_A2_HVWVGDSX2_S2_L002_I2_001.fastq.gz
    ├── k9_PBMC_OSA_3_CKDL210024467-1a-SI_TT_A2_HVWVGDSX2_S2_L002_R1_001.fastq.gz
    └── k9_PBMC_OSA_3_CKDL210024467-1a-SI_TT_A2_HVWVGDSX2_S2_L002_R2_001.fastq.gz

6 directories, 23 files
```
 
</details>

<br>

Next, let's make the directories:
```sh
mkdir -p ./project_scrna_01/02_scripts/
cd project_scrna_01
```

<br>

Now we are ready to get raw data into place. For today, we will be creating a symbolic link to the `.fastq` files in my `scratch` space.

It is worth noting, that in real life application the process of getting your raw data onto the server will vary based on where your data is stored. There are several mechanisms to get your data onto the server, but I tend to use an `rsync` command or [Globus](https://app.globus.org/file-manager). Regardless of transfer method you will want to put it in your 01_input directory in your scratch space with each sample in its own subdirectory within `01_input`.

<br>

<details open>
  <summary>Symbolic link to canine PBMC data</summary>

<br>

```sh
ln -sf /scratch/alpine/dyammons@colostate.edu/proj02_k9_pbmc/01_input/ 01_input
```

</details>


<details>
  <summary>Symbolic link to canine nasal lavage data (LC and TL only)</summary>

<br>

```sh
ln -sf /scratch/alpine/dyammons@colostate.edu/proj04_k9_nasal/01_input/ 01_input
```

</details>

<br>

<details>
  <summary>Symbolic link to bovine data (SR and GJ only) </summary>

<br>

```sh
ln -sf /scratch/alpine/dyammons@colostate.edu/proj08_bov_scrna/01_input/ 01_input
```

</details>

<br>

Once the files are properly linked we should now be able to see them, so let's check.
```sh
ls ./01_input/
#Healthy_PBMC_2  OSA  OSA_PBMC_2  healthy_PBMC_4  healthy_pbmc_6  k9_PBMC_Healthy_3  k9_PBMC_OSA_3

ls ./01_input/Healthy_PBMC_2/
#Healthy_PBMC_2_S7_L004_R1_001.fastq.gz  Healthy_PBMC_2_S7_L004_R2_001.fastq.gz

ls ./01_input/healthy_pbmc_6/
#MD5.txt
#healthy_pbmc_6_CKDL220005143-1a-SI_TT_H1_HGFT3DSX3_S5_L002_I1_001.fastq.gz
#healthy_pbmc_6_CKDL220005143-1a-SI_TT_H1_HGFT3DSX3_S5_L002_I2_001.fastq.gz
#healthy_pbmc_6_CKDL220005143-1a-SI_TT_H1_HGFT3DSX3_S5_L002_R1_001.fastq.gz
#healthy_pbmc_6_CKDL220005143-1a-SI_TT_H1_HGFT3DSX3_S5_L002_R2_001.fastq.gz
```

<details>
  <summary>If all else fails we will point directly to my scratch space</summary>

Canine PBMC
```sh
input_dir=/scratch/alpine/dyammons@colostate.edu/proj02_k9_pbmc/01_input/
```

Canine nasal lavage
```sh
input_dir=/scratch/alpine/dyammons@colostate.edu/proj04_k9_nasal/01_input/
```
</details>

<details>
  <summary>If you were transfering your own data to the server I would use the following rsync command</summary>

 <br>
 
Useful command to move (pull or push) data:
```sh
rsync -avzP -e 'ssh -p 22' <source path> <user name with "\" before the "@">@login.rc.colorado.edu:/scratch/alpine/<user name>/project_scrna_01/01_input/

#example where I am moving a .sif file onto my scratch space
#rsync -auzP -e 'ssh -p 22' dyammons_r-env_r.sif dyammons\@colostate.edu@login.rc.colorado.edu:/scratch/alpine/dyammons@colostate.edu
```

The above command will send all the files in the directory you are located in on a local terminal to the server, so just navigate to the directory containing your `.fastq` files then run the command. 

Again, if you do not want to use an `rsync` command I highly recommend using [Globus](https://app.globus.org/file-manager). Globus is a much better option than FileZilla as Globus will check file hashes and repeatedly try if transfer initially fails - this is not the case with FileZilla.

</details>

<br>

## Run Cell Ranger counts
Now that you have everything in place running the final step _should_ be a breeze!

Navigate to the `02_scripts` directory:
```sh
cd /scratch/alpine/$USER/project_scrna_01/02_scripts/
```

<br>

Create the sbatch script to run `cellranger count`:
```sh
touch cute_cellrngr_cnts.sbatch
```

<br> 

Print out the path to the reference index that we stashed ealier so we can ensure that we set the correct path in the `SBATCH` file:
```sh
#path to reference directory
echo $reference
```

Copy the contents of [cute_cellrngr_cnts.sbatch](./cute_cellrngr_cnts.sbatch) to the newly created file then customize the user preferences section to make sure all the paths/options match the needs of your run. Once everything looks good we will submit a test job to ensure we have a fair chance of getting the job to run first try.

Submit the job with this command...
```sh
sbatch cute_cellrngr_cnts.sbatch
```

and check the progress of the job by running:
```sh
squeue -u $USER
```

<br> 

Now that we know the script at least has the correct paths, let's update the file for a real run.

<details>
  <summary>Show real job script</summary>

Here is the script we will be using for the run. Carefully copy over the contents that are missing from the test `.sbatch` file.
```sh
#!/usr/bin/env bash

#SBATCH --job-name=cellrngr_cnt
#SBATCH --ntasks=24       # modify this number to reflect how many cores you want to use (up to 64)
#SBATCH --nodes=1         # this script is designed to run on one node
#SBATCH --time=06:00:00   # set time; default = 4 hours

#SBATCH --partition=amilan  # modify this to reflect which queue you want to use. Either 'amilan' or 'atesting'
#SBATCH --qos=normal      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=dyammons@colostate.edu ### change to your email ###

#SBATCH --output=cellrngr_cnt_%A_%a.log  #modify as desired - will output a log file where the "%A" inserts the job ID number and the %a

#SBATCH --array=0-7 #set this to 0-(# of samples - 1), so the example is for 8 samples -- if you are only running 1 sample, then you can set it to 0-0



##### USER SETTINGS #####

#specify path to reference
ref_pwd=/projects/dyammons@colostate.edu/references/canine/canine_ref_genome_cellranger_7_1_0_ROS_Cfam_1_0_110_base


#specify sample names to run
samples=$(ls -lh ../01_input/ | grep "^d" | awk '{print $9}')
declare -a StringArray=($samples)

#if you have extra dirs or only want to run select samples, then store the sample names in the StringArray variable
#declare -a StringArray=("sample1" "sample2" "sample3")

##### END SETTINGS #####


##### SET ENV #####

### Load cellranger
module purge
module load cellranger
cellranger --version


##### BEGIN CODE #####

### Excute cellranger count
sampleName=$(ls ../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ | grep "fastq.gz" | head -n1 | awk -F "_S" '{print $1}')

cmd1="cellranger count --id=${StringArray[${SLURM_ARRAY_TASK_ID}]} \
                       --fastqs=../01_input/${StringArray[${SLURM_ARRAY_TASK_ID}]}/ \
                       --sample=${sampleName} \
                       --transcriptome=${ref_pwd} \
                       --expect-cells=5000" ### MODIFY as needed
echo $cmd1
echo -e "\t$ ${cmd1}"
time eval $cmd1

##### END CODE #####
```

</details>

<br> 

Each job should take 2-24 hours to run and will create several files for downstream use.

<br>

<details>
  <summary>Click me if your job times out</summary>

If you do not request enough time for the job (default for the code we are using is 6 hours) you can fairly easily resume the run. All you have to do is go into the `cellranger count` output folder for the sample that did not finish and delete the `_lock` file.

Once that file is deleted you will be able to resume the run. Before submitting the job again you will want to change the `StringArray` variable to store the sample(s) that need to be resumed.  
(NOTE the changes in `#` useage in the code block below)
```sh
##### Load in sample names to run #####
#samples=$(ls -lh ../01_input/ | grep "^d" | awk '{print $9}')
#declare -a StringArray=($samples)

#if you have extra dirs or only want to run select samples, then store the sample names in the StringArray variable
declare -a StringArray=("sample1" "sample2" "sample3") #replace with sample names
```

You also will have to change the `SBATCH` directive `array` to match the number of sample you are running.
```sh
#if three need to be be resumed you would change as follows
#SBATCH --array=0-2
```

Make the above changes to the `.sbatch` file, then resume the run with `sbatch cute_cellrngr_cnts.sbatch`.

</details>

## Congratulations! You have completed the alignment of your data to the genome and can now get a sense of what the data look like.
