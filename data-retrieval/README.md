## Data retrieval

#### 1) Getting data from NCBI GEO

While you can get the information required to download publically available using web-based portals (NCBI SRA, BioProject etc.) there are also command line tools avalible that may make things easier.  

The below code will create a new `conda` environment called `ncbi_datasets` which will have the software required to query and NCBI resources.

<br>

```sh
#create the env
conda create -n ncbi_datasets

#activate and install software
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli
conda install -c bioconda entrez-direct
conda install -c bioconda perl-time-hires
```

<br>

```sh
#query to get metadata for a project
esearch -db sra -query PRJNA936588 | efetch -format runinfo
```

<details>
<summary>Output</Summary>

<br>

```text
Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
SRR23525335,2023-04-25 07:36:29,2023-02-19 21:24:10,381376085,115175577670,381376085,302,44917,,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos4/sra-pub-zq-3/SRR023/23525/SRR23525335/SRR23525335.lite.1,SRX19417772,GSM7050908,RNA-Seq,cDNA,TRANSCRIPTOMIC SINGLE CELL,PAIRED,0,0,ILLUMINA,Illumina NovaSeq 6000,SRP423368,PRJNA936588,3,936588,SRS16811979,SAMN33356351,simple,9615,Canis lupus familiaris,GSM7050908,,,,,,,no,,,,,"MIP, COLORADO STATE UNIVERSITY",SRA1592663,,public,31B9F3E1DF3D92B84FA6215103D88003,E95FCC036A0BD6713449E0916F3BDBAE
```
</details>

<br>

Let's take this a step further and query for the GEO dataset, as this is were the single-cell count matrices we want typically are found.

<br>

```sh
#query to get metadata for a project
esearch -db gds -query GSE225599 | efetch
```

<details>
<summary>Output</Summary>

<br>

```text
1. A single-cell RNA sequencing atlas of circulating leukocytes from healthy and osteosarcoma affected dogs
(Submitter supplied) We used single-cell RNA sequencing to characterize the heterogeneity of circulating leukocytes in dogs, then employed the dataset to investigate how primary osteosarcoma (OS) tumors impacted circulating leukocytes.
Organism:       Canis lupus familiaris
Type:           Expression profiling by high throughput sequencing
Platform: GPL25760 17 Samples
FTP download: GEO (MTX, RDS, TSV) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/
Series          Accession: GSE225599    ID: 200225599
```
</details>

<br>

From this output we can see the `ftp` path we need to retrieve the count matricies (`ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/`). While we could get all the data here, we might not want everything, so get just the count matricies (usually the only `*.tar` files in the `suppl` directory) we can run the following.

<br>

```sh
#get the tar ball -- note this will download in you current directory
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/*.tar
```

#### 2) Getting data from Novogene

[`getNovoData.sh`](./getNovoData.sh) is a `bash` script that provides a paralellized means of downloading data from Novogene servers and ensures data integretiy after transfer.  

The script requires the input of one `.txt` file that contains the `urls` to the files you wish to download.  

See the script documentation for usage/execution.

#### 3) Getting data from a local data source


The preferred method of getting data on the server is using a `rsync` command since this command should work universally.  

<br>

```sh
#general command structure
rsync -avzP -e 'ssh -p 22' <source path> <user name with "\" before the "@">@login.rc.colorado.edu:/scratch/alpine/<user name>/project_scrna_01/01_input/
```

<br>

```sh
#example where I am moving a .sif file onto my scratch space
#rsync -auzP -e 'ssh -p 22' dyammons_r-env_r.sif dyammons\@colostate.edu@login.rc.colorado.edu:/scratch/alpine/dyammons@colostate.edu
```