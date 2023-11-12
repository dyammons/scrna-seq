## Instructions for new canine reference build

### Get and prepare files

Pull down the files from ensembl (11/09/2023 @ 8:27 AM).
```sh

rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.110.gtf.gz .
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz .

```

Unzip the files
```sh
gunzip Canis_lupus_familiaris.ROS_Cfam_1.0.110.gtf.gz
gunzip Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz
```

Check for viable biotype options
```sh
grep -oP 'gene_biotype \K\S+' *.gtf | cut -d"\"" -f2 | sort -u
```

<details>
  <summary>Show output</summary>

<br> 

```sh
###output
# IG_C_gene
# IG_V_gene
# TR_C_gene
# TR_J_gene
# TR_V_gene
# Y_RNA
# lncRNA
# miRNA
# misc_RNA
# processed_pseudogene
# protein_coding
# pseudogene
# rRNA
# ribozyme
# scaRNA
# snRNA
# snoRNA
# vault_RNA
```
</details>


Check for chromosomes included in the FASTA
```sh
grep '^>' *.fa | head -n50
```
<details>
  <summary>Show output</summary>

<br> 

```sh
# >1 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:1:1:123313939:1 REF
# >2 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:2:1:86187811:1 REF
# >3 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:3:1:92870237:1 REF
# >4 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:4:1:89007665:1 REF
# >5 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:5:1:89573405:1 REF
# >6 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:6:1:78268176:1 REF
# >7 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:7:1:81039452:1 REF
# >8 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:8:1:75260524:1 REF
# >9 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:9:1:62002293:1 REF
# >10 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:10:1:70361000:1 REF
# >11 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:11:1:75541347:1 REF
# >12 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:12:1:73497294:1 REF
# >13 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:13:1:64037277:1 REF
# >14 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:14:1:61043064:1 REF
# >15 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:15:1:65200600:1 REF
# >16 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:16:1:62021213:1 REF
# >17 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:17:1:65471548:1 REF
# >18 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:18:1:56883407:1 REF
# >19 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:19:1:55265241:1 REF
# >20 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:20:1:58896461:1 REF
# >21 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:21:1:52140716:1 REF
# >22 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:22:1:62106979:1 REF
# >23 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:23:1:53282923:1 REF
# >24 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:24:1:48838997:1 REF
# >25 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:25:1:51941001:1 REF
# >26 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:26:1:40674351:1 REF
# >27 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:27:1:46248802:1 REF
# >28 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:28:1:41862212:1 REF
# >29 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:29:1:42049852:1 REF
# >30 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:30:1:40414903:1 REF
# >31 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:31:1:39518933:1 REF
# >32 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:32:1:39023732:1 REF
# >33 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:33:1:31649084:1 REF
# >34 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:34:1:42263871:1 REF
# >35 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:35:1:26942268:1 REF
# >36 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:36:1:31065185:1 REF
# >37 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:37:1:30932408:1 REF
# >38 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:38:1:24102048:1 REF
# >X dna:primary_assembly primary_assembly:ROS_Cfam_1.0:X:1:127069619:1 REF
# >Y dna:primary_assembly primary_assembly:ROS_Cfam_1.0:Y:1:3937623:1 REF
# >JAAUVH010000248.1 dna:primary_assembly primary_assembly:ROS_Cfam_1.0:JAAUVH010000248.1:1:3974031:1 REF
# ...
```

</details>

From that we can see that the MT and Y chomosomes are not present. This can be confirmed with:
```sh
grep '^>MT' *.fa
# no output

grep '^>Y' *.fa
# >Y dna:primary_assembly primary_assembly:ROS_Cfam_1.0:Y:1:3937623:1 REF

grep '^>X' *.fa 
# >X dna:primary_assembly primary_assembly:UU_Cfam_GSD_1.0:X:1:124992030:1 REF
```

So, to correct for this we will bring over the .gtf and .fa portions from the canFam3.1 MT chromosome data.
(For this the required file are `../canFam31/`)
```sh
grep "^MT" ../canFam31/Canis_lupus_familiaris.CanFam3.1.104.gtf > canFam_mt_append.gtf

cat Canis_lupus_familiaris.ROS_Cfam_1.0.110.gtf canFam_mt_append.gtf > canROS.gtf
```

Pull down the MT data fresh (11.12.2023 @ 8:50 AM)
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna.chromosome.MT.fa.gz .

gunzip Canis_lupus_familiaris.CanFam3.1.dna.chromosome.MT.fa.gz

cat Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa Canis_lupus_familiaris.CanFam3.1.dna.chromosome.MT.fa > canROS_toplevel_mt.fa
```


### mkgtf script


```sh

module load cellranger
cellranger --version
#cellranger cellranger-7.1.0

cellranger mkgtf canROS.gtf canROS_FILTERED.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:TR_C_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_V_gene \
--attribute=gene_biotype:lncRNA \
--attribute=gene_biotype:protein_coding

#Writing new genes GTF file (may take 10 minutes for a 1GB input GTF file)...
#...done
```

### mkref script

```sh
#add a shebang here to make executable

#SBATCH --job-name=cr_mkref

#SBATCH --nodes=1         # this script is designed to run on one node
#SBATCH --ntasks=8
#SBATCH --time=01:00:00   # set time; default = 4 hours; 0 = no restriction

#SBATCH --partition=amilan  # modify this to reflect which queue you want to use. Either 'shas' or 'shas-testing'
#SBATCH --qos=normal      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=dyammons@colostate.edu ### change to your email address ###

#SBATCH --output=cellrngr_mkref_%j.log  # this will capture all output in a logfile with %j as the job #

######### Instructions ###########

### Load cellranger
module purge
module load cellranger
cellranger --version

### Run the command

### This script will generate a reference to be used for cellranger counts ###

### NOTE: ###
# --genome is output directory (so change if you want a different name)
# --fasta is the path to whole genome fasta file
# --genes is a link to the .gtf (use primary assembly from ensembel - dont use all assemblys)
# --memgb specifices RAM; optional, but can be useful for optimization of each genome

#ensure all paths/file names are correct!  Canis_lupus_familiaris.ROS_Cfam_1.0.110.gtf
cellranger mkref --genome=canine_ref_genome_cellranger_7_1_0_ROS_Cfam_1_0_110_base \
                 --fasta=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/ros/canROS_toplevel_mt.fa \
                 --genes=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/ros/canROS_FILTERED.gtf

```
