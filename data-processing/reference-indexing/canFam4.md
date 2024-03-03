## Instructions for new canine reference build


Pull down the files from ensembl (11/08/2023 @ 5:36 PM).
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/canis_lupus_familiarisgsd/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.110.gtf.gz .
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/canis_lupus_familiarisgsd/dna/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa.gz .
```

Check for viable biotype options
```sh
grep -oP 'gene_biotype \K\S+' *.gtf | cut -d"\"" -f2 | sort -u
```

Check for chromosomes included in the FASTA
```sh
grep '^>' *.fa | head -n50
```

From that we can see that the MT and Y chomosomes are not present. This can be confirmed with:
```sh
grep '^>MT' *.fa
#no output
grep '^>Y' *.fa
#no output

grep '^>X' *.fa 
#>X dna:primary_assembly primary_assembly:UU_Cfam_GSD_1.0:X:1:124992030:1 REF
```

So, to correct for this we will bring over the .gtf and .fa portions from canFam3.1. Since Y isn't in CanFam3.1 either we will just bring over the MT chromosome data.
(For this the required file are `../canFam31/`)
```sh
grep "^MT" ../canFam31/Canis_lupus_familiaris.CanFam3.1.104.gtf > canFam_mt_append.gtf

cat Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.110.gtf canFam_mt_append.gtf > canGSD.gtf
```

```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna.chromosome.MT.fa.gz .

gunzip Canis_lupus_familiaris.CanFam3.1.dna.chromosome.MT.fa.gz

cat Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa Canis_lupus_familiaris.CanFam3.1.dna.chromosome.MT.fa > canGSD_toplevel_mt.fa
```


### mkgtf script


```sh

module load cellranger
cellranger --version
#cellranger cellranger-7.1.0

cellranger mkgtf canGSD.gtf canGSD_FILTERED.gtf \
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

#!/usr/bin/env bash

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

#ensure all paths/file names are correct!
cellranger mkref --genome=canine_ref_genome_cellranger_7_1_0_gsd_UU_Cfam_GSD_1_0_110_base \
                 --fasta=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/gsd/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa \
                 --genes=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/gsd/canGSD_FILTERED.gtf
                
```
