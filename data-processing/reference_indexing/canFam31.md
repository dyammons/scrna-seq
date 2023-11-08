## Instructions for canFam3.1 canine reference build

Downloaded on 11/04/2023 via
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/fasta/canis_lupus_familiaris/dna/*.dna.toplevel*.fa.gz .
```

Check for viable biotype options
```sh
grep -oP 'gene_biotype \K\S+' *.gtf | cut -d"\"" -f2 | sort -u
```

### mkgtf script used to filter
```sh
cellranger mkgtf Canis_lupus_familiaris.CanFam3.1.104.gtf Canis_lupus_familiaris.CanFam3.1.104_FILTERED.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:TR_C_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_V_gene
```

### mkref script used to index genome
```sh

### This script will generate a reference to be used for cellranger counts ###

### NOTE: ###
# --genome is output directory (so change if you want a different name)
# --fasta is the path to whole genome fasta file
# --genes is a link to the .gtf (use primary assembly from ensembel - dont use all assemblys)
# --memgb specifices RAM; optional, but can be useful for optimization of each genome

#ensure all paths/file names are correct!
cellranger mkref --genome=canine_ref_genome_cellranger-7.1.0_canFam3.1_base \
                 --fasta=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa \
                 --genes=/scratch/alpine/dyammons@colostate.edu/scRNA_references/canine/Canis_lupus_familiaris.CanFam3.1.104_FILTERED.gtf
                 
```
