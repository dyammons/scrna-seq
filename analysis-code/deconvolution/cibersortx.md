# Deconvolution using CIBERSORTx

The goal of this tutorial is to go over how to extract the data required to generate a deconvolution reference from a single-cell RNA seq Seurat object. Then use the CIBERSORTx docker container to evaluate the performance using single cell data as the ground truth.

### Create a CIBERSORTx account
Visit the [CIBERSORTx](https://cibersortx.stanford.edu/) website and create an account. You will have to use a `.edu` email account.  
Once created, request access to the docker container by visiting the MENU>Download page [link](https://cibersortx.stanford.edu/download.php).  
After a few day you should have permissions to use the container and installation instructions should be available.  
Below are quick instructions to retrieve and use the docker container, but refer to their instructions if ever in question.

### Set up the container
The container was made using docker, but on Alpine we only have access to a closely related contaerization tool called Singularity.  
The default instructions provided by the CIBERSORT team will not directly work with Singularity, so slight modifications are neccsiary.

First we will use `singularity pull` instead of `docker pull` to retieive the container.
```sh
# In the desired location run:
singularity pull docker://cibersortx/fractions
```

This command should pull down the container and create a file called `fractions_latest.sif`.  
If errors arise with the command, you may have to configure your Singularity cache path.

Once the container is availible, you now have all the software you need to run CIBERSORTx, so let's get some data to run!

### Prepare a single cell reference
The single cell reference used to generate the deconolution matrix can be derrived from any annotated single cell dataset.  
All we need is the count matrix in a tab-deliminated `.txt` file that has genes as rows and columns as cell type annotation for each cell (a gene X cell matrix).
To generate this we will use a publically available dataset, but the approach should work for any dataset processed using `Seurat`.
```r

library(Seurat)
seu.obj <- readRDS("../../bov_lav_scRNA/output/s3/allCells_clean_viral_S3.rds")

#Select cell type annotation level to use for deconvolution
ct_level <- "celltype.l2"
#Down sample to make data easier to work with (keep equal number of cells)
Idents(seu.obj) <- ct_level
seu.obj <- subset(
    x = seu.obj, downsample = min(table(seu.obj@meta.data[[ct_level]]))
)
#Get the data - this can use a lot of RAM if there are still a lot of cells
#Here we get log normalized counts
cnt_mat <- FetchData(
    seu.obj, vars = rownames(seu.obj), assay = "RNA", layer = "data"
)
#Transpose the data to get it in genes (rows) X cells (columns)
cnt_mat <- t(cnt_mat)
#Replace the cell barcodes with the cell type name
colnames(cnt_mat) <- as.character(seu.obj@meta.data[[ct_level]])
#Save the matrix
write.table(
    cnt_mat, quote = FALSE, file = paste0("../output/scrna_", 
                                          gsub("\\.", "_", ct_level), ".txt")
)
```

### Run the deconvolution algo

```sh
singularity exec -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_data:/src/data -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_output:/src/outdir /pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
--username dyammons@colostate.edu \
--token b451a793ee4a3357367ea2fb3d551fcc \
--single_cell TRUE \
--refsample scrna_celltype_l2.txt \
--mixture mixture.txt \
--rmbatchSmode TRUE 

```
Still working on this....
