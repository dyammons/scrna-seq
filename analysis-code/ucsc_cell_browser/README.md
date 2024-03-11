# Steps to generate and host an interactive UCSC cell browser on SharePoint

There are two steps to generating a cell browser:
0. Extract data from the Seurat object(s)
0. Use the cellbrowser module to build the interactive browser


## Extract data from the Seurat object(s)
For this you will need a conda envrionment, singularity container, or other means of working with a Seurat object.


For each object you want to create a browser for, you will need to load in the Seurat object and extract count matricies, metadata, etc....

To demo this, we will pull down data from NCBI GEO from the [canine PBMC atlas dataset](https://github.com/dyammons/Canine_Leukocyte_scRNA) using wget from a linux terminal.

```sh
#pull data down into working directory
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_bcell.rds.gz #b cells
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_cytotoxic.rds.gz #cd8 t, nk, dn t

#unzip the files
gunzip *.rds.gz
```

```r
#load in data and create output directory
seu.obj <- readRDS(file = "../output/s3/GSE225599_bcell.rds")
seu.obj <- UpdateSeuratObject(seu.obj)

outDir <- "./ouput/cb_input/bcell/"
dir.create(outDir)

# #after going through the code below using the bcell object, 
# #load the cytotoxic T cells in a repeat the export steps
# 
# seu.obj <- readRDS(file = "../output/s3/GSE225599_cytotoxic.rds")
# outDir <- "./ouput/cb_input/cytotoxic/"
# dir.create(outDir)
```

```r
Idents(seu.obj) <- "celltype.l3"
cluster.markers <- FindAllMarkers(seu.obj, only.pos = T, min.pct = 0.25)
cluster.markers <- cluster.markers %>% filter(p_val_adj < 0.01)
cluster.markers <- cluster.markers[c("cluster","gene","avg_log2FC","p_val_adj")]
write.table(markers,paste0(outDir,"markers.tsv"), quote=FALSE, sep='\t', row.names = F)
```

```r
meta <- seu.obj@meta.data
meta <- rownames_to_column(meta, "barcode")
meta <- meta[ ,colnames(meta) %in% c("barcode",colsTOkeep)]
write.table(meta, paste0(outDir,"meta.tsv"), quote=FALSE, sep='\t', row.names = F)
```

```r
quickGenes <- as.data.frame(feats[feats %in% rownames(seu.obj)])
colnames(quickGenes) <- "symbol"
write.csv(quickGenes, paste0(outDir, "quickGenes.csv"), quote=FALSE, row.names = F)
```

```r
data.df <- as.data.frame(seu.obj@assays$RNA@layers$data)
colnames(data.df) <- colnames(seu.obj)
rownames(data.df) <- rownames(seu.obj)
data.df <- rownames_to_column(data.df, "symbol")
write.table(data.df, paste0(outDir,"exprMatrix.tsv"), quote=FALSE, sep='\t', row.names = F)
R.utils::gzip(paste0(outDir,"exprMatrix.tsv"), overwrite = TRUE)
```

```r
data.df <- as.data.frame(seu.obj[[reduction]]@cell.embeddings)
data.df <- rownames_to_column(data.df, "barcode")
write.table(data.df,paste0(outDir,reduction,".coords.tsv"), quote=FALSE, sep='\t', row.names = F)
```





Create the env
```sh
conda create -n cb
conda activate cb
```

Install analysis software
```sh
conda install python=3.7
pip install cellbrowser
```




