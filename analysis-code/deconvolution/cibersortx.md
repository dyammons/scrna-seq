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
In this example we will use a single cell dataset that we have paired samples from bulk seq to assess accuaracy of deconvolution.

#### Get ground truth for paired samples
```r
#Load libraries
library(Seurat)
library(tidyverse)


### Load single cell reference

#Load Seurat object
seu.obj <- readRDS("../../bov_lav_scRNA/output/s3/allCells_clean_viral_S3.rds")
#Select cell type annotation level to use for deconvolution
ct_level <- "celltype.l2"
#Remove `T other` cell subtype and merge monocytes with macrophage
seu.obj <- subset(seu.obj, invert = TRUE, subset = celltype.l2 == "T Other")
Idents(seu.obj) <- "celltype.l2"
seu.obj <- RenameIdents(seu.obj, "Monocyte" = "Macrophage")
seu.obj$celltype.l2 <- droplevels(Idents(seu.obj))


### Extract cell type percentage data to help evaulate accuracy of results

#Extract and filter on samples that have bulk samples
ground_truth <- as.data.frame(
    table(seu.obj@meta.data[[ct_level]], seu.obj$orig.ident)
) %>% filter(
    Var2 %in% c("bov_lav_4651130", "bov_lav_4667180",
                "bov_lav_4630212", "bov_lav_4208263", "bov_lav_4296177")
)
#convert to percentages
ground_truth <- ground_truth %>% 
    group_by(Var2) %>%
    mutate(
        pct = round(Freq / sum(Freq), 2)
    )
```
#### Export the data as a matrix for import into CIBERSORTx

```r
#Down sample to make data easier to work with (keep equal number of cells)
Idents(seu.obj) <- ct_level
seu.obj <- subset(
    x = seu.obj, downsample = min(table(seu.obj@meta.data[[ct_level]]))
)
message(paste0(
    "Each cell type was downsampled to ",
    min(table(seu.obj@meta.data[[ct_level]])),
    " cells.\nConsider this value when selecting replicates ",
    "to build the deconvolution reference."
))
#Get the data - this can use a lot of RAM if there are still a lot of cells
#Here we get raw counts - CIBERSORTx should convert this to CPM
cnt_mat <- FetchData(
    seu.obj, vars = rownames(seu.obj), assay = "RNA", layer = "counts"
)
#Transpose the data to get it in genes (rows) X cells (columns)
cnt_mat <- t(cnt_mat)
#Replace the cell barcodes with the cell type name
colnames(cnt_mat) <- as.character(seu.obj@meta.data[[ct_level]])
cnt_mat <- cbind(rownames(cnt_mat), cnt_mat)
colnames(cnt_mat)[1] <- "GeneSymbol"
cnt_mat <- as.data.frame(cnt_mat)
#Save the matrix
write.table(
    cnt_mat, quote = FALSE, row.names = FALSE, sep = "\t",
    file = paste0("../output/cibersort_data/scrna_", gsub("\\.", "_", ct_level), ".txt")
)
```

### Prepare the mixtures (bulk seq counts) to enusre they are formated properly
This step will vary based on your input, but the approach should provide
some guidance. In this case there is an error in the format of the input data
such that the column numbers vary throughout the file. Although the file is a
comma seperated, read it in as a .tsv, then extract the values needed to build
the dataframe.

```r
#Read input counts data
mixture <- read.table("/pl/active/dow_lab/dylan/bov_lav_scRNA/input_bulk/count_matrix/filtered_counts.csv")
#Get column names
colNames <- unlist(strsplit(split = ",", mixture[1, ]))
colNames <- c("NAME", colNames[(length(colNames)-13):length(colNames)])
#Extract the data values
dat <- strsplit(split = ",", mixture[2:nrow(mixture), ])
clean_dat <- lapply(dat, function(x){
    matrix(
        x[c(1, (length(x)-13):length(x))], 
        ncol = length(colNames), byrow = T
    )
})
clean_dat <- do.call(rbind, clean_dat)
#Bring the column names over and write the file
colnames(clean_dat) <- colNames
write.table(
    clean_dat, quote = FALSE, row.names = FALSE, sep = "\t",
    file = "../output/cibersort_data/mixture.txt"
)
```

### Run CIBERSORTx in CLI via docker container
In a seperate terminal run CIBERSORT using the docker container.
The following is and example Singularity command:

```sh
singularity exec \
-B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_data:/src/data \
-B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_output/noBatch:/src/outdir \
/pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif \
/src/CIBERSORTxFractions \
--username dyammons@colostate.edu \
--token [insert token here] \
--single_cell TRUE \
--refsample scrna_celltype_l2.txt \
--mixture mixture.txt \
--rmbatchBmode TRUE \
--fraction 0 \
--replicates 50 \
--remake TRUE
```
Since we have ground truth values in this example we will run CIBERSORTx three different ways to assess which method has greatest accuracy.
```sh
#Run with no batch correction (need to insert valid token/username)
singularity exec -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_data:/src/data -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_output/noBatch:/src/outdir /pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions --username dyammons@colostate.edu --token  --single_cell TRUE --refsample scrna_celltype_l2.txt --mixture mixture.txt --fraction 0 --replicates 5 --remake TRUE

#Run with bulk mode batch correction (need to insert valid token/username)
singularity exec -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_data:/src/data -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_output/bBatch:/src/outdir /pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions --username dyammons@colostate.edu --token  --single_cell TRUE --refsample scrna_celltype_l2.txt --mixture mixture.txt --fraction 0 --replicates 5 --remake TRUE --rmbatchBmode TRUE

#Run with single cell batch correction, the reccomened correction method (need to insert valid token/username)
singularity exec -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_data:/src/data -B /pl/active/dow_lab/dylan/bov_lav_bulk/output/cibersort_output/sBatch:/src/outdir /pl/active/dow_lab/dylan/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions --username dyammons@colostate.edu --token  --single_cell TRUE --refsample scrna_celltype_l2.txt --mixture mixture.txt --fraction 0 --replicates 5 --remake TRUE --rmbatchSmode TRUE
```


### Evaluate performance
After running the CIBERSORTx command, return to the R session to load in the results. 
If restarting the R session, then you will need to load the `ground_truth` dataframe back in.


### Compare CIBERSORTx batch correction approaches
```r
#Load in the results
filez <- c(
    "../output/cibersort_output/noBatch/CIBERSORTx_Results.txt",
    "../output/cibersort_output/bBatch/CIBERSORTx_Adjusted.txt",
    "../output/cibersort_output/sBatch/CIBERSORTx_Adjusted.txt"
)
df.list <- lapply(filez, function(x){
    res <- read.table(file = x, sep = "\t", header = TRUE)
    res$method <- strsplit(x, "/")[[1]][4]
    return(res)
})
res <- do.call(rbind, df.list)
#Prep data for plotting with ggplot
res <- pivot_longer(res, cols = colnames(res)[2:(ncol(res) - 4)])
res$value <- round(res$value, 2)
res$name <- gsub("\\.", " ", res$name)
res <- left_join(res, meta[ , c(1,3)], by = c("Mixture" = "orig.ident"))
res$method <- factor(res$method, levels = c("noBatch", "bBatch", "sBatch"))
### Plot heatmap of estimated percentages
ggplot(data = res, aes(name, Mixture, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
        low = "white", high = "red", limit = c(0,1), space = "Lab", 
        name = "Estimated cell\nfraction"
    ) +
    theme_minimal() + 
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, 
        size = 12, hjust = 1)
    ) +
    facet_grid(.~method) +
    coord_fixed() + 
    geom_text(aes(name, Mixture, label = value), color = "black", size = 2) +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(paste0("../output/cibersort_output/cell_fractions.png"), width = 16, height = 5)
```

### Plot coorelation with ground truth to assess accuracy
```r
#Correct nomenclature discrepency between sample names
ground_truth$Var2 <- gsub("bov_lav_", "BovLav", ground_truth$Var2)
ground_truth$Var2 <- substr(ground_truth$Var2, 1, nchar(ground_truth$Var2) - 1)
#Filter res to only include the samples that have scrna data
res <- res[res$Mixture %in% c("BovLav465113", "BovLav466718",
                              "BovLav463021", "BovLav420826", 
                              "BovLav429617"), ]
#Join the two dataframes
plotData <- left_join(res, ground_truth, by = c("Mixture" = "Var2", "name" = "Var1"))
plotData$diff <- round(plotData$value - plotData$pct, 2)
#Plot heatmap of percentages diffference
ggplot(data = plotData, aes(name, Mixture, fill = diff)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", 
        limit = c(-1,1), midpoint = 0, space = "Lab", 
        name = "Estimate difference\n(estimated - truth)"
    ) +
    theme_minimal() + 
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, 
        size = 12, hjust = 1)
    ) +
    facet_grid(.~method) +
    coord_fixed() + 
    geom_text(aes(name, Mixture, label = diff), color = "black", size = 2) +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(paste0("../output/cibersort_output/cell_fractions_DIFF.png"), width = 14, height = 5)
#Plot heatmap of estaimtes for samples with percentages diffference
ggplot(data = plotData, aes(name, Mixture, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
        low = "white", high = "red", limit = c(0,1), space = "Lab", 
        name = "Estimated cell\nfraction"
    ) +
    theme_minimal() + 
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, 
        size = 12, hjust = 1)
    ) +
    facet_grid(.~method) +
    coord_fixed() + 
    geom_text(aes(name, Mixture, label = value), color = "black", size = 2) +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(paste0("../output/cibersort_output/cell_fractions_subset.png"), width = 14, height = 5)
#Calcuate the overall difference from "truth"
acc_summary <- plotData %>%
    group_by(method) %>%
    summarize(
        abs_diff = sum(abs(diff))
    )
print(acc_summary)
#Create initial scatter plot
p <- ggplot(plotData, aes(x = pct, y = value)) + 
    stat_smooth(method = "lm", se = FALSE, fullrange = FALSE, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = name, shape = method)) +
    labs(x = "Percentage by scRNA", y = "Percentage by CIBERSORTx") +
    guides(color = guide_legend(title = "Cell type", size = 3, override.aes = list(fill = NA))) +
    theme(
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_line(color = "gray"), 
        panel.grid.minor = element_line(color = "gray"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 14),
        plot.title = element_blank(),
        title = element_text(size = 16),
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
    ) + coord_cartesian(clip = "off") 
ggsave("../output/cibersort_output/scrnaCoorPlot.png", width = 7, height = 5)
#Facet by method and run correlation
pi <- p + facet_wrap("method", scales = "free", nrow = 1) + ggpubr::stat_cor() + theme(plot.background = element_rect(fill = "white"))
ggsave("../output/cibersort_output/scrnaCoorPlot_by_method.png", width = 12, height = 5)
#Facet by method X cell type and run correlation
pi <- p + facet_wrap(method ~ name, scales = "free", nrow = 3) + ggpubr::stat_cor() + theme(plot.background = element_rect(fill = "white"))
ggsave("../output/cibersort_output/scrnaCoorPlot_by_ct.png", width = 14, height = 6, scale = 2)
```

### Conclusions
- In this example B batch correction has the greatest accuracy
- Modification of the annotations stored in `celltype.l2` slot of the Seurat object had to be modified to enhance performance.






