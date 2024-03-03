Pending further annotation....

```r
seu.obj <- readRDS(file = "./output/s3/naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
dout <- "./output/copyKat/"

seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6.csv", groupBy = "clusterID", metaAdd = "majorID")

seu.sub.list <- SplitObject(seu.obj, split.by = "orig.ident")

df.list <- list()
cntr <- 0
lapply(seu.sub.list, function(seu.obj) {
    
    rm(copykat_run)
    gc()
    
    #check the distribution of cells that have more than 2000 UMIs
    table(seu.tumor@meta.data$nCount_RNA>2000,seu.tumor$clusterID)
    
    #select best cells
    seu.tumor.clean <- seu.tumor[,seu.tumor@meta.data$nCount_RNA>2000]
    
    #try with "normal" cells called
    norm_cells <- colnames(subset(seu.obj, majorID == "myeloid" | majorID == "tcell"))
    cntr <- cntr+1
    rawData <- as.matrix(seu.obj@assays$RNA@counts)

    outName <- unique(seu.obj$orig.ident)
    
    #run copyKAT
    copykat_run <- copykat(rawmat=rawData, 
                           id.type="S", 
                           ngene.chr=5, 
                           win.size=25, 
                           KS.cut=0.1, 
                           sam.name=outName, 
                           distance="euclidean", 
                           norm.cell.names=norm_cells, 
                           n.cores=12,
                           output.seg="FALSE")

    pred.score <- data.frame(copykat_run$prediction)
    pred.score$cell.names <- NULL 
    CNA.val <- data.frame(copykat_run$CNAmat)

    seu.obj <- AddMetaData(
      object = seu.obj,
      metadata = pred.score,
      col.name = 'copyKatPred'
      )
    
    #df <- as.data.frame(seu.obj@meta.data)
    #df.list[[seu.obj@meta.data]] <- df
    
    pi <- DimPlot(seu.obj, 
                  reduction = "umap", 
                  group.by = "copyKatPred",
                  pt.size = 0.5,
                  label = FALSE,
                  label.box = FALSE
                 )
    
    df <- data.frame(matrix(ncol = 0, nrow = dim(seu.obj)[2]))
    
    df[["code"]] <- rownames(seu.obj@meta.data)
    df[["cnvStat"]] <- seu.obj@meta.data$copyKatPred
    
    df.list[[cntr]] <- df
    
    outfile <- paste(dout, outName, "_copyKatPred_seu.rds", sep = "")
    saveRDS(seu.obj, file = outfile)

    outfile <- paste(dout, outName, "_uMAP_by_ploidy.png", sep = "")
    #save the final plot!
    ggsave(plot = pi, outfile, width = 10, height = 10)
    
})
       
cellCounts <- do.call(rbind, df.list)
outfile <- paste(dout,"cnvStatz.csv", sep = "")
write.csv(cellCounts, file = outfile)

files <- list.files(path = "./output/copyKat_2/", pattern="copyKatPred_seu.rds", all.files=FALSE,
                        full.names=T)

create.seu.call <- function(x) {
        readRDS(x)
    }

#load all seurat objects and store in a large list
katOuts.list <- mapply(create.seu.call, files)

df.list <- list()

for (obj in seq(1:length(katOuts.list))) {
    seu <- katOuts.list[[obj]]
    df <- data.frame(matrix(ncol = 0, nrow = dim(seu)[2]))
    
    
    
    df[["code"]] <- rownames(seu@meta.data)
    df[["cnvStat"]] <- seu@meta.data$copyKatPred
    
    df.list[[obj]] <- df
}
cellCounts <- do.call(rbind, df.list)
outfile <- paste("./output/copyKat_2/cnvStat.csv", sep = "")
write.csv(cellCounts, file = outfile)


### Fig 1F: Plot CopyKAT output
cellCounts <- read.csv("./output/copyKat_2/cnvStat.csv")
cellCounts$X <- NULL
cnts <- cellCounts$cnvStat
names(cnts) <- cellCounts$code

seu.obj <- AddMetaData(
      object = seu.obj,
      metadata = cnts,
      col.name = 'copyKatPred'
      )

#exculde low quality sample
seu.obj.sub <- subset(seu.obj,invert = T,
                  subset = 
                  orig.ident ==  "run_count_tumor_no_tx_7")


pi <- DimPlot(seu.obj.sub, 
                  reduction = "umap", 
                  group.by = "copyKatPred",
                  pt.size = 0.25,
                  label = FALSE,
                  label.box = FALSE,
               shuffle = T
                 )
pi <- formatUMAP(pi) + theme(legend.position = c(0.01, 0.95)) + theme(axis.title = element_blank(),
                                                 panel.border = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_uMAP_by_ploidy.png", sep = ""),width = 7,height=7)
```
