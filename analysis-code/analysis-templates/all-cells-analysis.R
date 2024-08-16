#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

#################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
###  Set variable and load data  ### <<<<<<<<<<<<<<
#################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output param
outName <- ""
exptName <- ""
reduction <- "umap.integrated.harmony"
clusMain <- "clusterID_integrated.harmony"
contrast <- c("", "") #code will test first vs second

#load in preprocessed data
seu.obj <- readRDS("../output/s3/_S3.rds") #modify as needed
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$clusterID_integrated.harmony <- factor(seu.obj$clusterID_integrated.harmony, levels = 0:max(as.numeric(names(table(seu.obj$clusterID_integrated.harmony)))))

################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
###  BEGIN allCells analysis  ### <<<<<<<<<<<<<<
################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, reduction = reduction, nrow = 1, ncol = 3, features = features, 
                    color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste0("../output/", outName, "/", outName, "_QC_feats.png"), width = 9, height = 3)

#generate viln plots using harmony clusters
vilnPlots(seu.obj = seu.obj, groupBy = clusMain, outName = outName,
          outDir = paste0("../output/viln/", outName, "/"), returnViln = T 
         )

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "../output/cb_input/", 
                markers = paste0("../output/viln/", outName, "/", outName, "_", clusMain, "_gene_list.csv"),
                reduction = reduction,  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", clusMain, "name", "cellSource"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                            "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                            "CD4", "MS4A1", "PPBP", "HBM")

)    

#use singleR to ID cells
singleR(seu.obj = seu.obj, clusters = clusMain, reduction = reduction, 
        outDir = "../output/singleR/", outName = outName)

#set colors - run after determining majorIDs -- MANUAL
colArray <- read.csv("./metaData/allCells_ID.csv")
# colArray <- colArray %>% arrange(majorID) %>% mutate(newCol = gg_color_hue(nrow(colArray)*3)[ c( rep(FALSE, 2), TRUE ) ] ) %>% arrange(clusterID_integrated.harmony)
# write.csv(colArray,"./metaData/allCells_ID.csv", row.names = F)


### Plot initial cluster UMAP
pi <- DimPlot(seu.obj, 
                reduction = reduction, 
                group.by = clusMain,
                # cols = colArray$colz, #uncomment after colors added
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = F) 
ggsave(paste0("../output/", outName, "/", outName, "_rawUMAP.png"), width = 7, height = 7)


#after setting majorID, run this code some code to load in colors
colArray.sub <- colArray[colArray$majCol == "yes",]

### Fig supp: umap by major ID
pi <- DimPlot(seu.obj, 
              reduction = reduction, 
              group.by = "majorID",
#               cols = colArray.sub$colz,  #uncomment after colors added
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              repel = TRUE
)
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("../output/", outName, "/", outName, "_majorUMAP.png"), width = 7, height = 7)


### Key feature plots
features <- c("PTPRC","CD3E","CTSW", 
                "DRA","CSF3R","S100A12", 
                "CD68","FLT3","FCER1A", 
                "GPNMB","VEGFB","CD34",
                "COL1A2","MS4A1","TOP2A")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, title.size = 14, features = features, order = F, legJust = "top", reduction = reduction) 
ggsave(paste0("../output/", outName, "/", outName, "_featPlots.png"), width = 9, height = 15)


### Key dot plot features -- this is best with majorID loaded in
p <- majorDot(seu.obj = seu.obj, groupBy = "majorID",
              features = c("ANPEP", "DLA-DRA", "FLT3", "IGHM", "JCHAIN",
                           "MS4A1", "S100A12", "SERPINA1", "CD4", "IL7R", 
                           "CD52", "CCL5", "GZMB", "KLRB1", "CSTB", "IL5RA", 
                           "IL17RB", "GATA3", "TOP2A", "CENPF", "CD34", "CD109")
) + theme(axis.title = element_blank(),
          axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_majorDot.png"), width = 8, height = 6)


### UMAP by sample -- if unequal sample size downsample by cellSource
Idents(seu.obj) <- "orig.ident"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
                reduction = reduction, 
                group.by = "name",
                cols = levels(seu.obj.ds$colz),
                pt.size = 0.25,
                label = FALSE,
                shuffle = TRUE
)
p <- formatUMAP(pi) + labs(colour="Cell source:") + theme(legend.position = "top", 
                                                            legend.direction = "horizontal",
                                                            legend.title=element_text(size=12)
                                                            ) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste0("../output/", outName, "/", outName, "_umap_bySample.png"), width =7, height = 7)


### Stacked bar graph by clusMain
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = clusMain) + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
                    values = levels(seu.obj$colz)) + 
    theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_stackedBar_cluster.png"), width =7, height = 5)


### Stacked bar graph by majorID -- preferred once variable is set
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = "majorID") + scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(labels = levels(seu.obj$name),
                    values = levels(seu.obj$colz)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_stackedBar_major.png"), width =7, height = 5)


### Frequency plots to run stats - cluster
freqy <- freqPlots(seu.obj, method = 1, nrow = 3, 
                   comp = "cellSource", groupBy = clusMain, legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
)
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_cluster.png"), width = 12, height = 8)


### Frequency plots to run stats - major
freqy <- freqPlots(seu.obj, method = 1, nrow = 1, 
                   comp = "cellSource", groupBy = "majorID", legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
)
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_major.png"), width = 8, height = 4)


### Complete linDEG in pseudobulk-type format by all cells -- linDEG is no longer reccomedned
seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparison = "cellSource", 
       outDir = paste0("../output/", outName, "/linDEG/"), outName = outName, colUp = "red", colDwn = "blue",subtitle = F)


### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)


### Or complete linDEG in each major group -- linDEG is no longer reccomedned
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "majorID", comparison = "cellSource", 
       outDir = paste0("../output/", outName, "/linDEG/"), outName = outName, colUp = "red", colDwn = "blue", subtitle = F)
          
############################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
###  END allCells analysis  ### <<<<<<<<<<<<<<
############################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
