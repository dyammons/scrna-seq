#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### set variables
mainOutName = ""

### create output dirs
dir.create("./output/")
dir.create("./output/s1")

#load in 10x data and qc filter eeach sample
load10x(din = "./input/", dout = "./output/s1/", outName = mainOutName, testQC = T, 
        nFeature_RNA_high = 3000, nFeature_RNA_low = 200, percent.mt_high = 80, 
        nCount_RNA_high = 25000, nCount_RNA_low = 100)

seu.obj <- readRDS("./output/s1/230706_duod_h3c6_NoIntrons_seu.integrated.obj_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = mainOutName, test_dims = c(50,45,40,35,30), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = mainOutName, final.dims = 50, final.res = 0.4, stashID = "clusterID_pre", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.4, n.neighbors = 30, assay = "integrated", saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#generate preliminary figures 
seu.obj <- readRDS("./output/s3/230706_duod_h3c6_NoIntrons_res0.4_dims50_dist0.5_neigh40_S3.rds")
#load in meta data after an intial run through
# colArray <- read.csv("./refColz.csv", header = T)
outName <- "allCells"

### create output dirs
dir.create(paste0("./output/",outName,"/"))
dir.create(paste0("./output/",outName,"/"))
dir.create(c(paste0("./output/",outName,"/"),"./output/viln/","./output/singleR/",paste0("./output/viln/",outName,"/"),paste0("./output/singleR/",outName,"/"))

#check QC params again
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = mainOutName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

#plot inital cluster umap
pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)

#use singleR to use human reference for cell classification
singleR(seu.obj = seu.obj, outName = mainOutName, clusters = "clusterID", outDir = paste0("./output/singleR/",outName,"/"))


# ### Optional: reference map using PBMC data
# reference <- readRDS(file = "../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds")
# reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")

# DefaultAssay(reference) <- "integrated"

# anchors <- FindTransferAnchors(
#     reference = reference,
#     query = seu.obj,
#     normalization.method = "SCT",
#     reference.reduction = "pca", #reference.reduction = "umap",
#     dims= 1:50 #dims= 1:2
# )

# predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l3,
#     dims = 1:50)
# seu.obj <- AddMetaData(seu.obj, metadata = predictions)

# pi <- DimPlot(seu.obj, 
#               reduction = "umap", 
#               group.by = "predicted.id",
#               #cols = levels(seu.obj.ds$colz), #check colorization is correct
#               pt.size = 0.25,
#               label = T,
#               label.box = T,
#               shuffle = F
# )
# pi <- formatUMAP(plot = pi)
# ggsave(paste("./output/", outName, "/", outName, "_umap_Predicted.png", sep = ""), width = 10, height = 7)


### Use conocial markers to ID cells
features <- c("PTPRC","CD3G","CD8A", "CD4",
              "GZMA", "S100A12", "DLA-DRA","FLT3",
              "MS4A1","JCHAIN","TOP2A","GATA3",
               "CD34", "ANPEP","EPCAM", "SI",
              "S100G","RELN", "BCAS1","CHGB"
              )

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 4, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_featPlots.png", sep = ""), width = 12, height = 15)

features <- c("PTPRC","CD68","IL5RA", "LGR5",
              "BMI1", "BMI1", "HOPX","LRIG1",
              "SELL","IL18BP","TERT","CXCL8",
               "CD34", "ESAM","KIT", "SI",
              "CX3CR1","ENOH", "LTF","CHGB"
              )


p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 4, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_featPlots_2.png", sep = ""), width = 12, height = 15)


### After intial run through -- with major IDs added can run:
#plot features
fig1c <- majorDot(seu.obj = seu.obj, groupBy = "majorID",
                  yAxis = c("T cell","Epithelial","Plasma cell","Myeloid","Mast cell","Cycling cell","Fibroblast"),
                  features = c("CD3E", "CCL4", "GZMA","GZMB", "CD8A",
                               "SI", "FABP1", "RBP2", "GUCA2A", "APOA1", 
                               "JCHAIN", "RARRES2", "IGHM", "MS4A1", 
                               "AIF1", "C1QC", "S100A12", "LYZ", "CXCL8", 
                               "KIT", "IGF1", "MS4A2",
                               "TOP2A","MKI67",
                              "ACAT2","TPM2","IGFBP7","SFRP1")
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_majorDot.png", sep = ""), width =8, height = 4)


### Create violin plots for key feats
features <- c(
              "MS4A2", "IL18BP",
              "SELL", "S100A12",
              "DLA-DRA",
              "CCL14", "C1QC",
              "MSR1","CSF1R","CCL3",
              "FLT3", "BATF3", "CADM1","AIF1")

pi <- VlnPlot(
    object = seu.obj,
    pt.size = 0,
    same.y.lims = F,
    group.by = "majorID_sub",
    combine = T,
    cols = majorColors.df$colz,
    stack = T,
    fill.by = "ident",
    flip = T,
    features = features
        ) + NoLegend() + theme(axis.ticks = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank())

#plot <- prettyViln(plot = pi, colorData = NULL, nrow = 2, ncol = 4)
ggsave(paste("./output/", outName, "/", outName, "selectViln.png", sep = ""), width = 5, height =6)

