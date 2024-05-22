#!/usr/bin/Rscript

library(Seurat)
library(tidyverse)
library(clustree)
library(stringr)
# install.packages("remotes")
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(patchwork)
library(scales)
library(cowplot)
library(ggrepel)
library(colorspace)
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
# remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
# BiocManager::install("SingleR")
library(SingleR)
# BiocManager::install("celldex")
library(celldex)
library(viridis)
# library(reshape)
# library(lemon)
# install.packages("devtools")
# devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(msigdbr)
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
# BiocManager::install("slingshot")
library(slingshot)
library(ggpubr)
# BiocManager::install("scRNAseq")
library(scRNAseq)
# BiocManager::install("scuttle")
library(scuttle)
library(ape)
library(lemon)
# BiocManager::install("ggtree")
library(ggtree)
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(reshape2)

############ gg_color_hue ############

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

############ load10x ############
#' Load 10x Cell Ranger count matrices and preprocess data using a Seurat + DoubletFinder protocol
#'
#' @param din string; Path to input directory which should be structured such that each sample has its own directory
#'   with the corresponding features.tsv.gz, matrix.mtx.gz, and barcodes.tsv.gz inside. Use relative path.
#' @param dout string; Path to put figures and processed Seurat objects in. Use relative path.
#' @param outName string; Short name that will be incorporated into the output files
#' @param testQC logical;
#'   * `TRUE` (default): only violin plots with standard QC parameters will be produced.
#'   * `FALSE`: run entire function.
#' @param nFeature_RNA_high numerical; high end cut off for the nFeature_RNA QC parameter.
#' @param nCount_RNA_low numerical; low end cut off for the nFeature_RNA QC parameter.
#' @param nCount_RNA_high numerical; high end cut off for the nCount_RNA QC parameter.
#' @param nCount_RNA_low numerical; low end cut off for the nCount_RNA QC parameter.
#' @param percent.mt_high numerical; high end cut off for the percent.mt QC parameter.
#' @param mt_pattern character; regex to find the mitochondrial features and determine the percentage of reads mapping
#'   to mitochondrial chromosomes. Default: "^MT-".
#' @param mt_feats vector of strings; list of features to use in the calculation of percent.mt
#' @param nfeatures numerical; number of variable features to use during dimension reduction. Default: 2000
#' @param removeDubs logical;
#'   * `TRUE` (default): remove suspected doublets identified using DoubletFinder.
#'   * `FALSE`: label doublets, but do not remove them.
#' @param removeRBC_pal logical; 
#'   * `TRUE`: predict and remove red blood cell and platelet clusters. Red blood cells are predicted based on
#'     the overexpression of HBM and platelets are predicted based on expression of PPBP.
#'   * `FALSE` (default): do not search for RBC or platelet clusters.
#' @param pal_feats vector of strings; list of features to use for estimating the platelet signature. Only used if 
#'   removeRBC_pal == TRUE.
#' @param isolatePalRBC logical; 
#'   * `TRUE`: keep only predicted red blood cell and platelet clusters.
#'   * `FALSE` (default): keep all cells.
#' @param featPlots vector of strings; list of features to plot when making feature plots.
#' 
#' @return no objects returned
#' 
#' @examples
#' \dontrun{
#' # Run initially with testQC set to `TRUE` then use output to set QC thresholds
#' load10x(din = "../input/", dout = "../output/", outName = "experiment1", testQC = TRUE)
#' 
#' # Run with QC thresholds specified
#' load10x(din = "../input/", dout = "../output/", outName = "experiment1", testQC = FALSE,
#'         nFeature_RNA_high = 4500, nFeature_RNA_low = 200,
#'         nCount_RNA_high = 20000, nCount_RNA_low = 100,
#'         percent.mt_high = 10, mt_pattern = "^MT-")
#' }
#' 
#' @export "_S1.png" and "_S1.rds" files into `dout` for downstream usage

load10x <- function(
    din = NULL, 
    dout = NULL, 
    outName = NULL, 
    testQC = TRUE,
    nFeature_RNA_high = 4500, 
    nFeature_RNA_low = 200,
    nCount_RNA_high = 20000, 
    nCount_RNA_low = 100,
    percent.mt_high = 10, 
    mt_pattern = "^MT-",
    mt_feats = NULL,
    nfeatures = 2000,
    removeDubs = TRUE, 
    removeRBC_pal = FALSE,
    pal_feats = NULL, 
    isolatePalRBC = FALSE,
    featPlots = c("PTPRC", "CD3E", "CD8A", "GZMA",
                  "IL7R", "ANPEP", "FLT3", "DLA-DRA",
                  "CD4", "MS4A1", "PPBP","HBM")
    ){

    #if running a real run, update the user on the specified QC parameters that will be applied to the samples
    if (testQC == FALSE){
        message("The QC parameters are: nFeature_RNA < ", nFeature_RNA_high, " & nFeature_RNA > ", nFeature_RNA_low, " & percent.mt < ", percent.mt_high, " & nCount_RNA < ", nCount_RNA_high," & nCount_RNA > ", nCount_RNA_low)
    }

    #get list of sub-dirs in the input dir and store as list to loop through
    fpath <-  paste0("./", din,"/")
    files <- list.files(path = fpath, pattern = NULL, all.files = FALSE,
                        full.names = F)

    df.list <- list()
    for (infile in files) {

        #set up df for export
        if (removeRBC_pal == TRUE){
            df <- data.frame(matrix(ncol = 4, nrow = 0))
            colnames(df) <- c("Initial","Filtered","Platelet_rbc_rm","Singlets")
        } else {
            df <- data.frame(matrix(ncol = 3, nrow = 0))
            colnames(df) <- c("Initial","Filtered","Singlets")
        }

        ### Load in the data

        #set import path
        pwd <- paste0("./", din,"/", infile)

        #read 10X data
        indata <- Read10X(pwd)

        #create seurat object
        seu.obj <- CreateSeuratObject(counts = indata,
                                      project = infile,
                                      min.cells = 3,
                                      min.features = 200)

        #add mitochondrial/hemoglobin/pal QC data to seurat metadata
        if(is.null(mt_feats)){
            seu.obj[["percent.mt"]] <- PercentageFeatureSet(seu.obj, pattern = mt_pattern)
        } else{
            seu.obj[["percent.mt"]] <- PercentageFeatureSet(seu.obj, features = mt_feats)
        }
        seu.obj[["percent.hbm"]] <- PercentageFeatureSet(seu.obj, pattern = "HBM")
        seu.obj[["percent.ppbp"]] <- PercentageFeatureSet(seu.obj, pattern = "PPBP")

        #if list of pal associated features is provided, use the sum of all features to calc "percent.pal"
        if(!is.null(pal_feats)){
            seu.obj[["percent.pal"]] <- PercentageFeatureSet(seu.obj, features = pal_feats)
        }

        #visualize QC metrics as a violin plot
        p <- VlnPlot(seu.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
        ggsave(paste0("./",dout,"/", infile,"_QC_S1.png"))

        #check if running through entire code
        if (testQC == TRUE){
            next
        }

        #stash initial #'s of cells in the seurat object
        df[1,1] <- dim(seu.obj)[2]

        ### Complete dimension reduction

        #subset based on QC cutoffs
        seu.obj <- subset(seu.obj,
                          subset = nFeature_RNA < nFeature_RNA_high & nFeature_RNA > nFeature_RNA_low &
                          percent.mt < percent.mt_high & nCount_RNA < nCount_RNA_high & nCount_RNA > nCount_RNA_low
                         )

        #next steps normalize, scale, and run UMAP
        seu.obj <- NormalizeData(seu.obj,
                                 normalization.method = "LogNormalize",
                                 Scale.factor = 10000)

        seu.obj <- FindVariableFeatures(seu.obj,
                                        selection.method = "vst",
                                        nfeatures = nfeatures)

        seu.obj <- ScaleData(seu.obj) # only scales variable features by default
        seu.obj <- RunPCA(seu.obj) # using the scaled expression of variable features

        seu.obj <- FindNeighbors(seu.obj,
                                 dims = 1:10
                                )

        seu.obj <- FindClusters(seu.obj,
                                resolution = 0.1
                               )

        seu.obj <- RunUMAP(seu.obj,
                           dims = 1:15
                          )

        #visulize featPlots
        p <- FeaturePlot(seu.obj, features = featPlots)
        ggsave(paste0("./", dout, "/", infile, "_featPlotDefault_S1.png"), width = 12, height = 8)

        #visulize UMAP
        p <- DimPlot(seu.obj, reduction = "umap")
        ggsave(paste0("./", dout,"/", infile, "_uMAP_S1.png"))

        #stash initial cluster IDs
        seu.obj[["int.clusID"]] <- Idents(object = seu.obj)

        #stash the number of remaining cells following removal of cells that fail QC
        df[1,2] <- dim(seu.obj)[2]

        ### Optionally remove platelets and red blood cells

        #if wanting to remove platelets and red blood cells, then run this block
        if (removeRBC_pal == TRUE){

            #find the rbc cluster
            if(length(AverageExpression(seu.obj, features = "HBM")) != 0){
                rbc.df <- as.data.frame(AverageExpression(seu.obj, features = "HBM"), header = TRUE)
                rbc.clus <- str_split(as.character(colnames(rbc.df)[max.col(rbc.df,ties.method="first")]),"[.]")[[1]][2]

                sus.rbc.size <- as.numeric(length(WhichCells(seu.obj, expression = HBM > 0)))
                sus.rbc.clus.size <- as.numeric(length(WhichCells(seu.obj, idents = rbc.clus )))

                sus.rbc.comp.pct <- sus.rbc.size/sus.rbc.clus.size
                rbc.clus <- ifelse(sus.rbc.comp.pct > 0.8, rbc.clus, "NULL")
                print(rbc.clus)

                if(is.null(rbc.clus)){
                    message("No red blood cells detected in sample!")
                }else{
                    message("Red blood cells detected in sample!")
                }

            } else{
                rbc.clus <- "NULL"
                message("No red blood cells detected in sample!")
                }

            #find the platelet cluster
            if(length(AverageExpression(seu.obj, features = "PPBP")) != 0){
                pal.df <- as.data.frame(AverageExpression(seu.obj, features = "PPBP"), header = TRUE)
                pal.clus <- str_split(as.character(colnames(pal.df)[max.col(pal.df,ties.method="first")]),"[.]")[[1]][2]

                sus.pal.size <- as.numeric(length(WhichCells(seu.obj, expression = PPBP > 0)))
                sus.pal.clus.size <- as.numeric(length(WhichCells(seu.obj, idents = pal.clus )))

                sus.pal.comp.pct <- sus.pal.size/sus.pal.clus.size
                pal.clus <- ifelse(sus.pal.comp.pct > 0.8, pal.clus, "NULL")
                print(pal.clus)

                if(is.null(pal.clus)){
                    message("No platelets detected in sample!")
                }else{
                    message("Platelets detected in sample!")
                }

            } else{
                pal.clus <- "NULL"
                message("No platelets detected in sample!")
                }

            #check if platelet or rbcs clusters were found, if one or more was found, remove them, else do not remove any clusters
            if(rbc.clus != "NULL" | pal.clus != "NULL"){
                if(isolatePalRBC == F){
                    seu.obj <- subset(seu.obj,
                                      subset = int.clusID != pal.clus & int.clusID != rbc.clus)
                }
                else{
                    seu.obj <- subset(seu.obj,
                                      subset = int.clusID == pal.clus | int.clusID == rbc.clus)
                }
            } else{message("No platelets or red blood cells removed from sample!")}

            #stash the number of remaining cells following removal of platelets and rbcs
            df[1,3] <- dim(seu.obj)[2]
        }

        ### Run DoubletFinder

        #next steps complete doublet identification using DoubletFinder
        #https://github.com/chris-mcginnis-ucsf/DoubletFinder
        sweep.res.list_seu.obj <- paramSweep(seu.obj, PCs = 1:10, sct = FALSE)
        sweep.stats_seu.obj <- summarizeSweep(sweep.res.list_seu.obj, GT = FALSE)
        bcmvn_seu.obj <- find.pK(sweep.stats_seu.obj)
        pk_val <- as.numeric(as.character(bcmvn_seu.obj$pK[bcmvn_seu.obj$BCmetric == max(bcmvn_seu.obj$BCmetric)]))
        annotations <- seu.obj$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)

        # Determine expected doublet rate

        # This formula assumes 0.5% doublet per 1000 cells (lower than 10x 
        # recommendation to account for homotypic doublets and cells removed 
        # during QC filtering)
        expected_dub_rate <- dim(seu.obj)[2]/1000*0.5/100

        nExp_poi <- round(expected_dub_rate*length(seu.obj$orig.ident))
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

        pN_value <- 0.25 #may want to MODIFY value, but this usually works just fine
        pANN_value <- paste0("pANN_",pN_value,"_",pk_val,'_',nExp_poi)

        #run doubletFinder
        seu.obj <- doubletFinder(seu.obj, PCs = 1:10, pN = pN_value, pK = pk_val, 
                                 nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

        #stash the cell calls under the meta.data slot "doublet"
        seu.obj[["doublet"]] <- seu.obj[[paste0("DF.classifications_", pN_value, "_", pk_val, '_', nExp_poi)]]

        #export UMAP highlighting doublet vs singlet cells
        DimPlot(seu.obj, pt.size = 1, label = TRUE, label.size = 5, reduction = "umap", group.by = "doublet" )
        ggsave(paste0("./", dout, "/", infile,"_DF_S1.png") )

        if (removeDubs == TRUE){
            #remove putative doublets
            seu.obj <- subset(seu.obj,
                              subset = doublet == "Singlet"
                             )
        }

        #stash the number of remaining cells following removal of doublets
        if (removeRBC_pal == TRUE){
            df[1,4] <- dim(seu.obj)[2]
        } else {df[1,3] <- dim(seu.obj)[2]
               }

        #stash the values in data frame list to be rbind later
        rownames(df) <- infile
        df.list[[which(infile == files)]] <- df

        ### Complete dimension reduction again to allow for visual inspection before proceeding 

        #next steps normalize, scale, and run UMAP
        seu.obj <- NormalizeData(seu.obj,
                                 normalization.method = "LogNormalize",
                                 Scale.factor = 10000)

        seu.obj <- FindVariableFeatures(seu.obj,
                                        selection.method = "vst",
                                        nfeatures = nfeatures)

        p <- VariableFeaturePlot(seu.obj)
        ggsave(paste0("./",dout,"/", infile,"_varFeatPlot_post_doubletFinder_S1.png"), 
               width = 7, height = 7)

        seu.obj <- ScaleData(seu.obj)

        #recluster the data with all unwanted cells removed
        seu.obj <- RunPCA(seu.obj, features = VariableFeatures(object = seu.obj))

        seu.obj <- FindNeighbors(seu.obj,
                                 dims = 1:10
                                )

        seu.obj <- FindClusters(seu.obj,
                                resolution = 0.1
                               )

        seu.obj <- RunUMAP(seu.obj,
                           dims = 1:15
                          )

        #plot features
        p <- FeaturePlot(seu.obj,features = featPlots)
        ggsave(paste0("./", dout, "/", infile, 
                      "_featPlotDefault_post_doubletFinder_S1.png"), width = 12, height = 8)

        #plot final umap
        DimPlot(seu.obj, reduction = "umap")
        ggsave(paste0("./",dout, "/", infile, "_uMAP_post_doubletFinder_S1.png") )

        #identify cycling cells
        seu.obj <- CellCycleScoring(
            object = seu.obj,
            s.features = cc.genes.updated.2019$s.genes,
            g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE
        )

        #export UMAP with cycling data
        DimPlot(seu.obj, reduction = "umap", group.by = "Phase")
        ggsave(paste0("./",dout,"/", infile,"_uMAP_cellCylce_post_doubletFinder_S1.png"))

        #store cell cycle state
        seu.obj[["clusters"]] <- Idents(object = seu.obj)

        #export processed seurat object as an .RDS file
        saveRDS(seu.obj, file = paste0("./",dout,"/", infile,"_S1.rds"))
    }

    #save the cell numbers throughout filtering
    cellCounts <- do.call(rbind, df.list)
    write.csv(cellCounts, file = paste0("./",dout,"/", outName, "_cell_counts_S1.csv"))
}

############ integrateData ############
#' Integrate multiple datasets into one using Seurat workflow
#'
#' @param din string; Path to directory containing the *_S1.rds files from the load10x() function. Use relative path.
#' @param pattern string; Suffix that is used to search for file in din to load.
#' @param dout string; Path to put figures and processed Seurat objects. Use relative path.
#' @param outName string; Short name that will be incorporated into the output files.
#' @param saveRDS logical; 
#'   * `TRUE` (default): save an integrated RDS file.
#'   * `FALSE`: Do not save an integrated RDS file.
#' @param normalization.method string; Options are "LogNormalize" or "SCT". This is the normalization approach used 
#'   prior to data integration. Default is "LogNormalize".
#' @param method string; Options are "CCAIntegration", "HarmonyIntegration", "JointPCAIntegration", or 
#'   "RPCAIntegration". The integrated data layer will be saved in the "integrated" slot. Value not used if 
#'   `runAllMethods` is set to `TRUE`. Default is "CCAIntegration".
#' @param runAllMethods logical; 
#'   * `TRUE`: run all the of the 4 main Seurat integration methods and stash the reductions into one object. Results of
#'     CCAIntegration = integrated.cca, HarmonyIntegration = integrated.harmony, JointPCAIntegration = integrated.joint,
#'     and RPCAIntegration = integrated.rcpa.
#'   * `FALSE` (default): Only run the 1 method specified in `method` argument.
#' @param indReClus logical; 
#'   * `TRUE`: complete integration on a previously integrated Seurat object. Must use the `seu.obj` argument to use 
#'     this option.
#'   * `FALSE` (default): intergate seperate objects into one.
#' @param seu.obj variable; Seurat object that is presubset on the cell type to compelte independent reclustering on. 
#'   Will split object by orig.ident metadata slot then reintegrate.
#' @param orig.reduction string; typically want pca, see Seurat documentation for IntegrateLayers().
#' @param read_h5 logical; 
#'   * `TRUE`: input should be readable using the Read10X_h5() function. 
#'   * `FALSE`: load data from previously processed .rds files.
#'   TO DO: This is still in testing and will likely be moved to load10x() function.
#'
#' @return An integrated Seurat object
#' 
#' @examples
#' \dontrun{
#' # Run integration using all 4 methods on the output of load10x()
#' seu.obj <- integrateData(din = "../output/s1/", pattern = "_S1.rds",
#'                          saveRDS = FALSE, runAllMethods = TRUE)
#' 
#' # Run integration on a previously integrated object that was subset on T cells
#' seu.obj.sub <- subset(seu.obj, idents = "tcells")
#' seu.obj <- integrateData(seu.obj = seu.obj.sub, indReClus = T,
#'                          saveRDS = FALSE, runAllMethods = TRUE)
#' }
#' 
#' @export optionally save "_S2.rds" file
#' @export save elbow plot under "_integrated_S2_elbow.png" in dout


integrateData <- function(
    din = "../output/s1/", 
    pattern = "_S1.rds",
    saveRDS = F, 
    outName = "",  
    dout = "../output/s2/",
    orig.reduction = "pca",
    normalization.method = "LogNormalize", 
    method = "CCAIntegration",
    new.reduction.name = "integrated",
    runAllMethods = FALSE,
    indReClus = FALSE, 
    seu.obj = NULL,
    read_h5 = FALSE,
    k = 100, 
    min.cell = 30,
    ...
    ) {
    
    options(future.globals.maxSize = 3e+09)
    
    if(!indReClus){
        #get seurat objects to process and integrate
        fpath <- paste0(din,"/") 
        files <- list.files(path = fpath, pattern = pattern, all.files=FALSE,
                            full.names=FALSE)
        sampleNames <-  unlist(lapply(files, function(x){strsplit(x, pattern)[[1]][1]}))
        
        if(read_h5){
            seu.list <- lapply(files, function(inFile){
                inFile_pwd <- paste0(fpath, inFile)
                projName <- paste(strsplit(inFile, "_")[[1]][1:3], collapse = "_")
                spar.matrix <- Read10X_h5(inFile_pwd, use.names = TRUE, unique.features = TRUE)
                seu.obj <- CreateSeuratObject(spar.matrix, project = projName)
                return(seu.obj)
            })
        } else {
            files <- paste0(fpath, files)
            seu.list <- lapply(files, readRDS)        
        }
        
        #merge object list into one Seurat object with mulitple layers
        seu.obj <- merge(seu.list[1][[1]], y = seu.list[2:length(seu.list)],
                          add.cell.ids = sampleNames, 
                          project = outName
                         )
    } else{
        #check input
        if(is.null(seu.obj)){
            message("To run indReClus == TRUE you must provide an integrated Seurat object for the seu.obj argument.")
            break()
        }
        #split objects to allow for re-integration
        try(seu.obj[['RNA']] <- split(seu.obj[["RNA"]], f = seu.obj$orig.ident), silent = T)
    }

    #if any sample has less than 30 cells. exlude the sample from downstream analysis
    if(min(table(seu.obj$orig.ident)) < min.cell){
        exclude <- as.data.frame(table(seu.obj$orig.ident)) %>% filter(Freq < min.cell) %>% pull(Var1)
        seu.obj <- subset(seu.obj, invert = T, 
                          subset = orig.ident %in% exclude)
        seu.obj$orig.ident <- factor(seu.obj$orig.ident)
    }

    #if any sample has less than 100 cells. modify k.weight to match the sample with the lowest number of cells
    if(k > min(table(seu.obj$orig.ident))){
        k <- min(table(seu.obj$orig.ident))-5
    }

    
    if(normalization.method == "LogNormalize"){
        seu.obj <- seu.obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% 
        RunPCA() %>% FindNeighbors(., dims = 1:30, reduction = "pca") %>% 
        FindClusters(., resolution = 1, cluster.name = "unintegrated_clusters")
    } else if(normalization.method == "SCT"){
        seu.obj <- SCTransform(seu.obj)
        seu.obj <- RunPCA(seu.obj) %>% FindNeighbors(., dims = 1:30, reduction = "pca") %>% 
        FindClusters(., resolution = 1, cluster.name = "unintegrated_clusters")
    } else{
        message("The normalization.method value is not recognized. The options are 'LogNormalize' or 'SCT'. Please update argument and try again.")
        break()
    }

    if(!runAllMethods){
        seu.obj <- IntegrateLayers(object = seu.obj, 
                                   method = noquote(method), 
                                   k.weight = k,
                                   orig.reduction = "pca", 
                                   new.reduction = new.reduction.name,
                                   normalization.method = normalization.method,
                                   verbose = TRUE,
                                   ...
                                  )
    
        # re-join layers after integration
        seu.obj[["RNA"]] <- JoinLayers(seu.obj[["RNA"]])
        seu.obj <- RunPCA(seu.obj)
        gc()
        
        } else{
        
        seu.obj <- IntegrateLayers(object = seu.obj, 
                                   method = CCAIntegration, 
                                   k.weight = k,
                                   orig.reduction = "pca", 
                                   new.reduction = "integrated.cca",
                                   normalization.method = normalization.method,
                                   verbose = FALSE)
        gc()
        seu.obj <- IntegrateLayers(object = seu.obj, 
                                   method = HarmonyIntegration, 
                                   k.weight = k,
                                   orig.reduction = "pca", 
                                   new.reduction = "integrated.harmony",
                                   normalization.method = normalization.method,
                                   verbose = FALSE)
        gc()
        seu.obj <- IntegrateLayers(object = seu.obj, 
                                   method = JointPCAIntegration, 
                                   k.weight = k,
                                   orig.reduction = "pca", 
                                   new.reduction = "integrated.joint",
                                   normalization.method = normalization.method,
                                   verbose = FALSE)
        gc()
        seu.obj <- IntegrateLayers(object = seu.obj, 
                                   method = RPCAIntegration,
                                   k.weight = k,
                                   orig.reduction = "pca", 
                                   new.reduction = "integrated.rcpa",
                                   normalization.method = normalization.method,
                                   verbose = FALSE)
        gc()
        
        # re-join layers after integration
        seu.obj[["RNA"]] <- JoinLayers(seu.obj[["RNA"]])
    }
    
    #run PCA to start the dimension reduction process
    p <- ElbowPlot(seu.obj, ndims = 50)
    ggsave(paste0(dout, outName, "_integrated_S2_elbow.png"))

    if(saveRDS){
        saveRDS(seu.obj, file = paste0(dout, outName, "_integrated_S2.rds"))
    }
    
    return(seu.obj)
}

############ integrateData_v4 ############
#' Integrate mulitple datasets into one using Seurat v3/v4 SCT+CCA integration workflow
#'
#' @param seu.list List of Seurat objects to integrate
#' @param subName Short name that will be incorporated into the output files
#' @param outDir the desired path to put figures and processed Seurat objects. Use realtive path
#' @param featTOexclude If provided these features will be forced to be excluded from integration and dimension reduction
#' @param varFeatsToUse If provided only these features used to complete integration and dimension reduction
#' @param nfeatures Number of features to use for integration. Default is Seurat's default of 2000
#' @param k Value passed to k.weight. Default is 100. This should only be used if working with a(n) object(s) with < 100 cells in at least one of the samples
#' @param saveRDS Logical; If TRUE then a .rds file will be saved in outDir using subName in the name
#' @param returnObj Logical; If TRUE then an integrated Seurat object will be returned
#' @param ndims Numeric; Number of dimensions to use for PCA. Default is 50
#' @param vars.to.regress String; Variable(s) to regress during integration. The provided value must be stored in the metadata of all the Seurat objects to be integrated. Default is "percent.mt"
#' @param z Numerical; Value passed to k.filter, k.score, and dims if integrating a(n) object(s) with < 200 cells in at least one of the samples. Manually modify the value to a number equal to the object with the fewest number of cells
#' 
#' @return An integrated Seurat object (if requested using returnObj)
#' @examples 
#' @export optionally save "_S2.rds" file
#' @export save elbow plot under "_integrated_S2_elbow.png" file
integrateData_v4 <- function(
    seu.list = NULL, 
    subName = "", 
    outDir = "", 
    featTOexclude = NULL, 
    varFeatsToUse = NULL, 
    nfeatures = 2000, 
    k = NULL, 
    saveRDS = F, 
    returnObj = T,
    ndims = 50,
    vars.to.regress = "percent.mt", 
    z = 201
) {
    
    message(paste0("The subclustered object will be output as: ", outDir, subName,"_S2.rds", sep = ""))
    
    options(future.globals.maxSize = 100000 * 1024^2)
    
    #read in data
    seu.sub.list <- seu.list
    k <- ifelse(is.null(k),100,k)
    
    if(is.null(varFeatsToUse)){
        seu.obj <- lapply(seu.sub.list,
                          SCTransform, 
                          vars.to.regress = vars.to.regress,
                          verbose = FALSE,
                          conserve.memory=TRUE)
    
        SelectedFeatures <- SelectIntegrationFeatures(object.list = seu.obj,
                                                      nfeatures = nfeatures) 
    
        if(!is.null(featTOexclude)){
            SelectedFeatures <- SelectedFeatures[!SelectedFeatures %in% featTOexclude]
            if(nfeatures != length(SelectedFeatures)){
                message <- paste0("NOTE: ", featTOexclude, 
                                  " was/were excluded from the variable features used in integration!")
                print(message)
                SelectedFeatures <- SelectIntegrationFeatures(object.list = seu.obj,
                                                              nfeatures = nfeatures+(nfeatures-length(SelectedFeatures))
                                                             )
                SelectedFeatures <- SelectedFeatures[!SelectedFeatures %in% featTOexclude]
            }else{
                message(paste0("NOTE: The features to exclude (", featTOexclude, 
                               ") was/were not included in the variable features used in integration, so the option was not used."))
            }
        }
        
        seu.integrated <- PrepSCTIntegration(object.list = seu.obj,
                                            anchor.features = SelectedFeatures,
                                            verbose = FALSE
                                            )
    }else{
        seu.obj <- lapply(seu.sub.list,
                          SCTransform, 
                          vars.to.regress = vars.to.regress,
                          verbose = FALSE,
                          variable.features.n = nfeatures,
                          return.only.var.genes = FALSE)
        
        seu.obj <- lapply(seu.obj, function(obj){
            varKeep <- varFeatsToUse[varFeatsToUse %in% rownames(obj@assays$SCT@scale.data)]
            obj@assays$SCT@scale.data = obj@assays$SCT@scale.data[varKeep, ]
            obj@assays$SCT@var.features = varFeatsToUse
            return(obj)
        })
        
         varFeats <- lapply(seu.obj, function(obj){
             varKeep <- varFeatsToUse[varFeatsToUse %in% rownames(obj@assays$SCT@scale.data)]
             return(varKeep)
         })
        
        SelectedFeatures <- Reduce(intersect, varFeats)
        seu.integrated <- PrepSCTIntegration(object.list = seu.obj,
                                             anchor.features = SelectedFeatures,
                                             verbose = FALSE
                                             )
        
    }

    gc()

    seu.integrated.anchors <- FindIntegrationAnchors(
        object.list = seu.integrated,
        normalization.method = "SCT",
        anchor.features = SelectedFeatures,
        dims = 1:ifelse(min(z)>30, 30, min(z)-1),
        k.filter = ifelse(min(z)>200, 200, min(z)-1),
        k.score = ifelse(min(z)>30, 30, min(z)-1)
    )

    #clean up environment a bit
    rm(seu.sub)
    rm(seu.sub.list)
    rm(seu.obj)
    rm(seu.integrated)
    gc()

    #integrate data and keep full gene set - still might not be retaining all genes
    seu.integrated.obj <- IntegrateData(
        anchorset = seu.integrated.anchors,
        normalization.method = "SCT",
        k.weight = k,
        verbose = FALSE
    )

    #clean up environment a bit
    rm(seu.integrated.anchors)
    gc()

    seu.integrated.obj <- RunPCA(seu.integrated.obj)

    outfile <- paste(outDir, subName,"_S2_elbow.png", sep = "")
    p <- ElbowPlot(seu.integrated.obj, ndims = ndims)
    ggsave(outfile)

    DefaultAssay(seu.integrated.obj) <- "integrated"

    if(saveRDS){
        outfile <- paste(outDir, subName,"_S2.rds", sep = "")
        saveRDS(seu.integrated.obj, file = outfile)
    }
    
    if(returnObj){
        return(seu.integrated.obj)
    }
}


############ clusTree ############
#' Cluster tree visualiziation to identify ideal clustering resolution and cluster stability
#'
#' @param seu.obj variable; Seurat object to plot.
#' @param dout string; Path to put figures and processed Seurat objects. Use relative path.
#' @param outName string; Short name that will be incorporated into the output files.
#' @param test.dims vector of numericals; Dimensions to run the visualization on.
#' @param resolution vector of numericals; Resolution to run the visualization on.
#' @param algorithm integer; Seurat clustering algorithm for modularity optimization. Options are 
#'   * 1 = original Louvain algorithm
#'   * 2 = Louvain algorithm with multilevel refinement
#'   * 3 (default) = SLM algorithm 
#'   * 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param prefix sting; prefix to use. Typically need to be set to "RNA_snn_res.".
#'
#' @return Nothing
#' 
#' @examples
#' \dontrun{
#' # Run clustree for dimensions set to 45 and 50
#' seu.obj <- clusTree(seu.obj = seu.obj, test_dims = c(45, 50))
#' }
#' 
#' #' @export "[dout]_clustree_[outName].png"

clusTree <- function(
    seu.obj = NULL,
    dout = "../output/clustree/", 
    outName = "",
    test_dims = c(50,45,40,35),
    resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)),
    algorithm = 3,
    prefix = "RNA_snn_res.",
    reduction = "integrated.harmony"
){

    for (dimz in test_dims){
        seu.test <- FindNeighbors(object = seu.obj, reduction = reduction, dims = 1:dimz)
        seu.test <- FindClusters(object = seu.test, algorithm = algorithm, resolution = resolution)
        
        p <- clustree::clustree(seu.test, prefix = prefix) + ggtitle(paste("The number of dims used:", dimz))

        png(file = paste0(dout, dimz , "_clustree_", outName, ".png") , height = 1100, width = 2000)
        print(p)
        dev.off()
    }
}

############ dataVisUMAP ############
#' Complete dimension reduction and visualize data
#'
#' @param seu.obj variable; Seurat object to process.
#' @param dout string; Path to put figures and processed Seurat objects. Use relative path.
#' @param outName string; Short name that will be incorporated into the output files.
#' @param final.dims numerical; Number of dimensions (PCs) to use for reduction. Typically between 15 and 50; default 40.
#' @param final.res numerical; Clustering resolution. Typically between 0.1 and 2; default 0.8.
#' @param min.dist numerical; Minmunm neighbor distance. Typically between 0.1 and 0.8; default 0.3.
#' @param n.neighbors numerical; Number of neighbors to use. Typically between 10 and 75; default 30.
#' @param stashID string; Name of metadata slot to store the results of unsupervised clustering. Will override existing
#'   metadata if the slot name already exists.
#' @param algorithm integer; Seurat clustering algorithm for modularity optimization. Options are 
#'   * 1 = original Louvain algorithm
#'   * 2 = Louvain algorithm with multilevel refinement
#'   * 3 (default) = SLM algorithm 
#'   * 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param prefix string; Typically set to "RNA_snn_res.".
#' @param assay string; Generally will want to use "RNA" assay.
#' @param reduction string; A valid reduction in the reductions slot.
#' @param saveRDS logical; if TRUE then a Seurat object will be saved with a "_S3.rds" suffix in dout
#'   * `TRUE` (default): saves the dimension reeduced Seurat object as "_S3.rds" in dout.
#'   * `FALSE`: does not save a .rds file.
#' @param features vector of strings; Optional list of features to plot after dimension reduction.
#'
#' @return A Seurat object with data in the UMAP reductions slot
#' 
#' @examples
#' \dontrun{
#' # Run clustree for dimensions set to 45 and 50
#' seu.obj <- dataVisUMAP(seu.obj = seu.obj, dout = "../output/s3/",
#'                        outName = "demoData", final.dims = 45,
#'                        final.res = 0.8, min.dist = 0.3,
#'                        n.neighbors = 30)
#' }
#' 
#' #' @export "_cluster_S3.png" "_sample_S3.png" "_featPlotDefault_S3.png" "_S3.rds"

dataVisUMAP <- function(file = NULL, seu.obj = NULL, 
                        outDir = "", outName = "", 
                        final.dims = NULL, final.res = NULL, min.dist = 0.6, n.neighbors = 75,
                        stashID = "clusterID", 
                        returnFeats = T,
                        algorithm = 3, 
                        prefix = "integrated_snn_res.",  
                        assay = "integrated", reduction = "RNA",
                        saveRDS = T, return_obj = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       ) {
    
    #can read in via .rds file if desired
    if(!is.null(file)){
        try(seu.integrated.obj <- readRDS(file), silent = T)
    }else if(!is.null(seu.obj)){
        seu.integrated.obj <- seu.obj
    }

    #perform clustering: Neighbor finding and resolution parameters
    DefaultAssay(seu.integrated.obj) <- assay
    seu.integrated.obj <- FindNeighbors(object = seu.integrated.obj, reduction = reduction, dims = 1:final.dims)
    seu.integrated.obj <- FindClusters(object = seu.integrated.obj, algorithm = algorithm, resolution = final.res)

    #choose appropriate clustering resolution
    res <- paste0(prefix, final.res) 
    seu.integrated.obj <- SetIdent(object = seu.integrated.obj, value = res)
    reduction.name <- paste0("umap.", reduction)
    #run UMAP
    seu.integrated.obj <- RunUMAP(seu.integrated.obj,
                                  reduction = reduction,
                                  reduction.name = reduction.name,
                                  dims = 1:final.dims,
                                  min.dist = min.dist,
                                  n.neighbors = n.neighbors
                                 )

    #save cluster number as "clusterID"
    stashID <- paste0(stashID, "_", reduction)
    seu.integrated.obj[[stashID]] <- seu.integrated.obj@active.ident

    #build hierarchical clustering based on clustering results
    DefaultAssay(seu.integrated.obj) <- "RNA"
    seu.integrated.obj <- NormalizeData(seu.integrated.obj)
    seu.integrated.obj <- BuildClusterTree(seu.integrated.obj, assay = "RNA", dims = 1:final.dims)

    #export UMAP colored by sample ID and cluster
    p <- DimPlot(seu.integrated.obj, label = TRUE, reduction = reduction.name, group.by = stashID)
    ggsave(paste0(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist",min.dist, "_neigh",n.neighbors, "_cluster_S3.png"))
    
    p <- DimPlot(seu.integrated.obj, reduction = reduction.name, group.by = "orig.ident")
    ggsave(paste0(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist",min.dist, "_neigh",n.neighbors, "_sample_S3.png"))

    if(returnFeats == T){
     
        #visulize UMAP featPlots
        p <- FeaturePlot(seu.integrated.obj, reduction = reduction.name, features = features)
        ggsave(paste0(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist", min.dist, "_neigh", n.neighbors, "_featPlotDefault_S3.png"), width = 12, height = 8)

    }

    #save seurat object as rds
    if(saveRDS == T){
        saveRDS(seu.integrated.obj, file = paste0(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist", min.dist, "_neigh", n.neighbors, "_S3.rds"))
    }
    return(seu.integrated.obj)
}


############ prettyFeats ############
#' Custamizable feature plots. Builds off of FeaturePlots() Seurat function. 
#'
#' @param seu.obj variable; Seurat object to plot data from.
#' @param features string; List of features to plot. Must be a gene symbol in the Seurat object or a metadata slot.
#' @param reduction string; A valid reduction in the reductions slot.
#' @param nrow integer; Number of rows of plots to create.
#' @param ncol integer; Number of columns of plots to create.
#' @param titles string; List of titles to use of each feature plot. Ensure order matches `features` order.
#' @param title.size integer; font size to use for plot tiles.
#' @param color string; List of font colors to use for the titles of each feature plot. Ensure order matches `features` 
#'   order.
#' @param noLegend logical;
#'   * `TRUE`: removes the legend from the plot. Plotting is done on variable scales.
#'   * `FALSE` (default): keeps the legend. Use `legJust`, `bottomLeg`, and `legInLine` to further custamize the legend
#' @param legJust string; Option "bottom" or NULL. If NULL will be plotted in the top row of plots.
#' @param bottomLeg logical; 
#'   * `TRUE`: plots legend horizontally below the feature plots.
#'   * `FALSE` (default): plots legend based on other paramters.
#' @param legInLine logical;
#'   * `TRUE`: plots the legend in a plot location.
#'   * `FALSE` (default): plots legend outside of the feature plots.
#' @param showAxis logical;
#'   * `TRUE`: plots axis that spans all plots.
#'   * `FALSE` (default): no axis is plotted.
#' @param smallAxis logical;
#'   * `TRUE`: plots small axis that spans half of the plot in the bottom left.
#'   * `FALSE` (default): no axis is plotted.
#' @param order logical; passed to Seurat FeaturePlot()
#'   * `TRUE`: plots non-zero values on top of 0 values.
#'   * `FALSE` (default): randomly plots data points.
#' @param min.cutoff numeric; Specify minimum expression theshold to be colored. Passed to Seurat FeaturePlot()
#' @param pt.size numeric; Specify size of dots in plot. Passed to Seurat FeaturePlot()
#'
#' @return A patchwork object with the custamized feature plots.
#' 
#' @examples
#' \dontrun{
#' p <- prettyFeats(seu.obj = seu.obj, 
#'                  reduction = "umap.integrated.harmony",
#'                  nrow = 3, ncol = 4, noLegend = T,
#'                  features = c("PTPRC", "CD3E", "CD8A", "GZMA",
#'                               "IL7R", "ANPEP", "FLT3", "DLA-DRA",
#'                               "CD4", "MS4A1", "PPBP","HBM"))
#' }
#' 
#' @export Nothing
#' 

prettyFeats <- function(
    seu.obj = NULL, 
    features = "",
    reduction = "umap.integrated",
    nrow = 3, 
    ncol = NULL, 
    titles = NULL, 
    title.size = 18, 
    color = "black", 
    noLegend = FALSE, 
    legJust = "bottom",
    bottomLeg = FALSE, 
    legInLine = FALSE,
    showAxis = FALSE,
    smallAxis = FALSE,
    order = FALSE, 
    min.cutoff = NA, 
    pt.size = NULL,
    ...
){

    DefaultAssay(seu.obj) <- "RNA"
    features <- features[features %in% c(unlist(rownames(seu.obj)),unlist(colnames(seu.obj@meta.data)))]

    if(is.null(ncol)){
        ncol = ceiling(sqrt(length(features)))
    }

    if(is.null(titles)){
        titles <- features #- add if statement
    }


    #strip the plots of axis and modify titles and legend -- store as large list
    plots <- Map(function(x,y,z) FeaturePlot(seu.obj,features = x, reduction = reduction, pt.size = pt.size, order = order, min.cutoff = min.cutoff, ...) + labs(x = "UMAP1", y = "UMAP2") +
                 theme(axis.text= element_blank(),
                       axis.ticks = element_blank(),
                       axis.title = element_blank(),
                       axis.line = element_blank(),
                       title = element_text(size= title.size, colour = y),
                       legend.position = "none"
                      ) +
                 scale_color_gradient(breaks = pretty_breaks(n = 3), limits = c(NA, NA), low = "lightgrey", high = "darkblue") +
                 ggtitle(z), x = features, y = color, z = titles)


    asses <- ggplot() + labs(x = "UMAP1", y = "UMAP2") +
    theme(axis.line = element_line(colour = "black",
                                   arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                                 ends = "last", type = "closed"),
                                  ),
          axis.title.y = element_text(colour = "black", size = 20),
          axis.title.x = element_text(colour = "black", size = 20),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
         )
    if(!noLegend){
        leg <- FeaturePlot(seu.obj,features = features[1], pt.size = 0.1) +
        theme(legend.position = 'bottom',
              legend.direction = 'vertical',
              legend.justification = "center",
              panel.border = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
             ) +
        scale_color_gradient(breaks = pretty_breaks(n = 1), labels = c("low", "high"), limits = c(0,1), low = "lightgrey", high = "darkblue") +
        guides(color = guide_colourbar(barwidth = 1))

        if(bottomLeg){
            leg <- leg + theme(legend.direction = 'horizontal') + guides(color = guide_colourbar(barwidth = 8))
        }

        legg <- get_legend(leg)
    }

    #nrow <- ceiling(length(plots)/ncol) - add if statement
    patch <- area()

    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
            counter = counter+1
            if (counter <= length(plots)) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
            }
        }
    }

    if(!noLegend){
        if(!bottomLeg & !legInLine){
            legPos <- ifelse(legJust == "bottom",ceiling(length(features)/ncol),1)
            patch <- append(patch, area(t = legPos, l = ncol+1, b = legPos, r = ncol+1))
        }else if(legInLine){
            patch <- append(patch, area(t = nrow, l = ncol, b = nrow, r = ncol))
        }else{
            patch <- append(patch, area(t = ceiling(length(features)/ncol)+1, l = ncol, b = ceiling(length(features)/ncol)+1, r = ncol))
        }
    }else{
        if(smallAxis){
            patch <- append(patch, area(t = nrow, l = 1, b = nrow, r = 1))
        }else{
            patch <- append(patch, area(t = 1, l = 1, b = nrow, r = ncol))
        }
    }

    p <- Reduce( `+`, plots ) +  {if(noLegend & showAxis){asses}} +
    {if(!noLegend){legg}} + plot_layout(guides = "collect") +
    {if(!noLegend & bottomLeg){plot_layout(design = patch, heights = c(rep.int(1, nrow),0.2))}else if(!noLegend & !bottomLeg){plot_layout(design = patch, widths = c(rep.int(1, ncol),0.2))}else{plot_layout(design = patch, widths = rep.int(1, ncol))}}

    return(p)
}


############ linDEG ############

#features to add
#error if more than 2 levels ORR make if so user specifies compatisipon

linDEG <- function(
    seu.obj = NULL, 
    threshold = 1, 
    thresLine = F, 
    groupBy = "clusterID",
    comparison = "cellSource", 
    outDir = "./output/", 
    outName = "", 
    cluster = NULL,
    labCutoff = 20, 
    noTitle = F, 
    contrast = NULL, 
    colUp = "red", 
    colDwn = "blue",
    subtitle = T, 
    returnUpList = F, 
    returnDwnList = F, 
    forceReturn = F, 
    useLineThreshold = F,
    pValCutoff = 0.01, 
    logfc.threshold = 0.58,
    saveGeneList = F, 
    addLabs = "", 
    labsHide = c("^MT-","^RPL","^ENS","^RPS"), 
    returnPlots = F,
    ...
    ){
    
    #set active.ident
    Idents(seu.obj) <- groupBy    
    
    if(is.null(contrast)){
        message(paste0("It is now reccomened that a contrast be provided. This should be two levels present in the comparison metadata slot.",
        "Please provide this information as c(",  "'comp1'", ",",  "'comp2'",") and the analysis will run comp1 VS comp2."))
        contrast <- levels(seu.obj@meta.data[[comparison]])[c(1:2)]
        message(paste0("Contrasting - (1) ", contrast[1]," - VS - (2) ", contrast[2]))
    }
    
    #loop through active.ident to make plots for each group
    lapply(if(is.null(cluster)){levels(seu.obj)}else{cluster}, function(x) {

        #subset on population of cells to analyze -- if all cells are included this will effectively not do anything
        seu.sub <- subset(seu.obj, idents = x)
        seu.sub@meta.data[[groupBy]] <- droplevels(seu.sub@meta.data[[groupBy]])
               
        #use FindMarkers buried in the the custom vilnSplitComp function to complete stats to ID DEGs
        geneList <- vilnSplitComp(seu.obj = seu.sub, groupBy = groupBy, refVal = comparison, 
                                  outDir = outDir, outName = outName, 
                                  saveOut = F, saveGeneList = saveGeneList, returnGeneList = T, 
                                  contrast = contrast, pValCutoff = pValCutoff, logfc.threshold = logfc.threshold, ...
                                 ) 
        geneList <- geneList[[1]]
        geneList$gene <- rownames(geneList)
        
        #extract average expression values for plotting (FYI: AverageExpression() is applied to non-logged data!!)
        Idents(seu.sub) <- comparison
        avg.seu.sub <- log1p(AverageExpression(seu.sub, verbose = FALSE)$RNA)
        avg.seu.sub <- as.data.frame(avg.seu.sub)
        
        #convert groups to "X" and "Y" for plotting where contrast[2] is on x-axis and contrast[1] is on y-axis
        avg.seu.sub <- avg.seu.sub[ ,contrast[c(2,1)]]
        colnames(avg.seu.sub) <- c("X","Y")

        #do some manupulation to label select data points
        #this is recomened approach
        if(!useLineThreshold){
            #join the coordinates df with the DEG results then run filtering to select data points to label
            #this labels the top n selected by the user, ranked by abs(X*Y) to find the data points furthest from n.s.

            avg.seu.sub <- avg.seu.sub %>% mutate(gene=rownames(avg.seu.sub)) %>% 
            left_join(geneList, by = "gene") %>% 
            mutate(direction=case_when(avg_log2FC > 0 ~ "up",
                                       avg_log2FC < 0 ~ "dwn",
                                       is.na(avg_log2FC) ~ "NA"),
                   direction=ifelse(gene %in% rownames(geneList)[grepl(paste(labsHide,collapse="|"), row.names(geneList))], "exclude", direction),
                   residual=case_when(direction == "up" ~ abs(Y-X),
                                      direction == "dwn" ~ abs(Y-X),
                                      direction == "NA" ~ 0,
                                      direction == "exclude" ~ 0)) %>% arrange(desc(residual)) %>% group_by(direction) %>% 
            mutate(lab=ifelse(row_number() <= labCutoff
                              & direction != "NA" 
                              & direction != "exclude"
                              & gene %in% c(rownames(geneList),addLabs), gene, ifelse(gene %in% addLabs, gene, NA)),
                   lab_col=case_when(direction == "up" ~ colUp,
                                     direction == "dwn" ~ colDwn,
                                     direction == "NA" ~ "black",
                                     direction == "exclude" & avg_log2FC > 0 ~ colUp,
                                     direction == "exclude" & avg_log2FC < 0 ~ colDwn)
                  )
        }else{
            avg.seu.sub <- avg.seu.sub %>% mutate(gene=rownames(avg.seu.sub),
                                                  up=threshold+1*X, 
                                                  dn=-threshold+1*X,
                                                  direction=case_when(Y > up ~ "up",
                                                                      Y < dn ~ "dwn",
                                                                      Y < up && Y > dn ~ "NA"),
                                                  residual=case_when(direction == "up" ~ Y-up,
                                                                     direction == "dwn" ~ dn-Y,
                                                                     direction == "NA" ~ 0),
                                                 ) %>% group_by(direction) %>% 
            arrange(desc(residual)) %>% 
            mutate(lab=ifelse(row_number() <= labCutoff 
                              & direction != "NA" 
                              & gene %in% c(rownames(geneList),addLabs), gene, ifelse(gene %in% addLabs, gene, NA)),  
                   lab_col=case_when(direction == "up" ~ colUp,
                                     direction == "dwn" ~ colDwn,
                                     direction == "NA" ~ "black")
                  )
        }
        
        #check if there was at least 1 DGE or force return == T before saving the plot
        if(length(na.omit(avg.seu.sub$lab)) > 0 | forceReturn == T){
            outfile <- paste0(outDir, outName, "_", gsub(" ", "_", x),"_linear_deg_by_", gsub(" ", "_", groupBy),".png")
            
            #make the plot
            p <- ggplot(data=avg.seu.sub, aes(x = X, y = Y, label=lab)) + 
            ggtitle(x, 
                    if(subtitle) {subtitle = paste("Average gene expression (", contrast[1]," vs ", contrast[2],")", sep = "")}
                   ) +
            geom_point(color = avg.seu.sub$lab_col) + 
            labs(x = contrast[2], y = contrast[1]) +
            {if(thresLine)geom_abline(intercept = threshold, slope = 1)} +
            {if(thresLine)geom_abline(intercept = -threshold, slope = 1)} + 
            geom_label_repel(max.overlaps = Inf, size=5, color = avg.seu.sub$lab_col, max.time = 2, max.iter = 10000000) + 
            theme_classic() + 
            {if(noTitle)theme(axis.title = element_text(size= 20),
                  axis.text = element_text(size= 14),
                  title = element_blank(),
                  plot.subtitle = element_text(size= 14)
                 )
             } + 
            {if(!noTitle)theme(axis.title = element_text(size= 20),
                  axis.text = element_text(size= 14),
                  title = element_text(size= 24),
                  plot.subtitle = element_text(size= 14)
                 )
                
            }
        
            #save the plot
            ggsave(outfile)
            
            #if request to return, then do so
            if(returnPlots){
                return(p)
            }
            
            #if request to return, then do so -- this returns features upregulated in contrast[1]
            if(returnUpList){
                up <- avg.seu.sub[avg.seu.sub$direction == "up",]
                upList <- up$gene
                return(upList)
            }
            
            #if request to return, then do so -- this returns features upregulated in contrast[2]
            if(returnDwnList){
                dwn <- avg.seu.sub[avg.seu.sub$direction == "dwn",]
                dwnList <- dwn$gene
                return(dwnList)
            }
            
            }
    })
}

############ formatUMAP ############
formatUMAP <- function(
    plot = NULL, 
    smallAxes = F
    ){
    
    plot <- plot + labs(x = "UMAP1", y = "UMAP2") +
        theme(
            axis.text = element_blank(), 
            axis.ticks = element_blank(),
            axis.title = element_text(size= 20),
            plot.title = element_blank(),
            title = element_text(size= 20),
            axis.line = element_blank(),
            panel.border = element_rect(color = "black",
                                        fill = NA,
                                        linewidth = 2)
           )
    
    if(smallAxes){
        axes <- ggplot(data.frame()) + labs(x = "UMAP1", y = "UMAP2") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
        scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + 
        theme_classic() + 
        theme(axis.line = element_line(colour = "black", 
                                       arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                                     ends = "last", type = "closed"),
                                       linewidth = 1.5
                                      ),
              axis.title = element_text(colour = "black", size = 16),
              axis.ticks = element_blank(),
              axis.text = element_blank(), 
              panel.border = element_blank(),
              panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()
             ) + coord_cartesian(expand = FALSE, xlim = c(0, NA), ylim = c(0, NA))
        
        plot <- plot + theme(axis.title = element_blank(),
                             panel.border = element_blank(),
                             plot.margin = unit(c(-7, -7, 14, 14), "pt")
                            ) + NoLegend()

        plot <- plot + inset_element(axes,left= 0,
                                bottom = 0,
                                right = 0.25,
                                top = 0.25,
                                align_to = "full"
                               )
    }
    
    return(plot)
}


############ cusLabels ############
#this function requires a UMAP plot gerneated using DimPlot with label = T, label.box = T
cusLabels <- function(
    plot = NULL, 
    shape = 21, 
    labCol = "black", 
    size = 8, 
    alpha = 1, 
    rm.na = T, 
    nudge_x = NULL, 
    nudge_y = NULL, 
    textSize = 4,
    ...
                  ) {
    
    #extract label coords and colors
    g <- ggplot_build(plot)
    
    labCords <- as.data.frame(g$data[2]) #add error if labels are not present
    labCordz <- labCords[,c("fill","x","y","label","colour")]

    colnames(labCordz) <- c("colour", "UMAP1", "UMAP2", "clusterID", "labCol")
    
    labCordz <- labCordz[order(labCordz$clusterID),]
    labCordz$labCol <- labCol
    
    if(rm.na == T){
        labCordz <- na.omit(labCordz)
    }
    
    if(!is.null(nudge_x)){
        labCordz$UMAP1 <- labCordz$UMAP1 + nudge_x
    }
    
    if(!is.null(nudge_y)){
        labCordz$UMAP2 <- labCordz$UMAP2 + nudge_y
    }

    #remove old labels
    plot$layers[2] <- NULL

    #add labels to the stripped plot to create final image
    plot <- plot + geom_point(data = labCordz, aes(x = UMAP1, y = UMAP2),
                            shape=shape,
                            size=size,
                            fill=labCordz$colour,
                            stroke=1,
                            alpha=alpha,
                            colour="black") +
    geom_text(data = labCordz, size = textSize, mapping = aes(x = UMAP1, y = UMAP2), label = labCordz$clusterID, color = labCordz$labCol)
    
    plot <- formatUMAP(plot, ...)
    return(plot)
}


############ freqPlots ############

freqPlots <- function(seu.obj = NULL, groupBy = "clusterID", refVal = "orig.ident", comp = "cellSource", colz = NULL, namez = NULL, method = 1, nrow = 3, title = F, legTitle = NULL, no_legend = F, showPval = T
                       ) {
    
    if(method == 1){
    fq <- prop.table(table(seu.obj@meta.data[[groupBy]], seu.obj@meta.data[,refVal]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c(groupBy, 
                                                               refVal)) 
        } else if(method == 2){
        
        Idents(seu.obj) <- refVal
        set.seed(12)
        seu.sub <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data[[refVal]])))

        fq <- prop.table(table(seu.sub@meta.data[[refVal]], seu.sub@meta.data[,groupBy]), 2) *100       
        df <- reshape2::melt(fq, value.name = "freq", varnames = c(refVal,
                                                                   groupBy))
    
    } else{print("Method is not recognized. Options are method = 1 (by sample) OR method = 2 (by cluster).")
          break()
          }
    
    df <- df %>% left_join(seu.obj@meta.data[ ,c(refVal,comp)][!duplicated(seu.obj@meta.data[ ,c(refVal,comp)][,1]),], by = refVal)
    #df <- df[grepl("\\b9\\b|12|28|33|41",df$clusterID),] #12 # add grep for only certain groupBy vals
    #colnames(df) <- c("clusterID","orig.ident","freq","cellSource")
    
    if(!is.null(colz)){
        if(length(colz) == 1){
            df <- df %>% left_join(seu.obj@meta.data[ ,c(refVal,colz)][!duplicated(seu.obj@meta.data[ ,c(refVal,colz)][,1]),], by = refVal)
            warn1 <- F
        }else{print("Warning: It is ideal if colors are stored in meta.data slot of the Seurat object. Colors may be mismatched.")
              warn1 <- T
             }
    }else{warn1 <- F}
    
    if(is.null(namez)){namez <- refVal}
    
    if(ifelse(length(namez) > 1, "", namez) != refVal){
         if(length(namez) == 1){
             df <- df %>% left_join(seu.obj@meta.data[ ,c(refVal,namez)][!duplicated(seu.obj@meta.data[ ,c(refVal,namez)][,1]),], by = refVal)
             warn2 <- F
         }else{print("Warning: It is ideal if names are stored in meta.data slot of the Seurat object. Names may be mismatched.")
               warn2 <- T
               }
    }else{warn2 <- F}
    
    p <- ggplot(df, aes_string(y = "freq", x = comp)) + 
    labs(x = NULL, y = "Proportion (%)") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), 
          strip.background = element_rect(fill = NA, color = NA), 
          strip.text = element_text(face = "bold"), 
          axis.ticks.x = element_blank(), 
          axis.text = element_text(color = "black")         )
    
    if(is.null(legTitle)){
        legTitle <- refVal
    }
    
    
    #note these are not corrected p-values
    pi <- p + facet_wrap(groupBy, scales = "free_y", nrow = nrow) + 
    guides(fill = "none") + 
    geom_boxplot(aes_string(x = comp), alpha = 0.25, outlier.color = NA) + 
    geom_point(size = 2, position = position_jitter(width = 0.25),
               aes_string(x = comp, y = "freq", color = refVal)) +
    labs(color = legTitle) +
    {if(showPval){ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3)}} + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme(panel.grid.major = element_line(color = "grey", linewidth = 0.25),
          #legend.position = "none",
          text = element_text(size = 12)
          ) +                    
    {if(!is.null(colz)){
        scale_color_manual(labels = if(warn1){
            namez
        }else{
            levels(df[[namez]])
        },
                           values = if(warn2){
                               colz
                           }else{
                               levels(df[[colz]])
                           })}} + {if(no_legend == T){NoLegend()}}
                                           
    return(pi)

}

############ vilnSplitComp ############

vilnSplitComp <- function(seu.obj = NULL, groupBy = "clusterID", refVal = "cellSource", outName = "", nPlots = 9, saveOut = T, saveGeneList = F, returnGeneList = F, contrast = NULL, filterOutFeats = c("^MT-","^RPL","^RPS"), logfc.threshold = 0.5, pValCutoff = 0.01, ...,
                          outDir = "./output/"
                       ) {
    
    #set the groups to contrast
    DefaultAssay(seu.obj) <- "RNA"
    seu.obj$cluster.condition <- paste(seu.obj@meta.data[[groupBy]], seu.obj@meta.data[[refVal]], sep = "_")
    Idents(seu.obj) <- "cluster.condition"
    
    p <- lapply(levels(seu.obj@meta.data[[groupBy]]), function(x) {

        #split the groups by condition
        comp1 <- paste(x, gsub(" ", "_", contrast[1]), sep = "_") 
        comp2 <- paste(x, gsub(" ", "_", contrast[2]), sep = "_")
        
        #use FindMarkers to complete DGE analysis -- default parameters used
        message(paste0("Contrasting - (1) ", comp1," - VS - (2) ", comp2))
        cluster.markers <- FindMarkers(seu.obj, ident.1 = comp1,  ident.2 = comp2, min.pct = 0, logfc.threshold = logfc.threshold, ...) 
        message(paste0("Number of DEGs prior to p_val_adj filter:"))
        print(cluster.markers %>% mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>% group_by(direction) %>% summarize(nRow = n()))
        
        #filter the results based on the user specified adjusted p-value
        cluster.markers <- cluster.markers %>% filter(p_val_adj < pValCutoff)
        message(paste0("Number of DEGs following p_val_adj filter:"))
        print(cluster.markers %>% mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>% group_by(direction) %>% summarize(nRow = n()))

        #save the output DGE results if requested
        if(saveGeneList){
            cluster.markers$cellType <- x
            write.csv(cluster.markers, file = paste0(outDir, outName, "_", gsub(" ", "_", x) ,"_geneList.csv"))
        }

        #return the gene list if requested
        if(returnGeneList){
            return(cluster.markers)
        }
        
        #filter out requested features using regex
        cluster.markers <- cluster.markers[!grepl(paste(filterOutFeats,collapse="|"), row.names(cluster.markers)), ]
        
        #plot the desired number of DEGs by subsetting the DEGs
        cluster.markers.toPlot <- cluster.markers %>% top_n(n = -nPlots, wt = p_val)
        if(length(rownames(cluster.markers.toPlot)) > 0){
            seurat.vlnplot <- VlnPlot(object = seu.obj,
                                      idents = c(comp1,comp2),
                                      features = rownames(cluster.markers.toPlot),
                                      split.by = refVal
            )
            
            if(saveOut){
                png(file = paste0(outDir, outName,"_", x,"_top",nPlots,"_comp_vlnPlot.png"), width=2520, height=1460)
                print(seurat.vlnplot)
                dev.off()
            }
        }
    })
}
############ loadMeta ############

loadMeta <- function(
    seu.obj = NULL, 
    seu.file = NULL, 
    metaFile = "", 
    groupBy = "clusterID", 
    metaAdd = "majorID",
    save = FALSE, 
    outName = "", 
    header = TRUE
){

    #make this an lapply incase you have a vector of meta data to add
    metaData <- read.csv(file = metaFile, header = header)
    new.cluster.ids <- metaData[[metaAdd]]
    seu.obj <- SetIdent(seu.obj, value = groupBy)

    names(new.cluster.ids) <- metaData[[groupBy]]
    
    if(is.factor(seu.obj@meta.data[ ,groupBy])){
        matchOrder <- levels(seu.obj@meta.data[ ,groupBy])
    } else{
        matchOrder <- unique(seu.obj@meta.data[ ,groupBy])
    }
    
    new.cluster.ids <- new.cluster.ids[match(matchOrder,
                          names(new.cluster.ids))]
    seu.obj <- RenameIdents(seu.obj, new.cluster.ids)
    seu.obj@meta.data[[metaAdd]] <- factor(Idents(seu.obj), levels = unique(unname(new.cluster.ids)))

    return(seu.obj)

    #make it so you can save the updated rds

}

############ vilnSplitCompxGene ############

vilnSplitCompxGene <- function(seu.obj = NULL, groupBy = "clusterID_sub", comp = "cellSource", metaAdd = "majorID", features = "", labelData = NULL,
                               cols = c("mediumseagreen","mediumpurple1"), save = FALSE, outName = "", outDir = "",
                               height = 4, width = 8
                              ){
    
    DefaultAssay(seu.obj) <- "RNA"
    Idents(seu.obj) <- groupBy

    rectangles <- data.frame(
        xmin = seq(1,length(levels(seu.obj)), by = 2)+0.99,
        xmax = seq(1,length(levels(seu.obj)), by = 2)+1.01,
        ymin = 6,
        ymax = 0
    )

    
    p <- VlnPlot(
        features = features,
        object = seu.obj,
        split.by = comp,
        cols = cols,
        stack = TRUE,
        flip = TRUE
    ) + theme(legend.position = "top",
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  plot.margin = margin(0, 0, 0, 3, "pt"),
                  legend.text = element_text(margin = margin(r = 10, unit = "pt"))
                 ) & theme(axis.text.x = element_blank())
        
        if(save == T){
            outfile <- paste0(outDir, outName,"_", x,"_comp_vlnPlot.png") 
            ggsave(plot = p, outfile, height = height, width = width) 
        }
        
    
#         #plot lablels for bottom of plot
#         g <- ggplot(colArray, aes(x = clusterID_sub, y = 1, fill = clusterID_sub, label = clusterID_sub)) + 
#     geom_tile(fill = "transparent") +
#     geom_point(shape = 21, stroke = 0.75, size = 7) +
#     geom_text(fontface = "bold", size = 3, color = colArray$labCol) + theme_bw(base_size = 12) +
#     scale_fill_manual(values = colArray$colour) + scale_y_discrete(expand = c(0, 0)) +
#     theme(legend.position = "none", panel.spacing = unit(0, "lines"),
#               panel.background = element_blank(), 
#               panel.border = element_blank(),
#               plot.background = element_blank(), 
#               plot.margin = margin(0, 0, 0, 3, "pt"),
#               axis.title.y = element_blank(),
#               axis.ticks.y = element_blank(),
#               axis.ticks.x = element_blank(),
#               axis.text = element_blank(),
#                   panel.grid.major = element_blank(), 
#                   panel.grid.minor = element_blank()) + xlab("Cluster") 
#      pi <- plot_grid(p, g, ncol = 1, rel_heights = c(0.85, 0.15), align = "v", axis = "lr") 
    
    return(p)
}

############ getPB ############

#where bioRep == annotation vector cell & orig.ident data (with levels for each patient)
getPb <- function(mat.sparse, bioRep) {

    #extract count sums for each cell within a biological replicate
    collapsedCnts <- lapply(levels(bioRep$cellSource), function(repicate){
        
        cells <- row.names(bioRep)[bioRep$cellSource == repicate]
        pseudobulk <- Matrix::rowSums(mat.sparse[ ,cells])

        return(pseudobulk)
        
        })

    #bind the values for each sample and add the column names
    mat.summary <- do.call(cbind, collapsedCnts)
    colnames(mat.summary) <- levels(bioRep$cellSource)

   return(mat.summary)
}

############ createPB ############

createPB <- function(seu.obj = NULL, groupBy = "clusterID_sub", comp = "cellSource", biologicalRep = "orig.ident",
                     clusters = NULL, outDir = "", min.cell = 5, lowFilter = TRUE, dwnSam = FALSE, featsTOexclude = NULL, cnts = TRUE
                    #  , grepTerm = NULL, grepLabel = NULL #improve - fix this so it is more functional -- DONE 11.27.2023
                    ){

    #set active.ident to the level to loop by
    Idents(seu.obj) <- groupBy
    
    #extract number of biological replicates in the grouping system
    ztest <- dim(table(seu.obj@meta.data[[biologicalRep]]))
    
    #if specific clusters no provied, loop through all levels in active.ident
    groupz <- ifelse(is.null(clusters), levels(seu.obj), list(clusters))

    #loop though all levels of active.ident
    test <- lapply(levels(seu.obj), function(x){

        #subset the data by group
        seu.sub <- subset(seu.obj, idents = x)

        #determine value to downsample by and extract metadata - uses bioRep with lowest number of cells > "5"
        z <- table(seu.sub@meta.data[[biologicalRep]])
        z <- z[z > min.cell]
        zKeep <- names(z)

        if(length(zKeep) >= 0.5*ztest){
        
            #remove samples that have insufficent cell numbers
            Idents(seu.sub) <- biologicalRep
            seu.sub.clean <- subset(seu.sub, idents = zKeep)
    
            seu.sub.clean@meta.data[[biologicalRep]] <- as.factor(seu.sub.clean@meta.data[[biologicalRep]])
            seu.sub.clean@meta.data[[biologicalRep]] <- droplevels(seu.sub.clean@meta.data[[biologicalRep]])
            
            z <- table(seu.sub.clean@meta.data[[biologicalRep]], seu.sub.clean@meta.data[[comp]]) %>% melt() %>% filter(value != 0)
            
            #extract and save metadata in a data.frame
            z <- as.data.frame(z)
            colnames(z) <- c("sampleID", "groupID", "nCell")
            z$clusterID <- toString(x) #add clusterID value
            # z$groupID <- ifelse(grepl(grepTerm, z$sampleID), grepLabel[1], grepLabel[2]) #fix this - hardcoded and will only work for PBMCs improve this

            #down sampling is not necissary, some people reccomend it, but I now reccomend skipping this (default)
            if(dwnSam){

                #identitfy sample with fewest number of cells
                ds <- min(z$nCell[z$nCell > min.cell])
                msg <- paste0("Downsampling cluster: ", x," at a level of ", ds, " cells per replicate.")
                message(msg)
        
                #randomly downsample the subset data
                Idents(seu.sub.clean) <- biologicalRep
                set.seed(12)
                seu.sub.clean <- subset(x = seu.sub.clean, downsample = ds)
            }
            
            #extract required data for pseudobulk conversion
            if(cnts){
                mat <- seu.obj@assays$RNA@layers$counts
            }else{
                mat <- seu.obj@assays$RNA@layers$data
            }

            #bring over gene sybmols and cell names
            colnames(mat) <- colnames(seu.obj)
            rownames(mat) <- rownames(seu.obj)
            
            #remove features that have less than 10 cells with a non-zero expression value
            if(lowFilter){
                mat <- mat[rowSums(mat > 1) >= 10, ]
            }
                        
            #extract data needed to make pseudobulk matrix
            bioRep <- as.data.frame(seu.sub.clean@meta.data[[biologicalRep]])
            colnames(bioRep) <- "cellSource"
            row.names(bioRep) <- colnames(seu.sub.clean)
            bioRep$cellSource <- as.factor(bioRep$cellSource)
            #use custom function to convert to pseudobulk
            pbj <- getPb(mat, bioRep)
            #optionally exclude features from conversion (HBM, PPBP, MT-, RPS- are feats to consider excluding)
            if(!is.null(featsTOexclude)){
                pbj <- pbj[!rownames(pbj) %in% featsTOexclude, ]
            }
            #log the number of reps included
            if(length(colnames(pbj))-1 != ztest){
                msg <- paste0("INFO: During pseudobluk conversion of ", x, " the following samples were included: ", paste(as.list(colnames(pbj)),collapse=" "),"\nINFO: min.cell was set to: ", min.cell,"\n")
                message(msg)
            } else {
                msg <- "All replicates were used for psudobluk conversion"
                message(msg)
            }
            #save the matrix
            write.csv(pbj, file = paste0(outDir, x, "_pb_matrix.csv"))
        } else {
            msg <- paste0("Unable to downsample cluster: ", x, " due to insufficent cell numbers")
            message(msg)
        } 
        return(z)
    })
 
    #collect the metadata
    df <- do.call(rbind, test)
    df <- na.omit(df)
    
    #save the metadata 
    csvOut <- paste0(outDir,groupBy ,"_deg_metaData.csv")
    write.csv(df, file = csvOut)
}

############ pseudoDEG ############
### BUG: appears if special characters are used in sample names
### fails at gathered_top20_sig <- meta %>% inner_join(gathered_top20_sig, by = c("sampleID" = "samplename"))

# contrast will be idents.1_NAME vs idents.2_NAME !!!
pseudoDEG <- function(
    metaPWD = "", 
    padj_cutoff = 0.1,
    lfcCut = 0.58,
    outDir = "",
    outName = "", 
    idents.1_NAME = NULL,
    idents.2_NAME = NULL,
    returnDDS = F,
    inDir = "",
    title = "", 
    fromFile = T,
    meta = NULL,
    pbj = NULL,
    returnVolc = F,
    paired = F,
    pairBy = "", 
    minimalOuts = F,
    saveSigRes = T,
    topn=c(20,20),
    filterTerm = NULL,
    addLabs = NULL,
    mkDir = F,
    test.use = "Wald",
    dwnCol = "blue",
    stblCol = "grey",
    upCol = "red",
    labSize = 3,
    strict_lfc = F,
    featuresToExclude = NULL
    ){

    if(fromFile){
        files <- list.files(path = inDir, pattern="pb_matrix.csv", all.files=FALSE,full.names=FALSE)
        clusters <- unname(sapply(files, function(x) {unlist(strsplit(x, split = "_pb_"))[1]}))
        outfileBase <- paste0(outDir, outName, "_cluster_")
    }else{
        clusters <- outName
        outfileBase <- outDir
    }

    lapply(clusters, function(x) {
        if(fromFile){
            inFile <- paste0(inDir, x,"_pb_matrix.csv")
            pbj <- read.csv(file = inFile, row.names = 1)
            pbj <- pbj[!apply(pbj==0, 1, all),]
            
            meta <- read.csv(file = metaPWD, row.names = 1)
            meta <- meta[meta$clusterID == x,]
            meta[,colnames(meta)] <- lapply(meta[,colnames(meta)] , factor)
        }
        if(mkDir){
            outDir <- paste0(outDir, "/", x, "/")
            dir.create(outDir)
            outfileBase <- paste0(outDir, outName, "_cluster_")
        }
        
        if(!is.null(featuresToExclude)){
            pbj <- pbj[!rownames(pbj) %in% featuresToExclude , ]
        }

        if(paired){
            dds <- DESeqDataSetFromMatrix(round(pbj), 
                                          colData = meta,
                                          design = formula(paste0("~ groupID + ",noquote(pairBy))))
            
            message("Using the following contrast formula:\n",
                    formula(paste0("~ groupID + ",noquote(pairBy))))
        }else{
            dds <- DESeqDataSetFromMatrix(round(pbj), 
                                          colData = meta,
                                          design = ~ groupID)
            
            message("Using the following contrast formula:\n",
                    formula(~ groupID))
        }
        if(returnDDS){
            return(dds)
        }

        #transform and plot the data with PCA
        rld <- varianceStabilizingTransformation(dds, blind=TRUE)

        if(!minimalOuts){
            outfile <- paste0(outfileBase, x,"_pca.png")
            p <- DESeq2::plotPCA(rld, intgroup = "groupID")
            ggsave(outfile, width = 7, height = 7)
            
            outfile <- paste0(outfileBase, x,"_pca2.png")
            p <- DESeq2::plotPCA(rld, intgroup = "sampleID")
            ggsave(outfile, width = 7, height = 7)
        }
        
        #set up QC to evaluate treatment seperation by heatmap & plot
        rld <- rlog(dds, blind=FALSE)
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)

        if(!minimalOuts){
            outfile <- paste0(outfileBase, x,"_pheatmap.png")
            p <- pheatmap::pheatmap(rld_cor)
            ggsave(p, file = outfile)
        }
        if(test.use == "LRT"){
            if(paired){
                dds <- DESeq(dds,
                             test = "LRT",
                             reduced = formula(paste0("~",noquote(pairBy)))
                            )
                res <- results(dds)
            } else{
                dds <- DESeq(dds,
                             test = "LRT",
                             reduced = ~1)
                res <- results(dds)
            }
        } else if (test.use == "Wald"){
            dds <- DESeq(dds)
            
            #extract Wald test results
            if(!is.null(idents.1_NAME) | !is.null(idents.2_NAME)){
                contrast <- c("groupID", idents.1_NAME, idents.2_NAME)
            } else{
                contrast <- c("groupID", unique(meta$groupID)[1], unique(meta$groupID)[2])
                print("Variables idents.1_NAME and/or idents.2_NAME not specify, this may impact directionality of contrast - confirm results are as expect and/or specify ident names.")
            }

            if(strict_lfc){
                lfcCut_strict <- lfcCut
            } else{
                lfcCut_strict <- 0
            }
            
            #perform logFoldChange and shrinkage
            message(paste("The lfc calculation will be based on", contrast[2], 
                          "versus", contrast[3]))
            res <- results(dds, 
                           contrast = contrast,
                           alpha = padj_cutoff,
                           lfcThreshold = lfcCut_strict, #this increases stringency
                           cooksCutoff = FALSE)
            summary(res)

            res <- lfcShrink(dds, 
                             contrast = contrast,
                             res = res,
                             type = "normal") #would prefer to use something other than normal

        } else{
            message("\nERROR: test.use variable is not valid. Please specify Wald or LRT then try again.\n")
            break
        }

        #extract the results and save as a .csv
        res_tbl <- res %>% data.frame() %>%
        rownames_to_column(var="gene") %>% as_tibble()
        sig_res <- res_tbl %>% filter(padj < padj_cutoff, abs(log2FoldChange) > lfcCut) %>% arrange(padj)
        if(saveSigRes){
            sig_res$gs_base <- toupper(x)
            write.csv(sig_res,
                      file = paste0(outfileBase, x,"_all_genes.csv"),
                      quote = FALSE,
                      row.names = FALSE)
        }
        
        if(test.use == "LRT"){
            return(res)
        } else{
            #get nomarlized counts and plot top 20 DEGs
            normalized_counts <- counts(dds, 
                                        normalized = TRUE)
            top20_sig_genes <- sig_res %>% arrange(padj) %>% pull(gene) %>% head(n=20)
            top20_sig_norm <- data.frame(normalized_counts) %>%
            rownames_to_column(var = "gene") %>% filter(gene %in% top20_sig_genes)
            gathered_top20_sig <- top20_sig_norm %>%
            gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
            gathered_top20_sig <- meta %>% inner_join(gathered_top20_sig, by = c("sampleID" = "samplename")) #need more metadata

            if(dim(gathered_top20_sig)[1] > 0){ 
                if(!minimalOuts){
                    outfile <- paste0(outfileBase, x,"_genePlot.png")
                    p <- ggplot(gathered_top20_sig) +
                    geom_point(aes(x = gene, 
                                   y = normalized_counts, 
                                   color = groupID), #need to fix samplename & change to groupID
                               position=position_jitter(w=0.1,h=0)) +
                    scale_y_log10() +
                    xlab("Genes") +
                    ylab("log10 Normalized Counts") +
                    ggtitle("Top 20 Significant DE Genes") +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    theme(plot.title = element_text(hjust = 0.5))
                    ggsave(p, file = outfile)
                }

                #extract sig DEGs based on normalized data and plot heat map
                sig_norm <- data.frame(normalized_counts) %>%
                rownames_to_column(var = "gene") %>% filter(gene %in% sig_res$gene)
                rownames(sig_norm) <- sig_norm$gene
                colAnn <- as.data.frame(meta[,"groupID"], colname = "groupID")
                rownames(colAnn) <- meta[,"sampleID"]
                colnames(colAnn) <- "groupID"
                heat_colors <- brewer.pal(6, "YlOrRd")

                #set threshold and flag data points then plot with ggplot
                if(!is.null(filterTerm)){
                    feats_to_hide <- res_tbl$gene[grepl(paste(filterTerm, collapse = "|"), res_tbl$gene)]
                } else{
                    feats_to_hide <- ""
                }
                
                # Clean results and add singnifiance tags
                res_table_thres <- res_tbl %>% 
                    mutate(
                        threshold = ifelse(padj < padj_cutoff & abs(log2FoldChange) >= lfcCut, 
                                           ifelse(log2FoldChange > lfcCut, 'Up', 'Down'), 'Stable')
                    ) %>%
                    na.omit() %>%
                    arrange(padj) %>%
                    filter(! gene %in% feats_to_hide)
                
                # Tag features to label
                cntUp <- nrow(res_table_thres[res_table_thres$threshold == "Up",])
                cntDwn <- nrow(res_table_thres[res_table_thres$threshold == "Down",])
                top20_up <- res_table_thres[res_table_thres$threshold == "Up",] %>% do(head(., n = topn[1]))
                top20_down <- res_table_thres[res_table_thres$threshold == "Down",] %>%  do(head(., n = topn[2]))
                res_table_thres <- res_table_thres %>% 
                    mutate(
                        label = ifelse(gene %in% c(top20_up$gene, top20_down$gene, addLabs), gene, NA)
                    )

                if(fromFile){
                    if(is.null(title)){
                        title <- paste0(idents.1_NAME, " vs ",idents.2_NAME, "within", x)
                    }
                }

                p <- ggplot(data = res_table_thres, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
                    geom_vline(xintercept = c(-lfcCut,lfcCut), lty = 4, col="black", lwd = 0.8) +
                    geom_hline(yintercept = -log10(padj_cutoff), lty = 4, col="black", lwd = 0.8) +
                    geom_point(alpha = 0.4, size = 3.5) +
                    geom_label_repel(max.overlaps = Inf, size = labSize, label = res_table_thres$label, 
                                     show.legend = FALSE) +
                    labs(
                        x = "log2(fold change)",
                        y = "-log10(padj)",
                        title = title
                    ) + 
                    scale_color_manual(
                        values = c("Down" = dwnCol, "Stable" = stblCol, "Up" = upCol), 
                        labels = c("Down" = paste0("Down (", cntDwn, ")"), 
                                   "Stable" = "Stable", 
                                   "Up" = paste0("Up (", cntUp, ")"))
                    ) +
                    theme_bw() +
                    theme(
                        plot.title = element_text(size = 20, hjust = 0.5), 
                        legend.position = "right", 
                        legend.title = element_blank()
                    )
                
                ggsave(p, file = paste0(outfileBase, x, "_volcano.png"))
                
                if(returnVolc){
                    return(p)
                }
            }
        }
    })
}
############ btwnClusDEG ############
#work in progress             - need to to fix doLinDEG option ### NOTE: cannot have special char in ident name
btwnClusDEG <- function(seu.obj = NULL,groupBy = "majorID_sub", idents.1 = NULL, idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, topn=c(20,20), strict_lfc = F,
                        minCells = 25, outDir = "", title = NULL, idents.1_NAME = "", idents.2_NAME = "", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = F, dwnSam = T, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3, filterTerm = NULL
                    ){
    
    seu.integrated.obj <- seu.obj
    DefaultAssay(seu.integrated.obj) <- "RNA"

    Idents(seu.integrated.obj) <- groupBy
    if(is.null(idents.2)){
        idents.2 <- levels(seu.integrated.obj)[!levels(seu.integrated.obj) %in% idents.1]
    }
    
    seu.sub <- subset(seu.integrated.obj, idents = as.list(append(idents.1,idents.2)))
    
    if(all(idents.1 %in% levels(seu.integrated.obj))){
        if(length(idents.1) > 1){
            grepTerm <- paste(idents.1,collapse="|")
        }else{
            grepTerm <- idents.1
        }
    }else{
        print("Idents not found in seurat object. Please check that they are correct and try again.")
        break()
    }

    if(is.null(title)){
        title <- paste0(gsub(" ", "_", idents.1_NAME), "_vs_",gsub(" ", "_",idents.2_NAME))
    }
    seu.sub$compare <- ifelse(grepl(grepTerm, seu.sub@meta.data[[groupBy]]), "idents.1", "idents.2")
    
    Idents(seu.sub) <- "compare"

    #extract number of biological replicates in the grouping system
    ztest <- dim(table(seu.sub@meta.data[[bioRep]]))

    seu.sub@meta.data$conditional <- paste(seu.sub@meta.data[[bioRep]], seu.sub@meta.data$compare, sep = "_")

    z <- table(seu.sub@meta.data$conditional)
    z <- z[z>minCells]
    zKeep <- names(z)

    Idents(seu.sub) <- "conditional"
    seu.sub.clean <- subset(seu.sub, idents = zKeep)
    #seu.sub.clean <- NormalizeData(seu.sub.clean)
    
    z <- table(seu.sub.clean@meta.data$conditional)
    
    #extract and save metadata in a data.frame
    z <- as.data.frame(z)
    colnames(z) <- c("sampleID", "nCell")
    z$groupID <- ifelse(grepl("idents.1", z$sampleID), idents.1_NAME, idents.2_NAME)
    
    ds <- min(z$nCell)

    if(dwnSam){
    message <- paste0("Downsampling at a level of ", ds, " cells per replicate")
    print(message)
        
    #randomly downsample the subset data
        Idents(seu.sub.clean) <- "conditional"
        set.seed(setSeed)
        seu.sub.clean <- subset(x = seu.sub.clean, downsample = ds)
    }

    #extract required data for pseudobulk conversion
    mat <- seu.sub.clean@assays$RNA$counts
    if(lowFilter){
        mat <- mat[rowSums(mat > 1) >= 10, ]
    }
    
    bioRep <- as.data.frame(seu.sub.clean@meta.data$conditional)
    colnames(bioRep) <- "cellSource"
    row.names(bioRep) <- colnames(seu.sub.clean)
    bioRep$cellSource <- as.factor(bioRep$cellSource)

    #use custom function to convert to pseudobulk
    pbj <- getPb(mat, bioRep)
    pbj <- pbj[rowSums(pbj[])>0,] # remove any rows that are all zeros

    ### pseudobulk DEG using DEseq2 backbone ###
    meta <- z
    
    meta$bioRepPair <- str_split_fixed(meta$sampleID, "_ident.", 2)[,1]
    if(doLinDEG == T){
        seu.sub.clean$compareLinDEG <- ifelse(grepl(grepTerm, seu.sub.clean@meta.data$clusterID), "idents.1", "idents.2")
        seu.sub.clean$groupz <- "cellz"
        linDEG(seu.obj = seu.sub.clean, threshold = 1, thresLine = T, groupBy = "groupz", comparison = "compareLinDEG", outDir = outDir, outName = paste0(gsub(" ", "_", idents.1_NAME), "_vs_",gsub(" ", "_",idents.2_NAME)), colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = T
              )
    }
    print(meta)
    p <- pseudoDEG(padj_cutoff = padj_cutoff, lfcCut = lfcCut, outName = paste0(gsub(" ", "_", idents.1_NAME), "_vs_",gsub(" ", "_",idents.2_NAME)), 
              outDir = outDir, title = title, fromFile = F, meta = meta, pbj = pbj, returnVolc = returnVolc, paired = paired, pairBy = "bioRepPair", strict_lfc = strict_lfc,
                   idents.1_NAME = idents.1_NAME, idents.2_NAME = idents.2_NAME, minimalOuts = T, saveSigRes = T, addLabs = addLabs, topn = topn, dwnCol = dwnCol, stblCol = stblCol,upCol = upCol, labSize = labSize, filterTerm = filterTerm
                     )    
    return(p)
}
############ volcFromFM ############
volcFromFM <- function(seu.obj = NULL, padj_cutoff = 0.01, lfcCut = 0.58, title = "",
                     clusters = NULL, outDir = "", preSub = F
                    ){    
    
    if(preSub == F & !is.null(clusters)){
        seu.obj <- subset(seu.obj,
                          subset = 
                          clusterID_sub == clusters) #this aint right fix
        }

    clusMarkers <- FindMarkers(seu.obj, ident.1 = "OS", group.by = "cellSource")

    clusMarkers <- tibble::rownames_to_column(clusMarkers, "gene")

    #set threshold and flag data points then plot with ggplot
    res_table_thres <- clusMarkers %>% mutate(threshold = ifelse(p_val_adj < padj_cutoff & abs(avg_log2FC) >= lfcCut, 
                                                                 ifelse(avg_log2FC > lfcCut ,'Up','Down'),'Stable')
                                             )

    res_table_thres <- res_table_thres[!is.na(res_table_thres$p_val_adj),]

    res_table_thres.sortedByPval = res_table_thres[order(res_table_thres$p_val_adj),]
    cntUp <- nrow(res_table_thres[res_table_thres$threshold == "Up",])
    cntDwn <- nrow(res_table_thres[res_table_thres$threshold == "Down",])
    res_table_thres.sortedByPval <- res_table_thres.sortedByPval[!grepl("^ENSCAF|^MT-|^RPS|^RPL", res_table_thres.sortedByPval$gene),] 
    top20_up <- res_table_thres.sortedByPval[res_table_thres.sortedByPval$threshold == "Up",] %>% do(head(., n=20))
    top20_down <- res_table_thres.sortedByPval[res_table_thres.sortedByPval$threshold == "Down",] %>%  do(head(., n=20))


    res_table_thres <- res_table_thres %>% mutate(label = ifelse(gene %in% top20_up$gene | gene %in% top20_down$gene, gene, NA),
                                                  shape = ifelse(p_val_adj == 0, 17, 19)
                                                 )
    
    
    outfile <- paste0(outDir, title,"_volcano.png")

    p <- ggplot(data = res_table_thres, 
                aes(x = avg_log2FC, 
                y = -log10(p_val_adj), 
                colour = threshold)) +
    geom_vline(xintercept = c(-0.58,0.58), lty = 4, col="black", lwd = 0.8) +
    geom_hline(yintercept = -log10(padj_cutoff), lty = 4, col="black", lwd = 0.8) +
    geom_point(alpha = 0.4, size = 3.5,
                shape = res_table_thres$shape) +
    geom_label_repel(max.overlaps = Inf, size=3, label = res_table_thres$label, show.legend = FALSE) +
    #xlim(c(-2.5, 2.5)) +
    #ylim(c(0, 3)) +
    labs(x="log2(fold change)",
         y="-log10(padj)",
         title=title) + #paste0("Differential expression (",idents.1_NAME ," vs ",idents.2_NAME ,")")
    #scale_colour_discrete(labels=c(paste0("Down (", cntDwn,")"), "Stable", paste0("Up (", cntUp,")"))) + 
    scale_color_manual(values=c("Down" = "blue", "Stable" = "grey","Up" = "red"), labels=c(paste0("Down (", cntDwn,")"), "Stable", paste0("Up (", cntUp,")"))) +
    theme_bw() +
    theme(plot.title = element_text(size = 20), 
          legend.position = "right", 
          legend.title = element_blank()
         )
    
    ggsave(plot = p, file = outfile)
    
}

############ vilnPlots ############
vilnPlots <- function(
    seu.obj = NULL, 
    inFile = NULL,
    groupBy = "clusterID",
    numOfFeats = 24,
    outName = "",
    outDir = "",
    outputGeneList = T,
    filterOutFeats = c("^MT-", "^RPL", "^RPS"),
    Assay = "RNA",
    returnViln = T,
    features = NULL,
    min.pct = 0.25,
    only.pos = T,
    resume = F,
    resumeFile = NULL
    ){ 
    
    if(!resume){
        if(!is.null(inFile)){
            try(seu.obj <- readRDS(file), silent = T)
        }else if(!is.null(seu.obj)){
            seu.obj <- seu.obj
        }

        DefaultAssay(seu.obj) <- Assay

        Idents(seu.obj) <- groupBy
        ident.level <- levels(seu.obj@active.ident)

        cluster.markers <- FindAllMarkers(seu.obj, only.pos = only.pos, min.pct = min.pct, features = features)
        cluster.markers <- cluster.markers[cluster.markers$p_val_adj < 0.01, ]
        
        if(length(filterOutFeats) > 0){
            cluster.markers <- cluster.markers[!grepl(paste(filterOutFeats,collapse="|"), rownames(cluster.markers)), ]
        }     
        
        if (outputGeneList == TRUE){
            dir.create(outDir)
            outfile <- paste0(outDir, outName, "_",  groupBy, "_gene_list.csv")
            write.csv(cluster.markers, file = outfile)
        }  
        
    }else{
        cluster.markers <- read.csv(file = resumeFile, header = T)
        Idents(seu.obj) <- groupBy
    }
    
    if(returnViln){
        
        lapply(unique(cluster.markers$cluster), function(x) {

            geneList <- cluster.markers[cluster.markers$cluster == x, ]
            geneList <- geneList$gene

            plotList <- head(geneList, n = numOfFeats)

            seurat.vlnplot <- VlnPlot(
                object = seu.obj,
                features = rev(plotList)
            )

            outfile <- paste0(outDir, outName, "_", x, "_top", numOfFeats, "_vlnPlot.png") 
            png(file = outfile, width=2520, height=1460)

            print(seurat.vlnplot)
            dev.off()
        })
    }
    
}
############ singleR ############

singleR <- function(seu.obj = NULL, reduction = NULL, clusters = "clusterID",
                    outName = "",  outDir = ""
                     ){
    
    cntData <- GetAssayData(seu.obj, slot = "data", assay = "RNA")
    Kotliarov <- scRNAseq::KotliarovPBMCData()
    refs <- list(human_ref1 = celldex::HumanPrimaryCellAtlasData(), 
                 human_ref2 = BlueprintEncodeData(),
                 human_ref3 = DatabaseImmuneCellExpressionData(),
                 human_ref4 = NovershternHematopoieticData(), 
                 human_ref5 = MonacoImmuneData()
                )
   
    
    #Perform annotation using each reference
    sapply(names(refs), FUN = function(ref_idx) {
        ref <- refs[[ref_idx]]
        ref_name <- paste("SingleR", ref_idx, sep = ".")
        
        #Annotation is performed on cluster-level rather than default single-cell level
        rst <- SingleR(test = cntData, ref = ref, clusters = seu.obj@meta.data[[clusters]], labels = ref$label.fine)
        
        # Assign predicted cell labels to seurat object
        # SingleR assigns labels for all clusters, even for those not in the reference
        # So use pruned.labels to remove uncertain labels 
        seu.obj@meta.data[[ref_name]] <- rst$pruned.labels[match(seu.obj@meta.data[[clusters]], rownames(rst))]

        #Visualize cell labels
        p <- DimPlot(seu.obj, reduction = reduction, group.by = ref_name, label = TRUE)
        ggsave(plot = p, filename = paste0(outDir,outName,"_",ref_name, ".png"),  width = 10, height = 7)
    })
}

############ cusLeg ############

cusLeg <- function(legend = NULL, colz = 2, rowz = NULL, clusLabel = "clusterID", legLabel = NULL, groupLabel = NULL, colorz = NULL, dotSize = 8,
                   groupBy = "", sortBy = "", labCol = "", headerSize = 6, #add if len labCol == 1 then add col to the df
                   cellHeader = T, bump = 2, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 0, compress_x = 0.75, titleOrder = NULL, spaceBtwnCols = NULL, breakGroups = F, horizontalWrap = F, returnData = F, overrideData = NULL
                     ){
    
    if(is.null(colorz)){
        legend$colorz <- gg_color_hue(length(legend[[clusLabel]]))
        colorz <- "colorz"
        
    }
    
    if(cellHeader == T){
        if(!is.null(titleOrder)) {
            legend[[groupLabel]] <- factor(legend[[groupLabel]], levels = titleOrder)
        }else{
            legend[[groupLabel]] <- factor(legend[[groupLabel]])
        }
        
        header <- legend  %>% group_by(!!as.name(groupLabel)) %>% summarize(n = n()) %>% arrange(factor(!!as.name(groupLabel), 
                           levels = groupBy)) %>% mutate(pos = n+1,
                                                         cum_sum = cumsum(pos),
                                                         dfNum = cum_sum-n
                                                        )  %>% as.data.frame()
        nTerms <- max(header$cum_sum) #add this var outsdide of if
        
        if(is.null(rowz)){
            rowz <- ceiling((length(legend[[clusLabel]])+dim(header)[1])/colz)
        }
    }
    
    if(length(spaceBtwnCols) < colz-1){
        message("Invalid spaceBtwnCols entry. Ensure list is correct length!")
    }
    
    #if there is one header...
    if(rowz < max(header$n) & nrow(header) == 1){
        if(!horizontalWrap){
            leg_y <- rep(seq(rowz, 1, by = -1),colz)[1:header$n]
        }else{
            leg_y <- rep(seq(rowz, 1, by = -1), each = colz)[1:header$n]
        }
        
        leg_y <- leg_y + compress_y

        if(is.null(spaceBtwnCols)){
            if(!horizontalWrap){
                leg_x <- rep(seq(1,colz)/2, each = rowz)[1:header$n]
            }else{
                leg_x <- rep(seq(1,colz)/2, rowz)[1:header$n]
            }
        }else{
            leg_x <- tryCatch(
                {
                    if(!horizontalWrap){
                        rep(cumsum(c(0.5,spaceBtwnCols)), each = rowz)[1:header$n]
                    }else{
                        rep(cumsum(c(0.5,spaceBtwnCols)), rowz)[1:header$n]
                    }
                }, 
                error = function(e) {
                    if(!horizontalWrap){
                        return(rep(seq(1,colz)/2, each = rowz)[1:header$n])
                    }else{
                        return(rep(seq(1,colz)/2, rowz)[1:header$n])
                    }
                    message("Invalid spaceBtwnCols entry. Using defeult spacing of 0.5")
                }
            )
        }
        
        header$header_x <- leg_x[header$dfNum]
        header$header_y <- leg_y[header$dfNum]+1
        headerBump = T
        
    #else there is more than one header
    }else{
        headerGroups <- seq(1,nrow(header), by = ceiling(nrow(header)/colz))
        if(length(headerGroups) != colz){
            headerGroups <- sort(c(headerGroups,seq(1,nrow(header), by = 1)[!seq(1,nrow(header), by = 1) %in% headerGroups][1]))
        }
        headData <- header[headerGroups,] %>% mutate(diff = cum_sum - lag(cum_sum, default = 0))
        
        leg_y <- do.call(c,(lapply(headData$diff,function(x){seq(max(headData$diff),max(headData$diff)-x+1, by = -1)})))
        
        if(is.null(spaceBtwnCols)){
            leg_x <- rep(seq(1,colz)/2, headData$diff)
        }else{
            leg_x <- tryCatch(
                {
                    rep(cumsum(c(0.5,spaceBtwnCols)), headData$diff)
                }, 
                error = function(e) {
                    return(rep(seq(1,colz)/2, headData$diff))
                    message("Invalid spaceBtwnCols entry. Using default spacing of 0.5")
                }
            )
        }
           
        leg_x <- leg_x[1:(dim(header)[1]+dim(legend)[1])]
        leg_y <- leg_y + compress_y
    
        header$header_x <- leg_x[header$dfNum]
        header$header_y <- leg_y[header$dfNum]
    
        leg_x <- leg_x[-c(header$dfNum)]
        leg_y <- leg_y[-c(header$dfNum)]
        headerBump = F
    }
        
    #legend <- legend %>% arrange(factor(!!as.name(groupLabel), levels = groupBy), !!as.name(sortBy))
    legend <- legend %>% arrange(!!as.name(groupLabel), !!as.name(sortBy))
    
    legend$leg_x <- leg_x
    legend$leg_y <- leg_y
    
    if(!is.null(overrideData)){
        legend <- as.data.frame(overrideData[1])
        header <- as.data.frame(overrideData[2])
    }

    p <- ggplot() + geom_point(data = legend, aes(x = leg_x, y = leg_y),
                    shape=21,
                    size=dotSize,
                    fill=legend[[colorz]],
                    stroke=1,
                    colour="black") +
          geom_text(data = legend, size = 4, mapping = aes(x = leg_x, y = leg_y), label = legend[[clusLabel]], colour = legend[[labCol]]) +
          geom_text(data = legend, size = 4, mapping = aes(x = leg_x+0.05, y = leg_y), label = legend[[legLabel]], hjust = 0) + #NoLegend() + 
          geom_text(data = header, size = headerSize, mapping = aes(x = header_x-0.03, y = header_y), label = header[[groupLabel]], hjust = 0, fontface =2) +
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_blank(),
            title = element_blank(),
            axis.line = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) + coord_cartesian(ylim = c(ymin, max(ifelse(headerBump == T, leg_y+1, leg_y))*topBuffer), expand = TRUE, clip = "off") + scale_x_continuous(limits = c(0.47,colz*compress_x))

    if(!returnData){
        return(p)
    }else{
        return(list(legend, header))
    }
        
}

############ majorDot ############

majorDot <- function(seu.obj = NULL, groupBy = "",
                     yAxis = NULL, scale = T,
                     features = "", split.by = NULL, cols = c("lightgrey", "blue"), cluster.idents = F
                    ){
    
    t <- try(head(seu.obj@assays$RNA@scale.data),silent = T)

    if("try-error" %in% class(t)){
        seu.obj <- ScaleData(seu.obj)
    }
                       
    p <- DotPlot(seu.obj,
                 assay = "RNA",
                 features = features,
                 group.by = groupBy,
                 cols = cols,
                 scale = scale,
                 split.by = split.by#,
                 #idents = levels(seu.obj@meta.data[[groupBy]]),
                 #cluster.idents = cluster.idents
                ) +
      geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
      labs(size='Percent\nexpression')  +
      theme(axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.justification='center',
            panel.background = element_rect(fill = "white",colour = NA),
            plot.background = element_rect(fill = "white",colour = NA),
            legend.key.size = unit(1, "line"),
            panel.border = element_rect(color = "black",
                                        fill = NA,
                                        linewidth = 1),
            ) +
      scale_colour_viridis(option="magma", name='Average\nexpression') +
      guides(size=guide_legend(override.aes = list(shape=21, colour="black", fill="white"),
                               label.position = "bottom")) +
      scale_size(range = c(0.5, 8), limits = c(0, 100)) +
      #annotate("rect", xmin = features_cnt$startVal, xmax = features_cnt$endVal, ymin = features_cnt$cluster-0.5, ymax = features_cnt$cluster+0.5, fill = NA, colour = "mediumpurple1", size = 1) +
      {if(!is.null(yAxis)){scale_y_discrete(limits=rev(yAxis))}} +
      guides(color = guide_colorbar(title = 'Scaled\nExpression')) 
      
    
    return(p)
}

############ autoDot ############

autoDot <- function(seu.integrated.obj = NULL, inFile = NULL, groupBy = "",
                     MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
                    filterTerm = "ENSCAFG", n_feat = 5
                    ){
    
    if(!is.null(file)){
        all.markers <- read.csv(inFile)
    }else if(!is.null(seu.obj)){
        all.markers = FindAllMarkers(seu.obj,
                                     min.pct = MIN_PCT_CELLS_EXPR_GENE,
                                     logfc.threshold = MIN_LOGFOLD_CHANGE,
                                     only.pos = TRUE)
    }

    key.genes <- all.markers[!grepl(paste(filterTerm,collapse="|"), all.markers$gene),] 
    key.genes.sortedByPval = key.genes[order(key.genes$p_val),]

    features <- key.genes.sortedByPval %>%  group_by(cluster) %>% do(head(., n=n_feat))
    features <- as.data.frame(features[!duplicated(features$gene),])

    #features_cnt <- features %>% count(cluster)
    #features_cnt$n <- rev(features_cnt$n)
    features_cnt <- as.data.frame(table(features$cluster))
    features_cnt$n <- rev(features_cnt$Freq)
    
    features_cnt <- features_cnt %>% mutate(endVal = cumsum(n)+0.5, startVal = endVal-n)

    features_cnt$endVal <- rev(features_cnt$endVal)
    features_cnt$startVal <- rev(features_cnt$startVal)
#     features_cnt$cluster <- as.numeric(features_cnt$cluster)
    features_cnt$Var1 <- as.numeric(features_cnt$Var1)

    
    p <- DotPlot(seu.integrated.obj,
             assay = "RNA",
             features = rev(features$gene),
             #col.min = 0,
             group.by = groupBy,
             scale = TRUE) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  labs(size='Percent\nExpression')  +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification='center',
        #legend.spacing.x = unit(1, "cm"),
        legend.key.size = unit(1, "line"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1),
       # ) +
        axis.ticks.x = element_blank()) +
  scale_colour_viridis(option="magma", name='Average\nExpression') +
  guides(size=guide_legend(override.aes = list(shape=21, colour="black", fill="white"),
                           label.position = "bottom")) +
  scale_size(range = c(0.5, 8), limits = c(0, 100)) +
  annotate("rect", xmin = features_cnt$startVal, xmax = features_cnt$endVal, ymin = features_cnt$Var1-0.5, ymax = features_cnt$Var1+0.5, fill = NA, colour = "grey20", size = 0.5) +
  geom_tile(aes(fill = id, x = 0), linewidth = 1, show.legend = FALSE) + 
  geom_tile(aes(fill = id, x = as.numeric(length(unique(features$gene)))+1), linewidth = 1, show.legend = FALSE) + #need to extract data then add colors from colArray to get this corrected
  #scale_fill_manual(values=c("0" = "#CDAD00","1" = "#FFC125", "2" = "#0288D1", "3" = "blue", "4" = "#00CDCD", "5" = "gold", "6" = "brown")) +
  coord_flip() +
  guides(color = guide_colorbar(title = 'Scaled\nExpression')) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_point(aes(x = 0, y=id, size = 100), shape = 21, stroke = 0.75) +
  geom_point(aes(x = as.numeric(length(unique(features$gene)))+1, y=id, size = 100), shape = 21, stroke = 0.75) +
  geom_text(aes(x = 0, label = id), size = 4.5) +
  geom_text(aes(x = as.numeric(length(unique(features$gene)))+1, label = id), size = 4.5)
    
    return(p)
}

############ stackedBar ############
#' Generate a stacked bar graph from a Seurat object
#'
#' @param seu.obj Seurat object to process
#' @param groupBy a valid metadata slot that the bar chart will be colorized by
#' @param clusters a valid metadata slot that will be the rows of the barchart
#' @param downSampleBy depreciated, will not be used. If you want to down sample, do so before hand
#' @param flipOrder logical if TRUE the ordering of the groupBy colorization will be filpped
#'        If using flipOrder, be careful when attempting to re-color the plot
#'
#' @return stack bar graph showing proportion of "groupBy" within each "clusters"
#'
#' @examples stackedBar(seu.obj, groupBy = "orig.ident", clusters = "clusterID")
#'
#' @export 

stackedBar <- function(
    seu.obj = NULL, 
    groupBy = "orig.ident", 
    clusters = "clusterID",
    flipOrder = F,
    downSampleBy = NULL
){
    
    #check input
    if(!is.factor(seu.obj@meta.data[[groupBy]]) & !is.factor(seu.obj@meta.data[[clusters]])){
        message(paste("\nINFO: The metadata slot(s)",
                      groupBy, "and/or", clusters, "are/is not a factor data type.",
                      "It is reccomended that the slots be converted to factors with", 
                      "the levels in the order that you want them to be plotted.",
                      "\n\nProceeding with plotting despite not using factors."))
    }
    
    #extract and clean dataframe
    cluster_freq.table <- as.data.frame(
        table(seu.obj@meta.data[[groupBy]], seu.obj@meta.data[[clusters]])) %>% melt()
    cluster_freq.table <- cluster_freq.table[,-3]
    colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")
    
    #convert to proportion within each "clusters" and remove 0 values
    cluster_freq.table <- cluster_freq.table %>% group_by(ClusterID) %>%
    mutate(pct = round(prop.table(Count), 2)) %>% filter(Count !=  0)

    if(flipOrder){
        if(is.factor(seu.obj@meta.data[[groupBy]])){
            cluster_freq.table$Sample <- factor(cluster_freq.table$Sample, 
                                                levels = rev(levels(cluster_freq.table$Sample))
                                               )
        } else{
            message(paste("\nINFO: The", groupBy, "metadata slot is not a factor.",
                          "To use the flipOrder option you must first convert the",
                          "metadata slot to a factor.\n\nIgnoring flipOrder option."
                         )
                   )
        }
    }
    
    #plot the data
    p <- ggplot(cluster_freq.table, aes(x = ClusterID, y = pct, fill = Sample)) +
    geom_bar(stat = "identity", position = "fill", width = 1, colour = "white") +
    theme_classic() +
    theme(title = element_text(size= 14),
          legend.title = element_blank(),
          legend.text = element_text(size= 12),
          legend.position = "top",
          legend.direction = "horizontal",
          plot.title = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key.size = unit(1,"line"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.margin = margin(t = 0, r = 21, b = 0, l = 0, unit = "pt")
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    ylab(label = "Fraction of cells") +
    xlab(label = "Cluster ID") + 
    scale_x_discrete(expand=c(0,0))
    return(p)
}

############ runMonocle ############
                              #work in progress
runMonocle <- function(seu.obj = NULL
                    ){
    
    cds <- as.CellDataSet(seu.obj) # make sure non-normlized data are imported
    
    cds <- clusterCells(cds)
    disp_table <- dispersionTable(cds)
    ordering_genes <- subset(disp_table, mean_expression >= 0.1)
    cds <- setOrderingFilter(cds, ordering_genes)
    cds <- reduceDimension(cds)
    cds <- orderCells(cds)
    
    diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~Media")
    ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
    
    cds <- setOrderingFilter(cds, ordering_genes)
    plot_ordering_genes(cds)
    
    cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
    plot_cell_trajectory(cds, color_by = "Hours")
    
    
    
    
    seu.obj.Export <- exportCDS(cds, 'Seurat')
}

### trajectory plot from slingshot
cleanSling <- function(plot = NULL, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, branchData = NULL
                  ) {
    
    pi <- formatUMAP(plot)

    #extract label coords and colors
    g <- ggplot_build(pi)
    
    labCords <- as.data.frame(g$data[2]) #add error if labels are not present
    labCordz <- labCords[,c(1:4,7)]

    colnames(labCordz) <- c("colour", "clusterID", "UMAP1", "UMAP2", "labCol")

    labCordz <- labCordz[order(labCordz$clusterID),]
    labCordz$labCol <- labCol
    
    if(rm.na == T){
        labCordz <- na.omit(labCordz)
    }
    
    #remove old labels
    pi$layers[2] <- NULL

    df.list <- list()
    cnt = 0 
    for (lin in 1:length(branchData)) {
        cnt <- cnt+1
        df <- labCordz[labCordz$clusterID %in% branchData[lin][[1]],]
        #df <- df[order(match(branchData[lin][[1]], labCordz$clusterID)),]
        df <- left_join(data.frame(clusterID = branchData[lin][[1]]),  
                         df,
                         by = "clusterID")
        df$lineage <- names(branchData[lin])
        df <- df %>% mutate(UMAP1_end = lead(UMAP1),
                            UMAP2_end = lead(UMAP2)
                           )

        df.list[[cnt]] <- df
    }
    lineData <- do.call(rbind, df.list)
    
    #add lines and dots to the stripped plot to create final image
    plot <- pi + geom_segment(data = lineData, aes(x = UMAP1,  xend = UMAP1_end,
                                                   y = UMAP2,  yend = UMAP2_end),
                              linewidth = 1) + 
    geom_point(data = lineData, aes(x = UMAP1, y = UMAP2), size = 3) + 
    geom_point(data = lineData[1,], aes(x = UMAP1, y = UMAP2), size = 1.5, colour = "green")

    
    return(plot)
}

############ createCIBERsort ############
#updated 07.14.23 to calc CPM instead of using log normalizaed data
createCIBERsort <- function(seu.obj = NULL, groupBy = NULL, downSample = F, outDir = "", outName = "", normMethod = "log"
                    ){

    Idents(seu.obj) <- groupBy

    if(downSample){
        #randomly downsample the subset data
        set.seed(12)
        seu.obj <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data[[groupBy]])))
        if(normMethod == "log"){
            seu.obj <- NormalizeData(seu.obj)
        }
        
    }
    
    #extract required data for pseudobulk conversion
    if(normMethod == "log"){
        mat <- seu.obj@assays$RNA@data
    }else{
        mat <- seu.obj@assays$RNA@counts
    }
    bioRep <- as.data.frame(seu.obj@meta.data[[groupBy]])
    colnames(bioRep) <- "cellSource"
    row.names(bioRep) <- colnames(seu.obj)
    bioRep$cellSource <- as.factor(bioRep$cellSource)
    
    #use custom function to convert to pseudobulk
    pbj <- getPb(mat, bioRep)

    if(normMethod == "cpm"){
        pbj <- t(t(pbj)/colSums(pbj))*1e6
    }
    
    outfile <- paste0(outDir,outName,"_ciberSort_matrix.csv")
    write.csv(pbj, file = outfile,quote=T)
    
}

############ sankeyPlot ############

sankeyPlot <- function(seu_obj = NULL, new.ident = NULL, old.ident = "clusterID", old.colorz = NULL,
                       new.colorz = NULL, old.labCol = NULL, new.labCol = NULL, flowCol = "grey"
                    ){
    
    Idents(seu_obj) <- new.ident

    #get node data
    new <- levels(seu_obj@active.ident)
    new <- paste0("S", new)
    nodeNum <- length(unique(seu.obj@meta.data[[old.ident]])) + length(levels(seu_obj@active.ident)) - 1
    nodes <- data.frame(node = c(0:nodeNum), 
                        name = c(as.character(sort(as.numeric(unique(seu.obj@meta.data[[old.ident]])))), new))

    #exctra data for plotting
    seu_obj_data <- as.data.frame(seu_obj@meta.data)
    seu_obj_data$barcode <- rownames(seu_obj_data)
    seu_obj_data <- seu_obj_data[c("barcode", old.ident,new.ident)] 
    colnames(seu_obj_data) <- c("barcode", "Initial","SubCluster")
    seu_obj_data[c("Initial","SubCluster")] <- lapply(seu_obj_data[c("Initial","SubCluster")], as.character)

    #use ggsankey package to get data in correct format
    data <- seu_obj_data %>% make_long(Initial, SubCluster)

    #prefix sub clusters with "S"
    data$next_node <- ifelse(!is.na(data$next_node),paste0("S",data$next_node),NA)
    data <- data %>% mutate(node = ifelse(x == "SubCluster",paste0("S",node, sep=""),node))

    #order the groups so they are colored appropriately
    data$node <- factor(data$node, levels = nodes$name)

    #make the plot
    p <- ggplot(data, aes(x = x,
                          next_x = next_x,
                          node = node,
                          next_node = next_node,
                          fill = factor(node),
                          label = node)) + geom_sankey() +
    geom_sankey(flow.alpha = 0.5, node.fill = c(old.colorz,new.colorz), flow.fill = flowCol) + 
    geom_sankey_label(size = 3.5, color = 1, fill = "white") +
    theme_sankey(base_size = 16) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1),
          plot.margin = margin(t = 7, r = 14, b = 7, l = 7, unit = "pt")) + scale_x_discrete(expand = c(0, 0))
    
    return(p)
}

############ dotPlotBY_TYPE ############

dotPlotBY_TYPE <- function(seu_obj = NULL, pwdTOvilnCSVoutput = NULL, groupBy = NULL, namedCols = NULL, database = "clfamiliaris_gene_ensembl", exlcude = "", boxColor = "black"
                          ){
    
    df <- read.csv(pwdTOvilnCSVoutput)
    
    ensembl <-  biomaRt::useMart("ensembl", dataset = database)
    
    cl_go_anno <- biomaRt::getBM(
        attributes = c("external_gene_name", "description", "go_id", "name_1006", "namespace_1003", "definition_1006"),
        filters = c("external_gene_name"), 
        values = unique(as.character(df$gene)),
        mart = ensembl) %>%
    dplyr::filter(namespace_1003 %in% c("cellular_component", "molecular_function"))
    
    surface <- sort(unique(cl_go_anno$external_gene_name[grep("component of membrane", cl_go_anno$name_1006)]))
    surface.neg <- sort(unique(cl_go_anno$external_gene_name[grep("mitochondrial|integral component of Golgi membrane|intracellular membrane-bounded organelle|nuclear membrane", cl_go_anno$name_1006)]))
    surface <- surface[!(surface %in% surface.neg)]
    
    secreted <- sort(unique(cl_go_anno$external_gene_name[
        c(grep("extracellular space", cl_go_anno$name_1006),
          grep("granzyme", cl_go_anno$description))]))
    
    tf <- sort(unique(cl_go_anno$external_gene_name[
        grep("transcription factor activity|transcription activator activity|transcription repressor activity|transcription coactivator activity|DNA binding|DNA-binding",cl_go_anno$name_1006)]))
    
    intersect(surface, secreted)
    intersect(surface, tf)
    intersect(secreted, tf)
    
    surface <- surface[!(surface %in% tf)]
    secreted <- secreted[(!(secreted %in% surface)) & (!(secreted %in% tf))]
    neg <- sort(unique(df$gene[!(df$gene %in% c(surface, secreted, tf))]))
    
    anno.df <- data.frame(
        gene = c(surface, secreted, tf, neg),
        anno = c(rep("cell surface", length(surface)),
                 rep("secreted", length(secreted)),
                 rep("transcription factor", length(tf)),
                 rep("", length(neg)))) %>%
    dplyr::filter(!(gene %in% exlcude))
    
    df.plot <- df %>% dplyr::left_join(anno.df, by = c("gene")) %>%
    dplyr::select(gene, cluster, anno, avg_log2FC) 
    
    df.plot$anno <- factor(df.plot$anno, levels = c("cell surface", "secreted", "transcription factor", ""))

    if(is.null(groupBy)){
        Idents(seu.obj) <- "clusterID_sub"
    }
    
    p <- mapply(function(i, j){
        x <- df.plot %>%
        distinct(gene,.keep_all= TRUE) %>%
        dplyr::filter(anno == i) %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(5, avg_log2FC) %>%
        #dplyr::group_by(cluster) %>%
        dplyr::arrange(gene, .by_group = T) %>% 
        dplyr::group_by(cluster) %>%
        mutate(n = n())
        
        features_cnt <- x[!duplicated(x$cluster),]#  x %>% unique(cluster)
        features_cnt$n <- rev(features_cnt$n)
        
        features_cnt <- as.data.frame(features_cnt) %>% mutate(endVal = cumsum(n)+0.5, startVal = endVal-n)
        
        features_cnt$endVal <- rev(features_cnt$endVal)
        features_cnt$startVal <- rev(features_cnt$startVal)
        features_cnt$cluster <- as.numeric(features_cnt$cluster)+1 
        
        p <- DotPlot(seu.obj, features = rev(x$gene), group.by = groupBy) +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        labs(size='Percent\nExpression', title = j)  +
        theme(axis.text.x = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.justification='center',
              legend.key.size = unit(1, "line"),
              panel.border = element_rect(color = "black",
                                          fill = NA,
                                          size = 1),
              axis.ticks.x = element_blank()) +
        scale_colour_continuous_diverging(palette = 'Blue-Red', mid = 0) + 
        guides(size=guide_legend(override.aes = list(shape=21, 
                                                     colour="black", fill="white"),
                                 label.position = "bottom")) +
        scale_size(range = c(0.5, 8), limits = c(0, 100)) +
        annotate("rect", xmin = features_cnt$startVal, xmax = features_cnt$endVal, ymin = features_cnt$cluster-0.5, ymax = features_cnt$cluster+0.5, fill = NA, colour = boxColor, size = 1) +
        geom_tile(aes(fill = id, x = 0), size = 1, show.legend = FALSE) + 
        geom_tile(aes(fill = id, x = as.numeric(length(unique(x$gene)))+1), size = 1, show.legend = FALSE) + 
        scale_fill_manual(values=namedCols) + 
        coord_flip() +
        guides(color = guide_colorbar(title = 'Scaled\nExpression')) +
        scale_y_discrete(expand = c(0, 0)) +
        geom_text(aes(x = 0, label = id), size = 4.5) + #hardcoded
        geom_text(aes(x = as.numeric(length(unique(x$gene)))+1, label = id), size = 4.5)
        return(p)
    }, i = c("cell surface", "secreted", "transcription factor", ""),
                j = c("Cell\nsurface","Secreted","Transcription\nfactor","Remaining\ngenes"),
                SIMPLIFY = F,
                USE.NAMES = F)
    
    p[[1]] <- p[[1]] + NoLegend()
    p[[2]] <- p[[2]] + theme(axis.title.y = element_blank())
    p[[3]] <- p[[3]] + theme(axis.title.y = element_blank()) + NoLegend()
    p[[4]] <- p[[4]] + theme(axis.title.y = element_blank()) + NoLegend()
    
    finalPlot <- plot_grid(plotlist = p, nrow = 1, labels = "none", label_size = 8, rel_widths = c(1,1,1,1)) 
    
    return(finalPlot)
}

############ prettyViln ############

prettyViln <- function(plot = NULL, colorData = NULL, nrow = 2, ncol = NULL){
    
    if(is.null(ncol)){
        ncol = ceiling(sqrt(length(plot)))
    }
    
    
    p <- lapply(plot, function(x) x + theme(axis.title = element_blank(),
                                          axis.text = element_blank(),
                                          axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =5 ),
                                          #axis.ticks.x = element_blank(),
                                          axis.ticks = element_blank(),
                                          title = element_text(size = 16),
                                          legend.position = "none",
                                          plot.margin = margin(5.5, 0, 0, 0, "pt")
                                         ) + coord_cartesian(expand = TRUE)
                )
    
    if(!is.null(colorData)){
        g <- ggplot(colArray, aes(x = clusterID_sub, y = 1, fill = clusterID_sub, label = clusterID_sub)) + 
        geom_tile(fill = "transparent") +
        geom_point(shape = 21, stroke = 0.75, size = 11) +
        geom_text(fontface = "bold", size = 6, color = colArray$labCol) + theme_bw(base_size = 12) +
        scale_fill_manual(values = colArray$colour) + scale_y_discrete(expand = c(0, 0)) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_blank(), 
              panel.border = element_blank(),
              plot.background = element_blank(), 
              plot.margin = margin(0, 0, 0, 0, "pt"),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()
             )
    }
        
        patch <- area()
        counter=0
        for (i in 1:nrow) {
            for (x in 1:ncol) {
                counter = counter+1
                if (counter <= ncol*nrow) {
                    patch <- append(patch, area(t = i, l = x, b = i, r = x))
                }
            }
        }
        
    
        pi <- Reduce( `+`, p ) + plot_layout(design = patch, 
                                                             widths = c(rep.int(1, ncol)), 
                                                             heights = c(rep.int(1, nrow))
                                                            ) 
        return(pi)
}

############ prettyVolc ############                                          
prettyVolc <- function(
    plot = NULL, 
    lfcCut = 0.58,
    rightLab = NULL, 
    leftLab = NULL, 
    rightCol = "red", 
    leftCol = "blue", 
    arrowz = T
    ){
    
    p <- plot + scale_x_symmetric(mid = 0) + 
    theme(
        legend.position = c(0.10, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title=element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        plot.title = element_blank()
       ) + 
    {if(arrowz){
        annotate("segment", x = lfcCut*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 xend = c(max(abs(plot$data$log2FoldChange)),-max(abs(plot$data$log2FoldChange)))[1], 
                 yend = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 lineend = "round", linejoin = "bevel", linetype ="solid", colour = rightCol,
                 linewidth = 1, arrow = arrow(length = unit(0.1, "inches"))
                ) 
        }} +
    {if(arrowz){
        annotate("segment", x = -lfcCut*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 xend = c(max(abs(plot$data$log2FoldChange)),-max(abs(plot$data$log2FoldChange)))[2],
                 yend = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 lineend = "round", linejoin = "bevel", linetype ="solid", colour = leftCol,
                 linewidth = 1, arrow = arrow(length = unit(0.1, "inches"))
                )
        }} + 
    {if(!is.null(rightLab)){
        annotate(geom = "text", x = (max(abs(plot$data$log2FoldChange))-lfcCut*1.5)/2+lfcCut*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.09,
                 label = rightLab,
                 hjust = 0.5,
                 size = 5)
        }} + 
    {if(!is.null(leftLab)){
        annotate(geom = "text", x = -(max(abs(plot$data$log2FoldChange))-lfcCut*1.5)/2-lfcCut*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.09,
                 label = leftLab,
                 hjust = 0.5,
                 size = 5)
    }} + NoLegend()
    
    return(p)
}


############ plotGSEA ############
#TO DO:L fix problem with trimTerm
plotGSEA <- function(
    pwdTOgeneList = NULL, 
    geneList = NULL, 
    geneListDwn = NULL,
    category = "C5", 
    species = "dog", 
    upCol = "red", 
    dwnCol = "blue", 
    saveRes = NULL,
    pvalueCutoff = 0.05, 
    subcategory = NULL, 
    termsTOplot = 8, 
    upOnly = F, 
    dwnOnly = F,
    trimTerm = T, 
    size = 4,
    trunkTerm = F,
    lolli = F
    ){
    
    if(!is.null(pwdTOgeneList)){
        geneLists <- read.csv(pwdTOgeneList)
        geneListUp <- geneLists %>% arrange(padj) %>% filter(log2FoldChange > 0) %>% .$gene
        geneListDwn <- geneLists %>% arrange(padj) %>% filter(log2FoldChange < 0) %>% .$gene
        
    }else{
        geneListUp <- geneList
    }
    
    
    can_gene_sets <- as.data.frame(msigdbr(species = species, category = category, subcategory = subcategory))
    msigdbr_list <- split(x = can_gene_sets$gene_symbol, f = can_gene_sets$gs_name)
    datas <- can_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

    if(!dwnOnly){
    enriched_up <- as.data.frame(enricher(gene = geneListUp, TERM2GENE = datas, pvalueCutoff = pvalueCutoff)
                                ) %>% arrange(p.adjust) %>% mutate(x_axis = -log10(p.adjust),
                                                                   direction = "UP")
    } else{
        enriched_up <- data.frame(matrix(ncol = 0, nrow = 0))
    }
    if(!upOnly){
        enriched_dwn <- as.data.frame(enricher(gene = geneListDwn, 
                                               TERM2GENE = datas, pvalueCutoff = pvalueCutoff)
                                     ) %>% arrange(p.adjust) %>% mutate(x_axis = log10(p.adjust),
                                                                        direction = "DOWN")
    
        enriched <- rbind(enriched_up[1:ifelse(nrow(enriched_up) > termsTOplot, termsTOplot, length(enriched_up)),], 
                          enriched_dwn[1:ifelse(nrow(enriched_dwn) > termsTOplot, termsTOplot, length(enriched_up)),])
        enriched$direction <- factor(enriched$direction, levels = c("UP", "DOWN"))
    } else {
        enriched <- enriched_up[1:ifelse(nrow(enriched_up) > termsTOplot, termsTOplot, nrow(enriched_up)),]

    }
    
    if(!is.null(saveRes)){
        if(!exists("enriched_dwn")){
            write.csv(enriched_up, file = saveRes)
        } else {
             write.csv(rbind(enriched_up,enriched_dwn), file = saveRes)
        }
       
    }

    enriched <- enriched %>% arrange(desc(x_axis))
    
    if(trimTerm){
        enriched$Description <- gsub("REACTOME_", "", enriched$Description)
        enriched$Description <- gsub("_AND_", "_&_", enriched$Description)
        enriched$Description <- gsub("OTHER_", "", enriched$Description)
        enriched$Description <- gsub("GOBP_", "", enriched$Description)
        enriched$Description <- gsub("GOCC_", "", enriched$Description)
        enriched$Description <- gsub("GOMF_", "", enriched$Description)
        enriched$Description <- gsub("HP_", "", enriched$Description)
        enriched$Description <- gsub("POSITIVE_", "POS_", enriched$Description)
        enriched$Description <- gsub("ANTIGEN_PROCESSING_&_", "", enriched$Description)
        enriched$Description <- gsub("ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS", "ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION", enriched$Description)
        enriched$Description <- gsub("PROTON_TRANSPORTING_ATP_SYNTHASE_ACTIVITY_ROTATIONAL_MECHANISM", "PROTON_TRANSPORTING_ATP_SYNTHASE_ACTIVITY", enriched$Description)
        if(trunkTerm){
            enriched$Description <- unlist(lapply(enriched$Description, str_trunc, width = 50, side = "right", ellipsis = "*")) ### TO DO: this causes problems if trunc results in the same string -- need to ensure that is handled before adding this back in
        }
    }
    
    #remove NAs and then plot
    enriched <- na.omit(enriched)
    enriched$Count <- as.numeric(enriched$Count)
    if(!lolli){
        p <- ggplot(data=enriched, aes(x=x_axis, y=Description, fill = direction)) +
        geom_bar(stat="identity") +  theme_classic() + scale_y_discrete(limits=rev(enriched$Description)
                                                                       ) + 
        geom_text(aes(x = ifelse(x_axis > 0, -0.05,0.05), label = Description), hjust = ifelse(enriched$x_axis > 0, 1,0), size = size) + 
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.line.y = element_blank(),
              legend.justification = "top"
             ) + coord_cartesian(clip = "off") + scale_x_symmetric(mid = 0, name = "Signed log10(padj)") + NoLegend() + scale_fill_manual(values = c(upCol,dwnCol))
        } else{
            p <- ggplot(data = enriched, aes(x = x_axis, y = Description, colour = direction, size = Count)) +
                geom_segment(aes(x = 0, xend = x_axis, y = Description, yend = Description), color = "black", size = 1) +
                geom_point() +
                theme_classic() +
                geom_vline(xintercept = 0) + 
                scale_y_discrete(limits = rev(enriched$Description)) +
                geom_text(aes(x = ifelse(x_axis > 0, -0.05,0.05), label = Description), 
                          hjust = ifelse(enriched$x_axis > 0, 1,0), size = size, color = "black") + 
                theme(axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.line.y = element_blank(),
                      legend.justification = "top"
                     ) +
                coord_cartesian(clip = "off") + 
                scale_x_symmetric(mid = 0, name = "Signed log10(padj)") + 
                scale_size_continuous(limits = c(1, max(enriched$Count))) +
                scale_colour_manual(values = c(upCol, dwnCol), guide = "none")
    }
    
    return(p)
}

############ crossConditionDEG ############
###TO DO: fix intercettion line
crossConditionDEG <- function(
    DE_res1 = NULL, 
    DE_res2 = NULL, 
    DE_name1 = "Canine", 
    DE_name2 = "Human",
    seu.obj1 = NULL,
    seu.obj2 = NULL,
    nlab = 10, 
    nlab_axis = 2, 
    overlapGenes = NULL, 
    colUp = "red",
    colDwn = "blue", 
    contrast = NULL, 
    seed = 12,
    hjustvar = c(1,0.5,0,0.5,1,0,1,0),
    vjustvar = c(0.5,1,0.5,0,0,1,1,0),
    outDir = "../output/", 
    saveGeneList = F
    ){

    if(!is.null(seu.obj1) & !is.null(seu.obj2)){
        overlapGenes <- intersect(rownames(seu.obj1), rownames(seu.obj2))
    } else(
        overlapGenes <- unique(c(DE_res1$gene, DE_res2$gene))
    )
    
    geneLists.spec1 <- DE_res1 %>% 
        mutate(
            condition = DE_name1,
            signed_pVal_spec1 = ifelse(log2FoldChange > 0, -log10(padj), log10(padj))
        ) %>% 
        filter(gene %in% overlapGenes)
    
    geneLists.spec2 <- DE_res2 %>% 
        mutate(
            condition = DE_name2,
            signed_pVal_spec2 = ifelse(log2FoldChange > 0, -log10(padj), log10(padj))
        ) %>% 
        filter(gene %in% overlapGenes)
    
    geneList <- as.data.frame(unique(c(geneLists.spec1$gene, geneLists.spec2$gene)))
    colnames(geneList) <- "gene"
    
    geneList <- geneList %>% 
        left_join(geneLists.spec1[,c("gene","condition","signed_pVal_spec1")], by = "gene") %>% 
        left_join(geneLists.spec2[,c("gene","signed_pVal_spec2")], by = "gene") %>% 
        replace(is.na(.), 0) %>% 
        mutate(
            condition = ifelse(condition == "0", DE_name2,DE_name1),
            condition = ifelse(signed_pVal_spec2 != "0" & signed_pVal_spec1 != "0", "Both",DE_name1),
            sig = ifelse(
                signed_pVal_spec2 * signed_pVal_spec1 == 0,
                signed_pVal_spec2 + signed_pVal_spec1, 
                signed_pVal_spec2 * signed_pVal_spec1
            ),
            direction = case_when(
                signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 > 0 & signed_pVal_spec2 > 0 ~ "up",
                signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 < 0 &  signed_pVal_spec2 < 0 ~ "dwn",
                signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec1 > 0 ~ "axis1",
                signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec2 > 0 ~ "axis2",
                signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec1 < 0 ~ "axis3",
                signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec2 < 0 ~ "axis4",                             
                signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 > 0 & signed_pVal_spec2 < 0 ~ "conflict1",
                signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 < 0 & signed_pVal_spec2 > 0 ~ "conflict2"
            ),
            label = ifelse(sig != 0, gene, NA)
        ) %>% 
        arrange(-abs(sig)) %>% 
        group_by(direction) %>%
        mutate(
            lab = ifelse(row_number() <= nlab & signed_pVal_spec2 * signed_pVal_spec1 != 0, 
                         gene, 
                         ifelse(signed_pVal_spec2 * signed_pVal_spec1 == 0 & row_number() <= nlab_axis, 
                                gene, NA
                               )
                        ), 
            lab_col = case_when(
                direction == "up" ~ colUp,
                direction == "dwn" ~ colDwn,
                direction == "axis1" ~ "black",
                direction == "axis2" ~ "black",
                direction == "axis3" ~ "black",
                direction == "axis4" ~ "black",
                direction == "conflict1" ~"hotpink",
                direction == "conflict2" ~"hotpink"
            )
        ) 
    
    anno.df <- as.data.frame(list(
        direction = c("axis1", "axis2", "axis3", "axis4", "conflict1", "conflict2", "up", "dwn"),
        xpos = c(Inf, 0, -Inf, 0, Inf, -Inf, Inf, -Inf),
        ypos =  c(0, Inf, 0, -Inf, -Inf, Inf, Inf, -Inf),
        hjustvar = hjustvar,
        vjustvar = vjustvar,
        colz = c("black", "black", "black", "black", "hotpink", "hotpink", colUp, colDwn)
    ))

    annotationz <- geneList %>% 
        group_by(direction) %>% 
        summarize(cntz = n()) %>% 
        left_join(anno.df, by = "direction")

    if(saveGeneList){
        write.csv(geneList, file = paste0(outDir,DE_name1,"_v_", DE_name2, "_", contrast[1],"_",contrast[2],".csv"), row.names = F)
    }
    
    p <- ggplot(geneList, aes(x = signed_pVal_spec1, y = signed_pVal_spec2)) + 
        geom_hline(yintercept = 0, linetype = 2) +
        geom_vline(xintercept = 0, linetype = 2) + 
        geom_point() + 
        geom_label(data = annotationz, aes(x = xpos,y = ypos, hjust = hjustvar, 
                                           vjust = vjustvar, label = cntz), 
                   color = annotationz$colz) +
        geom_label_repel(max.overlaps = Inf, size = 2, label = geneList$lab, 
                         color = geneList$lab_col, show.legend = F, seed = seed,
                         alpha = 0.8
                        ) + 
        theme_classic() + 
        scale_x_symmetric(mid = 0) + 
        scale_y_symmetric(mid = 0) + 
        labs(
            x = paste0(DE_name1, " " ,contrast[1]," vs ",contrast[2]," DE signed log10(p.adj)"), 
            y = paste0(DE_name2, " ",contrast[1]," vs ",contrast[2]," DE signed log10(p.adj)"), 
            title = paste0(DE_name1, " and ", DE_name2," " ,contrast[1]," vs ",contrast[2])
        ) + 
        coord_cartesian(clip = 'off') + 
        theme(plot.title = element_text(size = 12, hjust = 0.5))

    return(p)
    
}
    

############ skewPlot ############
skewPlot <- function(seu.obj = seu.obj
                    ){
    
    pct.df <- table(seu.obj$majorID_sub, seu.obj$name) %>% melt() %>% group_by(Var.2) %>% mutate(samN = sum(value))
    pct.df$Var.1 <- as.factor(pct.df$Var.1)
    
    pct.df$pct <- pct.df$value/pct.df$samN*100
    pct.df$cellSource <- ifelse(grepl("tils",pct.df$Var.2),"TILs","Blood")

    statz <- compare_means(pct ~ cellSource, group.by = "Var.1", pct.df)

    status <- pct.df %>% group_by(Var.1,cellSource) %>% summarize(med = median(pct)) %>% spread(cellSource,med) %>% left_join(statz[c("Var.1","p.adj")], by = "Var.1") %>% mutate(lfc = abs(log2(TILs/Blood)),
    lab_col=case_when(lfc >= 3 & p.adj < 0.05 ~ "Unique",
                      p.adj < 0.05 & lfc < 3 ~ "Skewed",
                      p.adj > 0.05 ~ "Common")) %>% select(Var.1,lab_col) %>% as.data.frame()

    pct.df <- pct.df %>% left_join(status, by = "Var.1")

    pct.df$lab_col <- as.factor(pct.df$lab_col)
    
    p <- ggplot(pct.df, aes(x = Var.1, y = pct, fill = lab_col)) +
    stat_summary(fun = mean, geom = "bar", width = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar",width = 0) + 
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) + 
    labs(x = "Cell type", y = "Percent total") +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.background = element_blank()
  ) + facet_wrap("cellSource", nrow = 2) + scale_fill_manual(labels = c("Common","Skewed","Unique"),
                                                           values = brewer.pal(n = 3, name = "Dark2"),
                                                           name = "Uniqueness\nclassification")
    return(p)
}

############ ExportToCB_cus ############

#need to update to get barcode in first column
ExportToCB_cus <- function(seu.obj = seu.obj, dataset.name = "", outDir = "./output/", markers = NULL, reduction = "umap", test = F, skipEXPR = F, colsTOkeep = NULL,
                           feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                          
                          ){
    
    outDir <- paste0(outDir,dataset.name, "/")
    dir.create(outDir)
    
    if(test){seu.obj <- subset(x = seu.obj, downsample = 500)}

    
    meta <- seu.obj@meta.data
    meta <- rownames_to_column(meta, "barcode")
    if(!is.null(colsTOkeep)){
        meta <- meta[ ,colnames(meta) %in% c("barcode",colsTOkeep)]
    }
    write.table(meta, paste0(outDir,"meta.tsv"), quote=FALSE, sep='\t', row.names = F)
    
    quickGenes <- as.data.frame(feats[feats %in% rownames(seu.obj)])
    colnames(quickGenes) <- "symbol"
    write.csv(quickGenes, paste0(outDir,"quickGenes.csv"), quote=FALSE, row.names = F)
    
    markers <- read.csv(markers)
    markers <- markers[c("cluster","gene","avg_log2FC","p_val_adj")]
    #colnames(markers)[3,4] <- c("avg_diff","p_val")
    write.table(markers,paste0(outDir,"markers.tsv"), quote=FALSE, sep='\t', row.names = F)
    
    if(!skipEXPR){
        data.df <- as.data.frame(seu.obj@assays$RNA@layers$data)
        colnames(data.df) <- colnames(seu.obj)
        rownames(data.df) <- rownames(seu.obj)
        data.df <- rownames_to_column(data.df, "symbol")
        write.table(data.df, paste0(outDir,"exprMatrix.tsv"), quote=FALSE, sep='\t', row.names = F)
        R.utils::gzip(paste0(outDir,"exprMatrix.tsv"), overwrite = TRUE)
    }
    
    data.df <- as.data.frame(seu.obj[[reduction]]@cell.embeddings)
    data.df <- rownames_to_column(data.df, "barcode")
    write.table(data.df,paste0(outDir,"umap.coords.tsv"), quote=FALSE, sep='\t', row.names = F)
    
}

############ convertTOclusID ############
#' convert string cell types to sorted numerical IDs  
#'
#' @param seu.obj Seurat object to process
#' @param metaSlot a valid metadata slot to convert to numerical ID
#' @param newMetaName optional name for the new metadata slot. Default =  paste0(metaSlot, "_clusID")
#' @return Seurat object with new metadata slot with numerical cluster ID with lowest value corresponding to largest cluster
#' @examples 
#' @export convertTOclusID(seu.obj, metaSlot = "majorID")

convertTOclusID <- function(
    seu.obj = seu.obj, 
    metaSlot = NULL,
    newMetaName = NULL
){
    
    #check input
    if(is.null(newMetaName)){
        newMetaName <- paste0(metaSlot, "_clusID")
    }
        
    #extract and order counts by string cluster
    clusterID_final <- table(seu.obj@meta.data[metaSlot]) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
    mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 
   
    #calculate numerical cluster ID (0 for largest cluster)
    newID <- clusterID_final$clusterID_final
    names(newID) <- clusterID_final[ ,metaSlot]
    Idents(seu.obj) <- metaSlot
    seu.obj <- RenameIdents(seu.obj, newID)
    seu.obj@meta.data[newMetaName] <- Idents(seu.obj)
    return(seu.obj)
}


############ cleanMeta ############
#' Remove extra metadata slots to tidy up Seurat object
#'
#' @param seu.obj Seurat object to process
#' @param metaSlot_keep Optional list of column names to keep. If specified this will be the only parameter used
#' @param metaSlot_remove Optional list of column names to discard.
#' @param grepTerms Optional list of terms to search for in metadata columns to then discard.
#'
#' @return Seurat object with cleaned metadata
#' @examples 
#' @export cleanMeta(seu.obj)

cleanMeta <- function(
    seu.obj = seu.obj, 
    metaSlot_keep = NULL,
    metaSlot_remove = NULL,
    grepTerms = c("DF", "pANN", "snn_res")
){
    if(!is.null(metaSlot_keep)){
        seu.obj@meta.data <- seu.obj@meta.data %>% 
            select(all_of(metaSlot_keep))
    } else{
        seu.obj@meta.data <- seu.obj@meta.data[ ,!grepl(paste(grepTerms, 
                                                          collapse = "|"), colnames(seu.obj@meta.data))]
        
        seu.obj@meta.data <- seu.obj@meta.data %>% 
            select(all_of(colnames(seu.obj@meta.data)[!colnames(seu.obj@meta.data) %in% metaSlot_remove]))
    }
    
    
    return(seu.obj)
}

############ daOR ############
#' Use odds ratio to evalute differential abundance at the cluster level
#'
#' @param seu.ob Seurat object to plot from
#' @param groupBy String; Seurat metadata slot to group rows by (typically clusters)
#' @param splitBy String; Seurat metadata slot to split columns by (typically cell source or biological replicates)
#' @param t_map Logical; If TRUE then transpose the plot
#' @param cluster_order List of string; order of clusters to plot
#' 
#' @return A ComplexHeatmap object
#' @examples 
#' @export optionally save heatmap as "[outDir][outDir]_OR_heat.png" file

daOR <- function(
    seu.obj = NULL,
    groupBy = "",
    splitBy = "",
    outName = "",
    outDir = "",
    t_map = FALSE,
    cluster_order = NULL
    ){
    
    freq.df <- table(seu.obj@meta.data[ , groupBy], seu.obj@meta.data[ , splitBy]) %>% melt()
    colnames(freq.df) <- c("Cluster", "Sample", "Freq")

    or.res <- lapply(unique(freq.df$Cluster), function(cluster){
        or.res.tp <- lapply(unique(freq.df$Sample), function(sam){

            tl <- freq.df %>% 
                filter(Cluster == cluster & Sample == sam) %>% 
                pull(Freq)

            tr <- freq.df %>% 
                filter(Cluster == cluster & Sample != sam) %>% 
                pull(Freq) %>% sum()

            bl <- freq.df %>% 
                filter(Cluster != cluster & Sample == sam) %>% 
                pull(Freq) %>% sum()

            br <- freq.df %>% 
                filter(Cluster != cluster & Sample != sam) %>% 
                pull(Freq) %>% sum()

            mat <- matrix(c(tl, tr,
                            bl, br),
                          ncol = 2)

            res.fisher <- fisher.test(mat)
            res.fisher <- data.frame(
                "Cluster" = cluster, 
                "Sample" = sam, 
                "OR" = unname(res.fisher$estimate), 
                "Pvalue" = res.fisher$p.value
            ) %>% mutate(
                "SE_OR" = sqrt(1/tl + 1/tr + 1/bl + 1/br),
                "lower_95CI" = exp(log(OR) - 1.96 * SE_OR),
                "upper_95CI" = exp(log(OR) + 1.96 * SE_OR)
            )

            return(res.fisher)
        })
        or.res.tp <- do.call(rbind, or.res.tp)
        return(or.res.tp)
    })

    or.res.df <- as.data.frame(do.call(rbind, or.res))

    or.res.df <- or.res.df %>% 
        mutate(
            "logOR" = log(OR),
            "Padj" = p.adjust(Pvalue),
            "sig" = case_when(
                Padj <= 0.01 & abs(logOR) >= log(2) ~ "*",
                Padj > 0.01 | abs(logOR) < log(2) ~ ""
            )
        )  
    
    or.res.mat <- or.res.df %>% 
        select(Cluster, Sample, logOR) %>% 
        pivot_wider(names_from = Sample, values_from = logOR) %>%
        as.data.frame() %>%
        column_to_rownames(var = "Cluster") %>%
        as.matrix()

    sig.mat <- or.res.df %>% 
        select(Cluster, Sample, sig) %>% 
        pivot_wider(names_from = Sample, values_from = sig) %>%
        as.data.frame() %>%
        column_to_rownames(var = "Cluster") %>%
        as.matrix()
       
    
    if(!is.null(cluster_order)){
        if(length(intersect(levels(or.res.df$Cluster), cluster_order)) == length(cluster_order)){
            or.res.mat <- or.res.mat[match(cluster_order, rownames(or.res.mat)), ]
            sig.mat <- sig.mat[match(cluster_order, rownames(sig.mat)), ]
        } else{
            message(paste(
                "One or more of the values in the requested order (specified",
                "with cluster_order) was not present in the data, so the order",
                "has not been altered.")
                   )
        }
    }
    
    if(t_map){
        ht <- Heatmap(t(or.res.mat),
                      name = "logOR",
                      cluster_rows = F,
                      col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                      cluster_columns = F,
                      show_column_names = TRUE,
                      row_title_side = "left",
                      column_title_side = "bottom",
                      column_names_rot = 90,
                      column_names_gp = gpar(fontsize = 10),
                      row_names_gp = gpar(fontsize = 10),
                      column_names_centered = FALSE,
                      cell_fun = function(j, i, x, y, w, h, fill) {
                              grid.text(t(sig.mat)[i, j], x, y, gp = gpar(fontsize = 14, col = "grey95"))
                          }
                     )

    } else{

        ht <- Heatmap(or.res.mat,
                      name = "logOR",
                      cluster_rows = F,
                      col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), #viridis(100),
                      cluster_columns = F,
                      show_column_names = TRUE,
                      column_title_side = "bottom",
                      column_names_rot = 90,
                      column_names_centered = TRUE,
                      cell_fun = function(j, i, x, y, w, h, fill) {
                              grid.text(sig.mat[i, j], x, y, gp = gpar(fontsize = 14, col = "grey95"))
                          }
                     )
    }

    return(ht)
}

############ runMilo ############
#' Wrapper to run miloR
#'
#' @param seu.ob Seurat object to plot from
#' @param da_design Data frame; Columns should include "Sample", "Condition", and optionally "Batch"
#' @param groupBy String; Seurat metadata slot to group rows by (typically samples)
#' @param splitBy String; Seurat metadata slot to split columns by (typically cell source or condition)
#' @param outName String; Name to use when saving the file
#' @param subName String; Name to use when saving the file. Defaults to outName if not specified
#' @param blocked Logical; If TRUE then use the "Batch" values in da_design in blocking strategy
#' 
#' @return A list containing a plot and the milo object
#' @examples 
#' @export saves plots as "../output/[outName]/[subName]_milo.png" and "*_NhoodSize.png"

runMilo <- function(
    seu.obj = NULL, 
    da_design = NULL, 
    groupBy = "",
    splitBy = "",
    outName = "",
    embedding = seu.obj@reductions$umap.integrated.harmony@cell.embeddings,
    subName = NULL,
    blocked = FALSE,
    ...
    ){
    
    if(is.null(subName)){
        subName <- outName
    }
    
    seu.obj$Sample <- seu.obj@meta.data[ , groupBy]
    seu.obj$Condition <- seu.obj@meta.data[ , splitBy]

    # Convert from Seurat to sce object
    sce <- as.SingleCellExperiment(seu.obj)
    reducedDim(sce, "PCA") <- seu.obj@reductions$pca@cell.embeddings
    reducedDim(sce, "UMAP") <- embedding

    # Preprocess using miloR to ID neighboorhoods
    milo.obj <- Milo(sce)
    milo.obj$Sample <- droplevels(factor(milo.obj$Sample))
    milo.obj <- buildGraph(milo.obj, k = 30, d = 40)
    milo.obj <- makeNhoods(milo.obj, prop = 0.2, k = 30, d = 40, refined = TRUE, refinement_scheme = "graph")
    p <- plotNhoodSizeHist(milo.obj)
    ggsave(paste0("../output/", outName, "/", outName, "_NhoodSize.png"), width = 7, height = 7)
    
    milo.obj <- countCells(milo.obj, meta.data = data.frame(colData(milo.obj)), samples = "Sample")

    # Set up metadata
    rownames(da_design) <- da_design$Sample
    da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]

    # Calc distance between neighborhoods and test for DA
    milo.obj <- calcNhoodDistance(milo.obj, d = 40)
    if(blocked){
        da_results <- testNhoods(milo.obj, design = ~ Batch + Condition, design.df = da_design)
    } else{
        da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = da_design)
    }
    
    n_diff <- da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.1) %>% nrow()
    if(n_diff == 0){
        message(
            paste(
                "No differentially abundant Nhoods at alpha = 0.1!",
                "The lowest spaitally adjusted P value is:", min(da_results$SpatialFDR),
                "\n Although not reccomended, you can increase alpha to the lowest",
                "spaitally adjusted P value to get an idea of which regoins of the",
                "UMAP are trendy."
            )
        )
    }

    # Plot the results (by neighborhood)
    milo.obj <- buildNhoodGraph(milo.obj)
    p <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods=!is.na(da_results$logFC), ...)
    ggsave(paste("../output/", outName, "/", outName, "_milo.png", sep = ""), width = 6, height = 6)
    return(list(p, milo.obj))
}



############ splitDot ############
#' Takes dge results and plots the data in split dotplots
#'
#' @param seu.ob Seurat object to plot from
#' @param groupBy String; Seurat metadata slot to group rows by (typically clusters)
#' @param splitBy String; Seurat metadata slot to split columns by (typically cell source or biological replicates)
#' @param namedColz Named list; values are colors and names are levels in "splitBy" metadata slot 
#' @param geneList_UP String list; list of enriched genes
#' @param geneList_DWN String list; list of downregulated genes
#' 
#' @return A ggplot2 object
#' @examples 
#' @export

splitDot <- function(
    seu.obj = NULL,
    groupBy = "",
    splitBy = "",
    namedColz = NULL,
    geneList_UP = NULL,
    geneList_DWN = NULL,
    geneColz = c("red", "blue")
    ){

    seu.obj$majorID_sub_split <- factor(paste0(
        as.character(seu.obj@meta.data[ , groupBy]), 
        "-_-", as.character(seu.obj@meta.data[ , splitBy])
    ), levels = paste0(
        rep(levels(seu.obj@meta.data[ , groupBy]), each = length(levels(seu.obj@meta.data[ , splitBy]))), 
        "-_-", levels(seu.obj@meta.data[ , splitBy])
    ))

    p <- DotPlot(seu.obj, assay = "RNA", features = c(geneList_UP, geneList_DWN),
                     group.by = "majorID_sub_split", scale = T
                )

    df <- separate(p$data, col = id, into = c(NA, "Cell source"), sep = "-_-", remove = F)
    
    labz.df <- as.data.frame(list(
        "y_pos" = seq(1.5, length(levels(seu.obj$majorID_sub_split)) - 0.5, 
                      by = length(levels(seu.obj@meta.data[ , splitBy]))),
        "labz" = levels(seu.obj@meta.data[ , groupBy])
    ))

    for (i in 1:nrow(labz.df)){
    p <- p + annotation_custom(
          grob = textGrob(label = labz.df$labz[i], hjust = 1, gp = gpar(cex = 1.5, fontsize = 9)),
          ymin = labz.df$y_pos[i],
          ymax = labz.df$y_pos[i],
          xmin = -0.6,
          xmax = -0.6)
     }
    
    ymax <- seq(2.5, length(levels(seu.obj$majorID_sub_split)) + 0.5, by = 4)
    p <- p + annotate("rect", xmin = 0.5, xmax = length(c(geneList_UP, geneList_DWN)) + 0.5, 
                      ymin = ymax - 2, 
                      ymax = ymax, 
                      alpha = 0.1, fill = "grey50") +
        theme(
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                       colour = c(rep("red", length(geneList_UP)), rep("blue", length(geneList_DWN)))),
            legend.box = "vertical",
            plot.margin = margin(7, 7, 7, 150, "pt"),
            legend.position = "bottom",
            legend.justification='center',
            legend.key = element_rect(fill = 'transparent', colour = NA),
            legend.key.size = unit(1, "line"),
            legend.background = element_rect(fill = 'transparent', colour = NA),
            panel.background = element_rect(fill = 'transparent', colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.border = element_rect(color = "black",
                                        fill = NA,
                                        size = 1)
        ) + 
        scale_colour_viridis(option="magma", name='Average\nexpression', breaks = c(-0.5, 1, 2), 
                             labels = c("-0.5", "1", "2")) +
        guides(color = guide_colorbar(title = 'Scaled\nExpression  '),
               size = guide_legend(override.aes = list(fill = NA, shape = 21),
                                   label.position = "bottom")) + 
        geom_tile(data = df, aes(fill = `Cell source`, x = 0), show.legend = T) + 
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_manual(values = namedColz) + 
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        labs(size='Percent\nexpression') +
        scale_size(range = c(0.5, 8), limits = c(0, 100)) +
        coord_cartesian(clip = 'off')
    
    return(p)
}

netVisual_aggregate_mod <- function(object, signaling, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                layout = c("circle","hierarchy","chord","spatial"),
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                alpha.image = 0.15, point.size = 1.5,
                                group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20,
                                ...) {
  layout <- match.arg(layout)
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }


  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  # prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
  }  else if (layout == "spatial") {
    prob.sum <- apply(prob, c(1,2), sum)
    if (vertex.weight == "incoming"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$indeg
    }
    if (vertex.weight == "outgoing"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$outdeg
    }
    coordinates <- object@images$coordinates
    labels <- object@idents
    gg <- netVisual_spatial(prob.sum, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)

  } else if (layout == "chord") {
    prob.sum <- apply(prob, c(1,2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                        group = group, cell.order = cell.order,
                                        lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                        scale = scale, reduce = reduce,
                                        title.name = NULL, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
  }

  return(gg)

}







sigDEG_heatmap <- function(
    seu.obj = NULL,
    groupBy = "",
    splitBy = "",
    dge_res = NULL,
    lfc_thres = 1,
    forceCleanPlot = F,
    cond_colz = NULL,
    clus_colz = NULL,
    saveName = "../output/",
    ht_height = 4000,
    ht_width = 4000,
    ...
    ){
    
    contrast <- levels(seu.obj@meta.data[ , splitBy])
    nClus <- length(levels(seu.obj@meta.data[ , groupBy]))
    nContrast <- length(contrast)

    seu.obj$type <- factor(paste0(as.character(seu.obj@meta.data[ , groupBy]), "-_-", as.character(seu.obj@meta.data[ , splitBy])),
                           levels = paste0(rep(levels(seu.obj@meta.data[ , groupBy]), 
                                               each = nContrast), "-_-", contrast))
    
    dge_res$gs_base <- factor(dge_res$gs_base, levels = toupper(levels(seu.obj@meta.data[ , groupBy])))
    res.df <- dge_res %>% 
        filter(abs(log2FoldChange) > lfc_thres) %>%
        mutate(
            direction = ifelse(log2FoldChange > 0, "Up", "Down")
        ) %>%
        arrange(gs_base, direction, padj)

    sig.mat <- matrix(nrow = length(unique(res.df$gene)), 
                      ncol = length(levels(seu.obj$type)),
                      dimnames = list(unique(res.df$gene), toupper(levels(seu.obj$type)))
                     )

    for(i in 1:nrow(sig.mat)){
        for(j in 1:ncol(sig.mat)){
            cellType <- strsplit(colnames(sig.mat)[j], "-_-")[[1]][1]
            condition <- strsplit(colnames(sig.mat)[j], "-_-")[[1]][2]
            if(cellType %in% res.df[res.df$gene == rownames(sig.mat)[i], ]$gs_base){
                lfc <- res.df[res.df$gene == rownames(sig.mat)[i] & res.df$gs_base == cellType, ]$log2FoldChange
                if(lfc > 1 & condition == toupper(contrast[2])){
                    sig.mat[i, j] <- "*"
                } else if(lfc < -1 & condition == toupper(contrast[1])){
                    sig.mat[i, j] <- "*"
                } else{
                    sig.mat[i, j] <- ""
                }
            } else{
                sig.mat[i, j] <- ""
            }
        }
    }

    res.df <- res.df[!duplicated(res.df$gene), ]

    #extract metadata and data
    metadata <- seu.obj@meta.data
    expression <- as.data.frame(FetchData(seu.obj, vars = res.df$gene, layer = "data"))
    expression$anno_merge <- seu.obj@meta.data[rownames(expression), ]$type

    #get cell type expression averages - do clus avg expression by sample
    triMean <- function(x, na.rm = TRUE) {
      mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
    }
    
    tuncMean <- function(x){
        mean(x, trim = 0.15, na.rm = TRUE)
    }
    
    seuMean <- function(x){
        mean(expm1(x))
    }

    clusAvg_expression <- expression %>% 
        group_by(anno_merge) %>% 
        summarise(across(where(is.numeric), seuMean)) %>% 
        as.data.frame()
    rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
    clusAvg_expression$anno_merge <- NULL

    #filter matrix for DEGs and scale by row
    mat_scaled <- t(apply(t(clusAvg_expression), 1, scale))
    colnames(mat_scaled) <- rownames(clusAvg_expression)
    mat_scaled <- mat_scaled[ , match(colnames(sig.mat), toupper(colnames(mat_scaled)))]
    mat_scaled <- mat_scaled[match(rownames(sig.mat), rownames(mat_scaled)), ]  

    #avoid confusing plotting
    if(forceCleanPlot){
        df1 <- as.data.frame(mat_scaled) %>%
            rownames_to_column() %>%
            pivot_longer(cols = colnames(.)[2:ncol(.)]) %>%
            mutate(
                join_name = toupper(paste0(rowname, "_", name))
            )
        
        df2 <- as.data.frame(sig.mat) %>%
            rownames_to_column() %>%
            pivot_longer(cols = colnames(.)[2:ncol(.)]) %>%
            mutate(
                join_name = paste0(rowname, "_", name)
            ) %>% 
            left_join(df1[ , c("join_name", "value")], by = "join_name") %>%
            separate(name, sep = "-_-", remove = T, into = c("name2", "condition")) %>%
            select(-join_name) %>%
            pivot_wider(names_from = condition, values_from = c(value.x, value.y)) %>%
            mutate(
                sig = case_when(
                    value.x_PRE == "*" & value.y_PRE > value.y_POST ~ "down",
                    value.x_POST == "*" & value.y_PRE < value.y_POST ~ "up",
                    .default = ""
                )
            ) %>% select(rowname, name2, sig)
        
        for(i in 1:nrow(sig.mat)){
            for(j in 1:ncol(sig.mat)){
                cellType <- strsplit(colnames(sig.mat)[j], "-_-")[[1]][1]
                condition <- strsplit(colnames(sig.mat)[j], "-_-")[[1]][2]
                gene <- rownames(sig.mat)[i]
                
                #if DE detected but means disagree remove sig star
                if(sig.mat[i, j] == "*"){
                    testVal <- df2[df2$rowname == gene & df2$name2 == cellType, ]$sig
                    if(testVal == ""){
                        sig.mat[i, j] <- ""
                    } 
                }
            }
        }

        genesToExclude <- rownames(sig.mat)[rowSums(sig.mat == "*") == 0]
        mat_scaled <- mat_scaled[!rownames(mat_scaled) %in% genesToExclude, ]
        sig.mat <- sig.mat[!rownames(sig.mat) %in% genesToExclude, ]
        
    }
    
    #set annotations
    ha <- HeatmapAnnotation(
        Cluster = factor(paste0("(c", rep(((1:nClus) - 1), each = nContrast), ") ", 
                                unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"-_-")[[1]][1]}))),
                         levels = paste0("(c", ((1:nClus) - 1), ") ", levels(seu.obj@meta.data[ , groupBy]))),
        Condition = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x, "-_-")[[1]][2]})),
        border = TRUE,
        col = list(Cluster = clus_colz, Condition = cond_colz),
        annotation_legend_param = list(
            Cluster = list(direction = "horizontal", nrow = 3),
            Condition = list(direction = "vertical", nrow = 2)
        ),
        show_annotation_name = FALSE,
        show_legend = FALSE
    )

    lgd1 <- Legend(labels = paste0("(c", ((1:nClus) - 1), ") ", levels(seu.obj@meta.data[ , groupBy])),
                   legend_gp = gpar(fill = clus_colz), 
                   title = "Cluster", 
                   direction = "vertical",
                   nrow = 3, 
                   gap = unit(0.6, "cm")
                  )

    lgd2 <- Legend(labels = names(cond_colz),
                   legend_gp = gpar(fill = cond_colz), 
                   title = "Condition", 
                   direction = "vertical",
                   nrow = 2
                  )

    pd <- packLegend(lgd1, lgd2, max_width = unit(45, "cm"), 
        direction = "horizontal", column_gap = unit(5, "mm"), row_gap = unit(0.5, "cm"))

    #plot the data
    ht <- Heatmap(
        mat_scaled,
        name = "mat",
        cluster_rows = F,
        row_title_gp = gpar(fontsize = 24),
        show_row_names = T,
        cluster_columns = F,
        top_annotation = ha,
        show_column_names = F,
        column_split = factor(unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"-_-")[[1]][1]})),
                              levels = levels(seu.obj@meta.data[ , groupBy])),
        row_title = NULL,
        column_title = NULL,
        heatmap_legend_param = list(
            title = "Average expression",
            direction = "horizontal"
            ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sig.mat[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))
        },
        ...
    )

    png(file = saveName, width = ht_width, height = ht_height, res=400)
    par(mfcol=c(1,1))   
    draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "bottom", 
         annotation_legend_list = pd, annotation_legend_side = "top")

    for(i in 1:nClus){
        decorate_annotation("Cluster", slice = i, {
            grid.text(paste0("c", (1:nClus) - 1)[i], just = "center")
        })
    }
    dev.off()
    return(ht)
}


