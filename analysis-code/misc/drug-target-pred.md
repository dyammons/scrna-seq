## Drug target prediction using scRNA datasets

This is an experimental/simplistic approach to predict which cells may be sussceptible to certain drugs. The analysis pulls on a database of genesets that are known or predicted gene targets of various drugs. The gene signatures were made publically avalible through Yoo et al. 2015 and the databse was coined DSigDB (drug signatures database for gene set analysis). The [webpage](http://tanlab.ucdenver.edu/DSigDB) is where the gene lists can be obtained. Some of the terms have been validated in vitro, while others were scraped from papers or otherwise informatically predicted to be targets. Check out the paper/database for more information, but be aware of how the geneset for each term was curated.

> Yoo, M., Shin, J., Kim, J., Ryall, K.A., Lee, K., Lee, S., Jeon, M., Kang, J. and Tan, A.C., 2015. DSigDB: drug signatures database for gene set analysis. Bioinformatics, 31(18), pp.3069-3071.

&nbsp;

NOTE: code not cleaned up
The analysis approach first loads in the dataset, then does some basic data wrangling. This approach requires a preprocessed Seurat object stored as `seu.obj`.
```r
#load in the data (downloaded from site)
drugTargs.df <- read.gmt("./DSigDB_All.gmt")

#set list of terms you wish to query - peruse the dataset to identify valid options
osDrugz <- c("Palbociclib_FDA","Dinaciclib_LINCS","Alisertib_Roche",
             "AZD1152_TTD_00002243","Sorafenib_FDA", "Cabozantinib_FDA","Lenvatinib_FDA",
             "Regorafenib_FDA","Sunitinib_FDA","Sirolimus","Everolimus","PD-153035",
             "MK-1775_TTD_00009340", "Ibrutinib_FDA")

#subet on terms of interest and stash as named list
drugTargs.df <- drugTargs.df[drugTargs.df$term %in% osDrugz,]
drugTargs.df$term <- droplevels(drugTargs.df$term)
modulez <- split(drugTargs.df$gene, drugTargs.df$term)
```
&nbsp;

The next step is to use module scoring to make inferences about cell type sensitivites (alternatively traditional GSEA could be applied).
```r
#use named list ot complete module scoring
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

#OPTIONAL: modify the names to shorten -- MUST ensure order matches; a named list might be better
names(modulez) <- c("Everolimus", "Rapamycin", "Alisertib", "Dinaciclib",
                    "Cabozantinib", "Lenvatinib", "Palbociclib", "Regorafenib",
                    "Sunitinib", "Sorafenib", "Ibrutinib", "Barasertib", "Adavosertib")

#use trick to update module score names
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
```
&nbsp;

The last bit is to visulize the data (will need to `source("customFunctions.R")` or otherwise load in the custom functions for this section).
```r
#set output name if not already set
outName <- "all_cells"

#set order of "features" if desired 
features <- c("Alisertib", "Barasertib", "Dinaciclib", "Palbociclib", 
              "Sorafenib", "Cabozantinib", "Lenvatinib", "Regorafenib", 
              "Sunitinib","Rapamycin", "Everolimus", "Adavosertib")

#plot the data in dotplot using custom function - can reorder y-axis with the yAxis option in the majorDot function
p <- majorDot(seu.obj = seu.obj, groupBy = "tumorO", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              #axis.ticks = element_blank(),
                              #legend.justification = "left",
                              #plot.margin = margin(7, 21, 7, 7, "pt")
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched'))

ggsave(paste("./output/", outName, "/", outName, "_dots_drugTargs.png", sep = ""),width = 6,height=6)

#plot the data in feature plot using custom function - can reorder y-axis with the yAxis option in the majorDot function
p <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18, noLegend = T) & scale_colour_viridis(option="magma", name='Expression')
ggsave(paste("./output/", outName, "/", outName, "_featsPlot.png", sep = ""), width = 15, height = 9)
```
