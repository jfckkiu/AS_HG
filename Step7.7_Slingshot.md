```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(doParallel)
```

    The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
    which was just loaded, will retire in October 2023.
    Please refer to R-spatial evolution reports for details, especially
    https://r-spatial.org/r/2023/05/15/evolution4.html.
    It may be desirable to make the sf package available;
    package maintainers should consider adding sf to Suggests:.
    The sp package is now running under evolution status 2
         (status 2 uses the sf package in place of rgdal)
    
    Attaching SeuratObject
    
    ── [1mAttaching core tidyverse packages[22m ──────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
    [32m✔[39m [34mdplyr    [39m 1.1.3     [32m✔[39m [34mreadr    [39m 2.1.4
    [32m✔[39m [34mforcats  [39m 1.0.0     [32m✔[39m [34mstringr  [39m 1.5.0
    [32m✔[39m [34mggplot2  [39m 3.4.3     [32m✔[39m [34mtibble   [39m 3.2.1
    [32m✔[39m [34mlubridate[39m 1.9.2     [32m✔[39m [34mtidyr    [39m 1.3.0
    [32m✔[39m [34mpurrr    [39m 1.0.2     
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    [36mℹ[39m Use the conflicted package ([3m[34m<http://conflicted.r-lib.org/>[39m[23m) to force all conflicts to become errors
    Loading required package: foreach
    
    
    Attaching package: ‘foreach’
    
    
    The following objects are masked from ‘package:purrr’:
    
        accumulate, when
    
    
    Loading required package: iterators
    
    Loading required package: parallel
    



```R
immune.combined <- readRDS(file =  '~/AS/AS_Mouse/AS_Mouse2/Final_result/Figure6_Monocle_Macro/Step2_Macro_mono_recluster.Rds')
```


```R
immune.combined@meta.data$Annotation_Formal <- fct_recode(immune.combined@meta.data$Annotation_Formal,
                                                      'Inflammatory Macrophages'   = 'IFN induced Macrophages' )
```


```R
p <- DimPlot(immune.combined, group.by = 'Annotation_Formal',pt.size = 0.1,) + scale_color_manual(values=c('Resident like Macrophages' = "#325A9B", 'Inflammatory Macrophages' = '#FA0087', 'IFN induced Macrophages' =  '#006400',
 'Foam cells' = '#DAA520', 'Monocytes' = '#16FF32'))  + labs(title = '') +  theme_void() + NoLegend()
ggsave(p, file =  'Monocle_celltype.pdf', height = 10, width = 10, units = 'cm')
```


```R
immune.combined@meta.data$`integrated_snn_res.0.5`
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'6'</li><li>'5'</li><li>'1'</li><li>'1'</li><li>'4'</li><li>'0'</li><li>'6'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'0'</li><li>'7'</li><li>'6'</li><li>'2'</li><li>'2'</li><li>'0'</li><li>'0'</li><li>'1'</li><li>'6'</li><li>'6'</li><li>'0'</li><li>'1'</li><li>'5'</li><li>'6'</li><li>'5'</li><li>'5'</li><li>'2'</li><li>'6'</li><li>'3'</li><li>'5'</li><li>'0'</li><li>'5'</li><li>'3'</li><li>'2'</li><li>'0'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'5'</li><li>'5'</li><li>'0'</li><li>'3'</li><li>'0'</li><li>'1'</li><li>'0'</li><li>'2'</li><li>'1'</li><li>'1'</li><li>'1'</li><li>'6'</li><li>'0'</li><li>'3'</li><li>'0'</li><li>'0'</li><li>'6'</li><li>'2'</li><li>'4'</li><li>'1'</li><li>'3'</li><li>'1'</li><li>'3'</li><li>'4'</li><li>'0'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'4'</li><li>'2'</li><li>'4'</li><li>'4'</li><li>'0'</li><li>'0'</li><li>'1'</li><li>'2'</li><li>'2'</li><li>'0'</li><li>'6'</li><li>'4'</li><li>'0'</li><li>'6'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'1'</li><li>'1'</li><li>'5'</li><li>'0'</li><li>'0'</li><li>'0'</li><li>'1'</li><li>'5'</li><li>'1'</li><li>'1'</li><li>'7'</li><li>'5'</li><li>'4'</li><li>'0'</li><li>'2'</li><li>'0'</li><li>'0'</li><li>'7'</li><li>'0'</li><li>'0'</li><li>'4'</li><li>'3'</li><li>'7'</li><li>'0'</li><li>'4'</li><li>'5'</li><li>'3'</li><li>'6'</li><li>'0'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'4'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'0'</li><li>'6'</li><li>'5'</li><li>'1'</li><li>'4'</li><li>'0'</li><li>'5'</li><li>'4'</li><li>'4'</li><li>'0'</li><li>'6'</li><li>'2'</li><li>'0'</li><li>'1'</li><li>'3'</li><li>'0'</li><li>'5'</li><li>'0'</li><li>'2'</li><li>'7'</li><li>'6'</li><li>'6'</li><li>'6'</li><li>'4'</li><li>'7'</li><li>'0'</li><li>'2'</li><li>'0'</li><li>'6'</li><li>'3'</li><li>'2'</li><li>'3'</li><li>'0'</li><li>'7'</li><li>'0'</li><li>'5'</li><li>'0'</li><li>'6'</li><li>'5'</li><li>'0'</li><li>'1'</li><li>'4'</li><li>'4'</li><li>'5'</li><li>'4'</li><li>'0'</li><li>'0'</li><li>'1'</li><li>'6'</li><li>'1'</li><li>'5'</li><li>'0'</li><li>'6'</li><li>'0'</li><li>'4'</li><li>'0'</li><li>'2'</li><li>'0'</li><li>'1'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'7'</li><li>'0'</li><li>'7'</li><li>'6'</li><li>'1'</li><li>'5'</li><li>'1'</li><li>'3'</li><li>'5'</li><li>'6'</li><li>'6'</li><li>'1'</li><li>'5'</li><li>'0'</li><li>'2'</li><li>⋯</li><li>'5'</li><li>'3'</li><li>'0'</li><li>'0'</li><li>'4'</li><li>'1'</li><li>'3'</li><li>'0'</li><li>'3'</li><li>'5'</li><li>'4'</li><li>'0'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'6'</li><li>'3'</li><li>'2'</li><li>'3'</li><li>'0'</li><li>'1'</li><li>'0'</li><li>'1'</li><li>'7'</li><li>'0'</li><li>'0'</li><li>'1'</li><li>'4'</li><li>'6'</li><li>'2'</li><li>'0'</li><li>'0'</li><li>'1'</li><li>'6'</li><li>'4'</li><li>'3'</li><li>'1'</li><li>'0'</li><li>'6'</li><li>'4'</li><li>'5'</li><li>'1'</li><li>'3'</li><li>'5'</li><li>'1'</li><li>'2'</li><li>'5'</li><li>'0'</li><li>'4'</li><li>'0'</li><li>'1'</li><li>'4'</li><li>'1'</li><li>'2'</li><li>'2'</li><li>'0'</li><li>'2'</li><li>'1'</li><li>'3'</li><li>'1'</li><li>'0'</li><li>'1'</li><li>'0'</li><li>'2'</li><li>'1'</li><li>'2'</li><li>'5'</li><li>'5'</li><li>'1'</li><li>'2'</li><li>'2'</li><li>'2'</li><li>'1'</li><li>'2'</li><li>'3'</li><li>'0'</li><li>'3'</li><li>'0'</li><li>'2'</li><li>'4'</li><li>'0'</li><li>'5'</li><li>'2'</li><li>'2'</li><li>'2'</li><li>'7'</li><li>'1'</li><li>'2'</li><li>'1'</li><li>'0'</li><li>'2'</li><li>'2'</li><li>'4'</li><li>'6'</li><li>'1'</li><li>'2'</li><li>'3'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'3'</li><li>'1'</li><li>'1'</li><li>'3'</li><li>'1'</li><li>'6'</li><li>'3'</li><li>'1'</li><li>'2'</li><li>'2'</li><li>'2'</li><li>'2'</li><li>'3'</li><li>'1'</li><li>'3'</li><li>'0'</li><li>'1'</li><li>'2'</li><li>'3'</li><li>'0'</li><li>'2'</li><li>'7'</li><li>'7'</li><li>'2'</li><li>'4'</li><li>'2'</li><li>'3'</li><li>'4'</li><li>'0'</li><li>'1'</li><li>'0'</li><li>'4'</li><li>'1'</li><li>'2'</li><li>'5'</li><li>'7'</li><li>'0'</li><li>'6'</li><li>'2'</li><li>'3'</li><li>'1'</li><li>'2'</li><li>'1'</li><li>'2'</li><li>'1'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'2'</li><li>'1'</li><li>'2'</li><li>'1'</li><li>'2'</li><li>'2'</li><li>'4'</li><li>'3'</li><li>'2'</li><li>'0'</li><li>'6'</li><li>'1'</li><li>'1'</li><li>'4'</li><li>'4'</li><li>'3'</li><li>'1'</li><li>'4'</li><li>'4'</li><li>'0'</li><li>'4'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'3'</li><li>'7'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'0'</li><li>'3'</li><li>'0'</li><li>'3'</li><li>'3'</li><li>'2'</li><li>'2'</li><li>'2'</li><li>'1'</li><li>'1'</li><li>'1'</li><li>'0'</li><li>'0'</li><li>'2'</li><li>'0'</li><li>'2'</li><li>'2'</li><li>'4'</li><li>'3'</li><li>'2'</li><li>'6'</li><li>'0'</li></ol>




```R
###subset Seurat related to monocytes
immune.combined <- subset(immune.combined,integrated_snn_res.0.5 %in% c(4,6,7))
immune.combined <- subset(immune.combined,Annotation_Formal %in% c('Inflammatory Macrophages','IFN induced Macrophages','Foam cells', 'Monocytes' ))
```


```R
DefaultAssay(immune.combined) <- "integrated"
```


```R
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 900)
```

    Warning message in FindVariableFeatures.Assay(object = assay.data, selection.method = selection.method, :
    “selection.method set to 'vst' but count slot is empty; will use data slot instead”
    Warning message in eval(predvars, data, env):
    “NaNs produced”
    Warning message in hvf.info$variance.expected[not.const] <- 10^fit$fitted:
    “number of items to replace is not a multiple of replacement length”



```R
features <- VariableFeatures(immune.combined)
length(features)
```


900



```R
#Find cell cycle related features
source("~/钱教授数据分析v2.0/tools/scTools.R")
cycle_related <- extract_cellcycle(immune.combined, features, cores = 1, cutoff = 0.4, assay = "integrated")
length(cycle_related)
features_subtr <- features[!features %in% cycle_related]
length(features_subtr)
```


91



809



```R
#stress score
mm.stree.genes <- c("Cebpd","Nfkbia","Fos","Egr1","Dnajb1","Ier2","Hspa1b","Hspa1a","Junb","Hspa8", 
                  "Dusp1","Hsph1","Mt1","Nr4a1","Hsp90ab1","Ppp1r15a","Jund","Btg2","Id3","Zfp36",
                  "Jun","Fosb","Atf3","Ubc","Socs3","Hspb1","Dnaja1","Cebpb")
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(mm.stree.genes),
                 ctrl = length(mm.stree.genes),
                 name = 'Stress.Score',
                 assay = "integrated")
```

    Warning message:
    “The following features are not present in the object: Cebpd, Hsp90ab1, Dnaja1, not searching for symbol synonyms”



```R
###IFN score
mm.IFN.genes <- read_csv('~/钱教授数据分析v2.0/tools/IFN_Mouse.csv') %>% 
                pull(gene)
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(mm.IFN.genes),
                 ctrl = length(mm.IFN.genes),
                 name = 'IFN.Score',
                 assay = "integrated")
```

    [1mRows: [22m[34m83[39m [1mColumns: [22m[34m1[39m
    
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m (1): gene
    
    
    [36mℹ[39m Use [30m[47m[30m[47m`spec()`[47m[30m[49m[39m to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set [30m[47m[30m[47m`show_col_types = FALSE`[47m[30m[49m[39m to quiet this message.
    
    Warning message:
    “The following features are not present in the object: Il15, Plscr2, Gm11127, H2-M10.3, H2-M1, Cyth1, Pml, Oas1e, Oas1h, Kif5c, Myd88, Ifi44l, Mx2, Il15ra, 1700057G04Rik, Gm7030, H2-M10.4, H2-M10.5, Trex1, Oas1c, Oas1g, Nmi, Casp1, Lmo2, Gmpr, Themis2, Rabgap1l, H2-M10.2, H2-M11, H2-M10.6, Stx11, Trim14, Oas1b, Oas1a, Tdrd7, Rbck1, Trim75, Tymp, Wars, H2-T24, H2-M10.1, H2-M9, 2410137M14Rik, Samhd1, Oas2, Oas1f, Oas1d, Trim21, not searching for symbol synonyms”



```R
#find stress related genes
stress_related <- extract_stress(immune.combined, features_subtr, cores = 1, cutoff = 0.3, assay = "integrated")
length(stress_related)
features_subtr_2 <- features_subtr[!features_subtr %in% stress_related]
length(features_subtr_2)
```


39



770



```R
#find IFN related genes 
IFN_related <- extract_ifn(immune.combined, features_subtr_2, cores = 1, cutoff = 0.2, assay = "integrated")
length(IFN_related)
features_subtr_3 <- features_subtr_2[!features_subtr_2 %in% IFN_related]
length(features_subtr_3)
```


50



720



```R
#scale intergratedData
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
```


```R
###run without extract genes 
immune.combined <- RunPCA(immune.combined, npcs = 50, do.print = TRUE, verbose = FALSE, pcs.print = 1:5,genes.print = 5,set.seed=666, features = features_subtr_3)
```


```R
VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca")
DimPlot(immune.combined, reduction = "pca")
DimHeatmap(immune.combined, dims = 1:30, cells = 500, balanced = TRUE)
#immune.combined <- JackStraw(immune.combined,dims=50,num.replicate = 50)
#immune.combined <- ScoreJackStraw(immune.combined, dims = 1:50)
#JackStrawPlot(immune.combined, dims = 1:50)
ElbowPlot(immune.combined,ndims=50)
```


    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_17_0.png)
    



    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_17_1.png)
    



    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_17_2.png)
    



    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_17_3.png)
    



```R
#Run PCA, DIM 30
dims.use <- 10
immune.combined <- FindNeighbors(immune.combined, dims = 1:dims.use)
#res.use <- c(seq(0.1,3,by=0.2),5)
#res.use
#immune.combined <- FindClusters(object =immune.combined,resolution = res.use)
immune.combined  <- RunUMAP(immune.combined,reduction.key = "UMAP_", seed.use=666,verbose = TRUE,dims  = 1:dims.use)
```

    Computing nearest neighbor graph
    
    Computing SNN
    
    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    16:49:13 UMAP embedding parameters a = 0.9922 b = 1.112
    
    16:49:13 Read 1655 rows and found 10 numeric columns
    
    16:49:13 Using Annoy for neighbor search, n_neighbors = 30
    
    16:49:13 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    |
    
    16:49:13 Writing NN index file to temp file /tmp/RtmpNQp4o0/file40c2b70da0d
    
    16:49:13 Searching Annoy index using 1 thread, search_k = 3000
    
    16:49:14 Annoy recall = 100%
    
    16:49:14 Commencing smooth kNN distance calibration using 1 thread
    
    16:49:14 Initializing from normalized Laplacian + noise
    
    16:49:14 Commencing optimization for 500 epochs, with 64122 positive edges
    
    16:49:16 Optimization finished
    



```R
DimPlot(immune.combined,label=T, group.by ='Annotation_Brief') 
```


    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_19_0.png)
    



```R
library(scater)
library(slingshot)
library(scales)
library(RColorBrewer)
```

    Loading required package: SingleCellExperiment
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following object is masked from ‘package:dplyr’:
    
        count
    
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: GenomicRanges
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:lubridate’:
    
        intersect, setdiff, union
    
    
    The following objects are masked from ‘package:dplyr’:
    
        combine, intersect, setdiff, union
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
        table, tapply, union, unique, unsplit, which.max, which.min
    
    
    Loading required package: S4Vectors
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following objects are masked from ‘package:lubridate’:
    
        second, second<-
    
    
    The following objects are masked from ‘package:dplyr’:
    
        first, rename
    
    
    The following object is masked from ‘package:tidyr’:
    
        expand
    
    
    The following object is masked from ‘package:utils’:
    
        findMatches
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    
    Attaching package: ‘IRanges’
    
    
    The following object is masked from ‘package:lubridate’:
    
        %within%
    
    
    The following objects are masked from ‘package:dplyr’:
    
        collapse, desc, slice
    
    
    The following object is masked from ‘package:purrr’:
    
        reduce
    
    
    Loading required package: GenomeInfoDb
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: ‘Biobase’
    
    
    The following object is masked from ‘package:MatrixGenerics’:
    
        rowMedians
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        anyMissing, rowMedians
    
    
    
    Attaching package: ‘SummarizedExperiment’
    
    
    The following object is masked from ‘package:SeuratObject’:
    
        Assays
    
    
    The following object is masked from ‘package:Seurat’:
    
        Assays
    
    
    Loading required package: scuttle
    
    Loading required package: princurve
    
    Loading required package: TrajectoryUtils
    
    
    Attaching package: ‘scales’
    
    
    The following object is masked from ‘package:purrr’:
    
        discard
    
    
    The following object is masked from ‘package:readr’:
    
        col_factor
    
    



```R
package.version('slingshot')
```


'2.8.0'



```R
#slingshot重新编码
Idents(immune.combined)  <- "Annotation_Formal"
cl <- Idents(object = immune.combined)
rd <- Embeddings(immune.combined,reduction = "umap")
```


```R
lineage <- getLineages(rd, cl,start.clus = "Monocytes")
ss  <- getCurves(lineage, shrink = TRUE,
       extend = "y", reweight = TRUE, reassign = TRUE, thresh = 0.001,
       maxit = 15, stretch = 2, approx_points = FALSE,
       smoother = "smooth.spline", shrink.method = "cosine",
       allow.breaks = TRUE)
```

    Using full covariance matrix
    



```R
immune.combined@misc[['slingshot']] <-  list("dr"="pca","data"=ss)
```


```R
pst <- slingPseudotime(immune.combined@misc[['slingshot']]$data)
immune.combined@meta.data[,colnames(pst)] <- as.data.frame(pst)
for(cv in colnames(pst)) {
  maxtime <- max(immune.combined@meta.data[,cv], na.rm = T)
  cv_pct  <- paste0(cv,"_pct")
  immune.combined@meta.data[,cv_pct]  <- (immune.combined@meta.data[,cv] / maxtime) * 100
}
```


```R
source("~/钱教授数据分析v2.0/tools/scTools.R")
```


```R
p1 <- plotPseudoTime(immune.combined, group.by = "Annotation_Formal", reduction = "umap",dims = 1:2, pt.size = 0.5) +
  theme_void() + NoLegend() +scale_color_manual(values=c('Resident like Macrophages' = "#325A9B", 'Inflammatory Macrophages' = '#FA0087', 'IFN induced Macrophages' =  '#006400',
 'Foam cells' = '#DAA520', 'Monocytes' = '#16FF32'))
p1
```


    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_27_0.png)
    



```R
p0 <- plotPseudoTime(immune.combined, group.by = "Annotation_Formal", reduction = "umap",dims = 1:2, pt.size = 0) +
  theme_void() + NoLegend() +scale_color_manual(values=c('Resident like Macrophages' = "#325A9B", 'Inflammatory Macrophages' = '#FA0087', 'IFN induced Macrophages' =  '#006400',
 'Foam cells' = '#DAA520', 'Monocytes' = '#16FF32'))
p0
ggsave(p0, file = 'slingshot_celltype.pdf',height = 4 , width = 4, units = 'cm')
```


    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_28_0.png)
    



```R
p3 <- plotPseudoTime(immune.combined, group.by = "group", reduction = "umap",dims = 1:2, pt.size = 0.5) +
  theme_void() + NoLegend() 
p3
ggsave(p3, file = 'slingshot_group.pdf',height = 4 , width = 4, units = 'cm')
```


    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_29_0.png)
    



```R
p2 <- FeaturePlot(immune.combined, pt.size = 0.1, features = 'curve1_pct') +viridis::scale_colour_viridis(option = 'C') + theme_void() + labs(title = '') + NoLegend()
p2
ggsave(p2, file = 'slingshot_pseudotime.pdf',height = 4 , width = 4, units = 'cm')
```

    Scale for 'colour' is already present. Adding another scale for 'colour',
    which will replace the existing scale.
    



    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_30_1.png)
    



```R
library(patchwork)
p <- p0 + p2 + plot_annotation(tag_levels  = 'A')
ggsave(p, file = 'integrated_slingshot.pdf', height = 5, width = 8, units = 'cm')
```


```R
saveRDS(immune.combined, file = 'Step4_Macro_pseudotime.Rds')
```


```R
immune.combined@meta.data$curve1_pct
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>83.1698798756554</li><li>80.3188338858554</li><li>100</li><li>11.6681064102858</li><li>93.2459762098385</li><li>100</li><li>82.3388696112757</li><li>87.4533139484008</li><li>87.8215988564925</li><li>99.9630138648991</li><li>100</li><li>36.356744047555</li><li>66.3099187293569</li><li>43.9265934478355</li><li>31.4300970590558</li><li>100</li><li>57.5588580510426</li><li>100</li><li>3.92503577100484</li><li>79.1498549095412</li><li>8.47376568603206</li><li>33.322037478411</li><li>9.6640271347349</li><li>57.4042136624184</li><li>99.154191345033</li><li>36.7253948084532</li><li>84.079723287497</li><li>46.8514275688472</li><li>41.1786150660419</li><li>68.2813339497998</li><li>92.6948206402077</li><li>0</li><li>26.2920568127165</li><li>84.9666222030216</li><li>80.771625336452</li><li>78.0041491719995</li><li>11.9674889405371</li><li>85.5526386580926</li><li>12.3302013747824</li><li>80.9807747555687</li><li>70.9964816161365</li><li>59.0962129612718</li><li>34.682685640015</li><li>88.3375742904952</li><li>44.6980493883841</li><li>56.4506535792369</li><li>3.08941298689614</li><li>8.30933913080583</li><li>96.8721437553988</li><li>97.0991979844668</li><li>100</li><li>72.6695211478543</li><li>49.686062716898</li><li>99.5762014180478</li><li>95.79517756355</li><li>83.5597797847224</li><li>49.6936166134114</li><li>52.6547670199824</li><li>62.5008550763768</li><li>83.4086449303506</li><li>87.3695532562364</li><li>49.684951968644</li><li>41.0026192405074</li><li>88.316417561159</li><li>83.2913998728962</li><li>88.9762694870278</li><li>74.4272432166907</li><li>6.91447966603846</li><li>48.1731752157026</li><li>56.8115490767512</li><li>54.1391402416942</li><li>68.7916848232689</li><li>70.1197106026948</li><li>86.3161856926318</li><li>41.4717067886695</li><li>32.6790622642978</li><li>30.4994121308763</li><li>47.3407474896837</li><li>34.254880744834</li><li>100</li><li>34.3944625795831</li><li>100</li><li>92.898580773604</li><li>87.6754907914006</li><li>43.3483579256902</li><li>65.4137824173747</li><li>38.651688412116</li><li>70.380432031638</li><li>97.8208735295386</li><li>56.460526284221</li><li>56.0927709760177</li><li>61.0191117439567</li><li>68.6317666204816</li><li>6.45827588883351</li><li>27.3754141765888</li><li>55.5323754133471</li><li>90.3385730902236</li><li>100</li><li>15.6607808554757</li><li>84.8626694600489</li><li>87.2679439531967</li><li>78.9036452120286</li><li>88.7095570835375</li><li>29.1617509212439</li><li>54.6095638543937</li><li>49.1167963978787</li><li>92.1208579971721</li><li>83.2763856887699</li><li>89.3782597512759</li><li>40.195360239989</li><li>100</li><li>56.4880478748234</li><li>32.2142134441348</li><li>84.5248262270645</li><li>32.1375872153344</li><li>94.4046711979252</li><li>82.9754389150742</li><li>62.3243189539365</li><li>77.8009409107128</li><li>41.4900171437311</li><li>46.2216471512193</li><li>58.3191362316641</li><li>6.89194199666754</li><li>85.087824989222</li><li>68.8983585021369</li><li>90.0682096794245</li><li>56.5194193310938</li><li>100</li><li>35.1211473181043</li><li>93.6489939273152</li><li>98.6785720881894</li><li>81.1776758513387</li><li>38.0592330667377</li><li>100</li><li>91.3792095222239</li><li>8.98007968483443</li><li>56.9860087575429</li><li>71.7376115723423</li><li>86.5479623481125</li><li>27.659401342897</li><li>82.9454206332915</li><li>80.1198403393852</li><li>32.4330284579141</li><li>92.180986593263</li><li>86.2912695876817</li><li>47.1465720386757</li><li>30.781574774227</li><li>30.785363618497</li><li>79.553703801738</li><li>97.4154965495357</li><li>38.0792702516152</li><li>75.260689805256</li><li>82.7078238358534</li><li>61.2543125397672</li><li>67.9316776600643</li><li>41.3677910374371</li><li>92.6322328095954</li><li>67.7152653552782</li><li>78.2769162116763</li><li>35.4543455465927</li><li>46.6595765826831</li><li>76.9880735905292</li><li>50.4544141509326</li><li>89.5872358949177</li><li>97.500081041958</li><li>91.1492661795253</li><li>4.48274534877871</li><li>35.0163189314813</li><li>55.4299457659327</li><li>42.8043419066821</li><li>41.3535726748151</li><li>27.3754141765888</li><li>9.05254083427811</li><li>82.3562043482761</li><li>7.63467194186024</li><li>93.1894686821944</li><li>29.2408795870423</li><li>55.7937841958199</li><li>9.13870356162002</li><li>84.6912747969068</li><li>8.47158584732402</li><li>100</li><li>57.350885075014</li><li>57.8434115795778</li><li>42.9652851779233</li><li>81.5473768196149</li><li>74.157688983168</li><li>65.5442161429089</li><li>12.0853242199178</li><li>10.2740590348599</li><li>38.1640983209702</li><li>70.0221838219398</li><li>64.0426195772399</li><li>88.4207087586512</li><li>82.3050610084018</li><li>4.52576146975175</li><li>89.8740791566236</li><li>68.4764903097867</li><li>96.6275555715741</li><li>33.840453613964</li><li>⋯</li><li>58.340558170809</li><li>92.9064813935827</li><li>74.3726004299642</li><li>67.8984761971307</li><li>1.911819452712</li><li>34.9879837117335</li><li>37.9127608788042</li><li>56.1785126158485</li><li>38.0724349967366</li><li>37.5705242489128</li><li>58.9221246439757</li><li>71.172165291334</li><li>85.464745527689</li><li>8.49560336005987</li><li>74.8665393385221</li><li>31.1473944987609</li><li>70.5695565034115</li><li>86.8727213365253</li><li>56.1933379219563</li><li>100</li><li>1.74173723252131</li><li>92.4535657277301</li><li>58.6857622804592</li><li>41.13062169788</li><li>81.9129477362253</li><li>46.4392475074368</li><li>47.4939099267539</li><li>100</li><li>9.51582521575719</li><li>100</li><li>64.3621963657511</li><li>9.01970281363987</li><li>87.6129725605106</li><li>8.11349681886937</li><li>28.8225644802306</li><li>68.5081460881331</li><li>59.808661772232</li><li>41.3974020686003</li><li>54.0823298124223</li><li>87.3105586332588</li><li>82.1868467212816</li><li>51.1458182670535</li><li>100</li><li>53.739829394546</li><li>91.0015151319554</li><li>5.82386681932729</li><li>50.6348262848287</li><li>52.2959675561132</li><li>37.348438961846</li><li>47.5129351541597</li><li>2.27687554011922</li><li>2.76277526859162</li><li>85.491877650742</li><li>90.0959763506478</li><li>1.87766549401851</li><li>60.681422176779</li><li>64.280796385301</li><li>91.4592091066067</li><li>65.45522395173</li><li>99.0469870901043</li><li>91.7662926469417</li><li>92.0009868652792</li><li>29.8821603086906</li><li>48.1330854719368</li><li>82.5235050979516</li><li>88.538057972759</li><li>54.8863953919746</li><li>81.9942903723592</li><li>67.4008087597689</li><li>55.6866338146545</li><li>44.3482750491121</li><li>34.7655894404809</li><li>36.5559271581292</li><li>83.3115371786429</li><li>3.9313657307193</li><li>58.8371274255436</li><li>34.5665226563011</li><li>8.22940939237612</li><li>34.0618600264144</li><li>91.0012424935343</li><li>29.714464321251</li><li>86.6745429885157</li><li>55.9674115562729</li><li>36.3605232404433</li><li>50.4545179362619</li><li>49.5671176934114</li><li>49.8727572610698</li><li>34.7871862653538</li><li>9.9957201898401</li><li>100</li><li>57.9146088675855</li><li>83.4634520319623</li><li>93.2703954098004</li><li>84.5981732317393</li><li>43.1513978185616</li><li>74.4796074759958</li><li>31.785325764594</li><li>100</li><li>57.9430225528165</li><li>58.5329544574328</li><li>64.2096882273303</li><li>57.3382403073744</li><li>6.91238447721564</li><li>65.9477614704797</li><li>100</li><li>86.9215433431163</li><li>8.79314746900186</li><li>2.43384388501453</li><li>100</li><li>52.8325122797074</li><li>83.156282721089</li><li>100</li><li>91.2795753646253</li><li>57.7866959287298</li><li>56.0282954086777</li><li>33.3251813650564</li><li>63.690058979092</li><li>100</li><li>45.9899921499918</li><li>99.8076402382107</li><li>7.87291704335186</li><li>33.3333980463121</li><li>3.4706292680691</li><li>41.0509056578423</li><li>91.5198953953987</li><li>51.2958435078196</li><li>100</li><li>57.5497862810572</li><li>8.41721492559429</li><li>69.8680619071918</li><li>2.66235159677205</li><li>92.1947247237914</li><li>100</li><li>88.3454021276455</li><li>4.38941418452736</li><li>93.5307163078136</li><li>96.2658450730889</li><li>46.1592621100168</li><li>44.3740089358115</li><li>43.5864976286931</li><li>33.6834264473631</li><li>59.5237031487483</li><li>9.33117572107641</li><li>100</li><li>38.7200693450717</li><li>9.55483194586265</li><li>46.2174137510709</li><li>1.57505838413161</li><li>86.0479225838253</li><li>96.7971511305497</li><li>80.4358999551558</li><li>57.9921658043531</li><li>49.5754055420852</li><li>52.067098909806</li><li>29.9850543192581</li><li>4.43173004507058</li><li>48.9402029253544</li><li>96.6970414910458</li><li>7.64657886740779</li><li>78.5522852480455</li><li>89.9637798059107</li><li>1.13070249959413</li><li>92.1446745336699</li><li>34.9165420133292</li><li>46.3817474391342</li><li>1.22771084598734</li><li>100</li><li>26.8986494419651</li><li>8.30602045663152</li><li>40.5278877118138</li><li>93.6694358054856</li><li>90.1305471246797</li><li>4.86904751164255</li><li>35.9760923320122</li><li>26.6263885400373</li><li>92.4243763230338</li><li>87.8360335156958</li><li>35.3848491630268</li><li>27.3754141765888</li><li>47.2657180244636</li><li>52.7354267515984</li><li>7.07410525136364</li><li>44.0750792714911</li><li>89.8829023305063</li><li>88.699602697134</li><li>6.91680787086019</li><li>0.859943014340876</li><li>48.3196673493636</li><li>7.64182555998193</li><li>91.853641344441</li><li>57.629307269003</li><li>86.0232071951995</li><li>37.2120602605815</li><li>70.4528370014987</li><li>75.5531228317976</li><li>55.6390745368363</li><li>52.0659720909821</li><li>7.44225608605195</li><li>41.1627243027755</li><li>97.8654363212237</li></ol>




```R
p2 <- FeaturePlot(immune.combined, features = 'curve1_pct') +viridis::scale_colour_viridis(option = 'C') + theme_void() + labs(title = '')
p2
```

    Scale for 'colour' is already present. Adding another scale for 'colour',
    which will replace the existing scale.
    



    
![png](Step7.7_Slingshot_files/Step7.7_Slingshot_34_1.png)
    

