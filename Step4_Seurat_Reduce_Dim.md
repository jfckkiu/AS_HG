```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(doParallel)
library(RColorBrewer)
library(clustree)
```

    Attaching SeuratObject
    
    Warning message in system("timedatectl", intern = TRUE):
    “running command 'timedatectl' had status 1”
    Registered S3 method overwritten by 'cli':
      method     from         
      print.boxx spatstat.geom
    
    ── [1mAttaching packages[22m ─────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
    
    [32m✔[39m [34mggplot2[39m 3.3.5     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.5     [32m✔[39m [34mdplyr  [39m 1.0.7
    [32m✔[39m [34mtidyr  [39m 1.1.4     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 1.4.0     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    
    Loading required package: foreach
    
    
    Attaching package: ‘foreach’
    
    
    The following objects are masked from ‘package:purrr’:
    
        accumulate, when
    
    
    Loading required package: iterators
    
    Loading required package: parallel
    
    Loading required package: ggraph
    



```R
main_path <- '~/AS/AS_Mouse/AS_Mouse2/'
```


```R
immune.combined <- readRDS(paste0(main_path,'Step3_AddDoublet_Seurat.Rds'))
```


```R
#cell cycle score
load("~/钱教授数据分析v2.0/tools/mouse_cycle.rda")
DefaultAssay(immune.combined) <- 'RNA'
immune.combined <- CellCycleScoring(immune.combined,g2m.features = g2m_genes,s.features = s_genes)
```

    Warning message:
    “The following features are not present in the object: Jpt1, Pimreg, not searching for symbol synonyms”



```R
features <-  rownames(immune.combined@assays$integrated@data)
length(features)
```


3000



```R
#Find cell cycle related features
source("~/钱教授数据分析v2.0/tools/scTools.R")
cycle_related <- extract_cellcycle(immune.combined, features, cores = 1, cutoff = 0.15, assay = "integrated")
length(cycle_related)
features_subtr <- features[!features %in% cycle_related]
length(features_subtr)
```


281



2719



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

    
    [36m──[39m [1m[1mColumn specification[1m[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    cols(
      gene = [31mcol_character()[39m
    )
    
    
    Warning message:
    “The following features are not present in the object: Il15, Plscr2, Gm11127, H2-M10.3, H2-M1, Cyth1, Pml, Oas1e, Oas1h, Kif5c, Myd88, Ifi44l, Mx2, Il15ra, 1700057G04Rik, Gm7030, H2-M10.4, H2-M10.5, Trex1, Oas1c, Oas1g, Nmi, Casp1, Lmo2, Gmpr, Themis2, Rabgap1l, H2-M10.2, H2-M11, H2-M10.6, Stx11, Trim14, Oas1b, Oas1a, Tdrd7, Rbck1, Trim75, Tymp, Wars, H2-T24, H2-M10.1, H2-M9, 2410137M14Rik, Samhd1, Oas2, Oas1f, Oas1d, Trim21, not searching for symbol synonyms”



```R
#find stress related genes
stress_related <- extract_stress(immune.combined, features_subtr, cores = 1, cutoff = 0.3, assay = "integrated")
length(stress_related)
features_subtr_2 <- features_subtr[!features_subtr %in% stress_related]
length(features_subtr_2)
```


50



2669



```R
#find IFN related genes 
IFN_related <- extract_ifn(immune.combined, features_subtr_2, cores = 1, cutoff = 0.3, assay = "integrated")
length(IFN_related)
features_subtr_3 <- features_subtr_2[!features_subtr_2 %in% IFN_related]
length(features_subtr_3)
```


44



2625



```R
#scaleData
DefaultAssay(immune.combined) <-'integrated'
immune.combined<- ScaleData(immune.combined, verbose = FALSE)
```


```R
immune.combined <- RunPCA(immune.combined, npcs = 50, do.print = TRUE, verbose = FALSE, pcs.print = 1:5,genes.print = 5,set.seed=666, features =  features_subtr_3)
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


    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_13_0.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_13_1.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_13_2.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_13_3.png)
    



```R
#PCA SET 20
immune.combined <- FindNeighbors(immune.combined, dims = 1:20)
res.use <- c(seq(0.1,3,by=0.2),3,5)
res.use
  immune.combined <- FindClusters(object =immune.combined,resolution = res.use)
clust.tree.out <- clustree(immune.combined)+
  theme(legend.position = "bottom")+
  scale_color_brewer(palette = "Set1")
clust.tree.out
```

    Computing nearest neighbor graph
    
    Computing SNN
    



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0.1</li><li>0.3</li><li>0.5</li><li>0.7</li><li>0.9</li><li>1.1</li><li>1.3</li><li>1.5</li><li>1.7</li><li>1.9</li><li>2.1</li><li>2.3</li><li>2.5</li><li>2.7</li><li>2.9</li><li>3</li><li>5</li></ol>



    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9305
    Number of communities: 7
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8696
    Number of communities: 12
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8381
    Number of communities: 17
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8179
    Number of communities: 20
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8005
    Number of communities: 21
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7855
    Number of communities: 23
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7720
    Number of communities: 24
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7589
    Number of communities: 25
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7498
    Number of communities: 26
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7380
    Number of communities: 29
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7280
    Number of communities: 32
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7188
    Number of communities: 32
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7109
    Number of communities: 37
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7039
    Number of communities: 38
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.6966
    Number of communities: 39
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.6927
    Number of communities: 41
    Elapsed time: 1 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 11983
    Number of edges: 397268
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.6416
    Number of communities: 56
    Elapsed time: 1 seconds


    Warning message:
    “The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
    Please use the `.add` argument instead.
    [90mThis warning is displayed once every 8 hours.[39m
    [90mCall `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.[39m”
    Warning message in RColorBrewer::brewer.pal(n, pal):
    “n too large, allowed maximum for palette Set1 is 9
    Returning the palette you asked for with that many colors
    ”
    Warning message:
    “Removed 304 rows containing missing values (geom_point).”



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_14_4.png)
    



```R
work_path <- paste0(main_path,'Seurat_result/')
```


```R
ggsave(clust.tree.out, file = paste0(work_path, 'clustree.pdf'), width = 20)
```

    Saving 20 x 6.67 in image
    
    Warning message in RColorBrewer::brewer.pal(n, pal):
    “n too large, allowed maximum for palette Set1 is 9
    Returning the palette you asked for with that many colors
    ”
    Warning message:
    “Removed 304 rows containing missing values (geom_point).”



```R
immune.combined <- RunUMAP(immune.combined,reduction.key = "UMAP_", seed.use=666,dims  = 1:20)
```

    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    20:06:03 UMAP embedding parameters a = 0.9922 b = 1.112
    
    20:06:03 Read 11983 rows and found 20 numeric columns
    
    20:06:03 Using Annoy for neighbor search, n_neighbors = 30
    
    20:06:03 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    20:06:04 Writing NN index file to temp file /tmp/RtmpqmiLVP/file6f354248c7a
    
    20:06:04 Searching Annoy index using 1 thread, search_k = 3000
    
    20:06:07 Annoy recall = 100%
    
    20:06:07 Commencing smooth kNN distance calibration using 1 thread
    
    20:06:08 Found 2 connected components, 
    falling back to 'spca' initialization with init_sdev = 1
    
    20:06:08 Initializing from PCA
    
    20:06:08 PCA: 2 components explained 43.83% variance
    
    20:06:08 Commencing optimization for 200 epochs, with 501400 positive edges
    
    20:06:12 Optimization finished
    



```R
DimPlot(immune.combined, group.by = 'integrated_snn_res.0.5',label = T) +NoLegend()
DimPlot(immune.combined, group.by = 'integrated_snn_res.1.1',label = T) +NoLegend()
DimPlot(immune.combined, group.by = 'integrated_snn_res.3',label = T) +NoLegend()
DimPlot(immune.combined, group.by = 'integrated_snn_res.2.1',label = T) +NoLegend()
```


    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_18_0.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_18_1.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_18_2.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_18_3.png)
    



```R
DimPlot(immune.combined, group.by = 'sample',label = T)
DimPlot(immune.combined, group.by = 'group',label = T)
DimPlot(immune.combined, group.by = 'DoubletFinder') 
```


    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_19_0.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_19_1.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_19_2.png)
    



```R
#Findmarker
DefaultAssay(immune.combined) <- 'RNA'
Idents(immune.combined) <- 'integrated_snn_res.1.1'
immune.combined.markers <- FindAllMarkers(immune.combined,min.pct = 0.05,logfc.threshold = 0.25, only.posT=T,
                                          min.diff.pct=0.1, only.pos = T, random.seed = 666)
write.table(immune.combined.markers,file=paste(work_path,'Findmarkers_res_11_2.csv'),sep=',',quote=F)
```

    Calculating cluster 0
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Calculating cluster 8
    
    Calculating cluster 9
    
    Calculating cluster 10
    
    Calculating cluster 11
    
    Calculating cluster 12
    
    Calculating cluster 13
    
    Calculating cluster 14
    
    Calculating cluster 15
    
    Calculating cluster 16
    
    Calculating cluster 17
    
    Calculating cluster 18
    
    Calculating cluster 19
    
    Calculating cluster 20
    
    Calculating cluster 21
    
    Calculating cluster 22
    



```R
DimPlot(immune.combined, group.by = 'integrated_snn_res.1.1',label = T) +NoLegend()
DimPlot(immune.combined, group.by = 'group') 
```


    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_21_0.png)
    



    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_21_1.png)
    



```R
DefaultAssay(immune.combined) <- 'RNA'
FeaturePlot(immune.combined, features = c('Kit'))
```


    
![png](Step4_Seurat_Reduce_Dim_files/Step4_Seurat_Reduce_Dim_22_0.png)
    



```R
saveRDS(immune.combined, file = 'Step4_Seurat_Reduce_Dim.Rds')
```
