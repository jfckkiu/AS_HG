```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(doParallel)
library(RColorBrewer)
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
    



```R
immune.combined <- readRDS(file = '~/AS/AS_Mouse/AS_Mouse2/Step4_Cluster_Annotatoin.Rds')
```


```R
#Findmarker
DefaultAssay(immune.combined) <- 'RNA'
Idents(immune.combined) <- 'Annotation_Formal'
immune.combined.markers <- FindAllMarkers(immune.combined,min.pct = 0.05,logfc.threshold = 0.25, only.posT=T,
                                          min.diff.pct=0.1, only.pos = T, random.seed = 666)
write.table(immune.combined.markers,file='~/AS/AS_Mouse/AS_Mouse2/Findmarkers_celltype.csv',sep=',',quote=F)
```

    Calculating cluster Resident like Macrophages
    
    Calculating cluster Inflammatory Macrophages
    
    Calculating cluster Foam cells
    
    Calculating cluster IFN induced Macrophages
    
    Calculating cluster Dendritic cells_1
    
    Calculating cluster T cells
    
    Calculating cluster Monocytes
    
    Calculating cluster Dendritic cells_2
    
    Calculating cluster Fibroblasts
    
    Calculating cluster Endothelial cells
    
    Calculating cluster Hematopoietic stem cells
    

