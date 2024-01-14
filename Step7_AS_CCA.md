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
main_path <- '~/AS/AS_Mouse/'
```


```R
GSE116240_12W <- Read10X(data.dir = paste0(main_path,'GSE116240_WD12W/'))
```


```R
GSE154817_C <- Read10X(data.dir = paste0(main_path,'GSE154817/C57_ND/'))
```


```R
GSE154817_3W <- Read10X(data.dir = paste0(main_path,'GSE154817/HFD_3W/'))
```


```R
Count.list <- list(GSE116240_12W, GSE154817_C, GSE154817_3W)
```


```R
colnames(Count.list[[1]]) <- str_c(colnames(Count.list[[1]]), "_GSE116240_12W")
colnames(Count.list[[2]]) <- str_c(colnames(Count.list[[2]]), "_GSE154817_C")
colnames(Count.list[[3]]) <- str_c(colnames(Count.list[[3]]), "_GSE154817_3W")
```


```R
#create_metadata
sam.info.list <- vector("list",3) 
```


```R
sample_type <- c("GSE116240_12W","GSE154817_C","GSE154817_3W")
```


```R
group_type <- c('12W','C','3W')
```


```R
for (i in 1:3){
sam.info.list[[i]] <- as.data.frame(matrix(data=NA, nrow=ncol(Count.list[[i]]),ncol=2))
rownames(sam.info.list[[i]]) <- colnames(Count.list[[i]])
colnames(sam.info.list[[i]]) <- c("sample","group")
sam.info.list[[i]][,1] <- rep(sample_type[[i]],ncol(Count.list[[i]]))
sam.info.list[[i]][,2] <- rep(group_type[[i]], ncol(Count.list[[i]]))
}
```


```R
Seurat.list <- vector("list",3)
```


```R
for (i in 1:3) {
  Seurat.list[[i]] <- CreateSeuratObject(
  counts = Count.list[[i]], 
  meta.data=sam.info.list[[i]],
  project = "Mouse_AS_vessels", 
  min.cells = 3, min.features = 50,
  names.field = 2,names.delim = "_")
 }
```


```R
immune.combined.list <- lapply(X = Seurat.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
```


```R
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = immune.combined.list,nfeatures = 3000)
```


```R
immune.anchors <- FindIntegrationAnchors(object.list = immune.combined.list, anchor.features = features)
```

    Scaling features for provided objects
    
    Finding all pairwise anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 11406 anchors
    
    Filtering anchors
    
    	Retained 3974 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 6937 anchors
    
    Filtering anchors
    
    	Retained 3513 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 8289 anchors
    
    Filtering anchors
    
    	Retained 3985 anchors
    



```R
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
```

    Merging dataset 3 into 2
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    
    Merging dataset 1 into 2 3
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    



```R
immune.combined[['Dataset']] <- Idents(immune.combined)
```


```R
saveRDS(immune.combined, file = 'Step1_CCA_Integration.Rds')
```
