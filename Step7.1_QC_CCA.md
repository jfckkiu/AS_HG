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
main_path <- '~/AS/AS_Mouse/AS_Mouse2/'
```


```R
immune.combined <- readRDS(paste0(main_path,'Step1_CCA_Integration.Rds'))
```


```R
DefaultAssay(immune.combined) <- 'RNA'
```


```R
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern = "^mt-")
```


```R
Idents(immune.combined) <- 'sample'
```


```R
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


    
![png](Step7.1_QC_CCA_files/Step7.1_QC_CCA_7_0.png)
    



```R
quantile(immune.combined@meta.data$nCount_RNA,  c(.25, .50,  .75, .90, .95,.99))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>25%</dt><dd>5490</dd><dt>50%</dt><dd>7186</dd><dt>75%</dt><dd>9662</dd><dt>90%</dt><dd>13538.8</dd><dt>95%</dt><dd>16705</dd><dt>99%</dt><dd>24466.08</dd></dl>




```R
quantile(immune.combined@meta.data$nFeature_RNA,  c(.25, .50,  .75, .90, .95,.99))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>25%</dt><dd>1857</dd><dt>50%</dt><dd>2232</dd><dt>75%</dt><dd>2715</dd><dt>90%</dt><dd>3312</dd><dt>95%</dt><dd>3781</dd><dt>99%</dt><dd>4730.96</dd></dl>




```R
plot1 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```


    
![png](Step7.1_QC_CCA_files/Step7.1_QC_CCA_10_0.png)
    



```R
immune.combined <- subset(immune.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA <20000
               & percent.mt < 25 )
```


```R
plot1 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```


    
![png](Step7.1_QC_CCA_files/Step7.1_QC_CCA_12_0.png)
    



```R
saveRDS(immune.combined, file = 'Step2_QC_CCA_Seurat.Rds')
```
