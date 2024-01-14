```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(MySeuratWrappers)  
source("~/钱教授数据分析v2.0/tools/scTools.R")
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
    
    Loading required package: reticulate
    
    
    Attaching package: ‘MySeuratWrappers’
    
    
    The following objects are masked from ‘package:Seurat’:
    
        DimPlot, DoHeatmap, LabelClusters, RidgePlot, VlnPlot
    
    



```R
immune.combined <- readRDS(file = '~/AS/AS_Mouse/AS_Mouse2/Step4_Cluster_Annotatoin.Rds')
```


```R
levels(immune.combined@meta.data$Annotation_Brief)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Res Mφs'</li><li>'Inflam Mφs'</li><li>'Foam cells'</li><li>'IFN Mφs'</li><li>'DCs_1'</li><li>'TCs'</li><li>'Monocytes'</li><li>'DCs_2'</li><li>'FBs'</li><li>'ECs'</li><li>'HSCs'</li></ol>




```R
immune.combined@meta.data$Annotation_Brief <- factor(immune.combined@meta.data$Annotation_Brief,
                                        levels =  c("Res Mφs",'Inflam Mφs','IFN Mφs','Foam cells','Monocytes',
                                            'DCs_1','DCs_2','TCs','FBs','ECs',
                                             "HSCs"))
```


```R
Idents(immune.combined) <- 'Annotation_Brief'
```


```R
markers = c('Ptprc','Cd14','Fcgr3','Fcgr1','Cd68','Lyve1','Nlrp3',
            'Lgals3','Cx3cr1','Ly6c2','Itgax','Cd209a','Xcr1',
           'Cd3d','Cd4','Cd8a','Acta2','Pecam1','Kit')
```


```R
#原来的cols
colors <- c("#325A9B","#FA0087","#DAA520","#006400",'#00BFFF',
           "#EE82EE","#16FF32","#782AB6","#7B68EE","#3283FE",
           "#FF00FF")
```


```R
VlnPlot(immune.combined, features = markers,  
        stacked=T,pt.size=0,  
        cols = colors,#颜色  
       direction = "horizontal", #水平作图  
        x.lab = '', y.lab = '')+#横纵轴不标记任何东西  
  theme(axis.text.x = element_blank(),   
      axis.ticks.x = element_blank(),
       axis.text.y = element_text(size = 6, hjust = 1),
       plot.title=element_text(hjust=0.5, size = 6)) 
```


    
![png](Step7.4_Vln_Plot_files/Step7.4_Vln_Plot_8_0.png)
    



```R
p <- VlnPlot(immune.combined, features = markers,  
        stacked=T,pt.size=0,  
        cols = colors,#颜色  
       direction = "horizontal", #水平作图  
        x.lab = '', y.lab = '')+#横纵轴不标记任何东西  
  theme(axis.text.x = element_blank(),   
      axis.ticks.x = element_blank(),
       axis.text.y = element_text(size = 6, hjust = 1),
       plot.title=element_text(hjust=0.5, size = 6)) 
ggsave(p, file = 'Selected_gene_Vnplot.pdf', height = 10, width = 10, units = 'cm')
```

    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”

