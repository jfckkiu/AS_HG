```R

```


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
library(patchwork)
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
    
    Loading required package: ggraph
    



```R
plan("multicore", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2) #50G
```


```R
immune.combined <- readRDS('~/AS/AS_Mouse/AS_Mouse2/Final_result/Figure8_Human_Plaques/Step1_AddDoublet_Seurat.Rds')
```


```R
DimPlot(immune.combined,group.by ='celltype_Brief',label=T)+NoLegend()
```


    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_5_0.png)
    



```R
modules_ME_yellow  <- data.table::fread('~/AS_HG/08.module_result/yellow-module-gene.txt')
modules_ME_brown  <- data.table::fread('~/AS_HG/08.module_result/brown-module-gene.txt')
```


```R
yellow_top100  <- modules_ME_yellow  %>% 
  arrange(-connectivity)  %>% 
   .[1:100,] %>% 
    pull(gene)  %>% 
     str_to_upper()
brown_top100  <- modules_ME_brown  %>% 
  arrange(-connectivity)  %>% 
   .[1:100,] %>% 
   pull(gene) %>% 
     str_to_upper()
```


```R
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(yellow_top100),
                 ctrl = length(yellow_top100),
                 name = 'Yellow_Module',
                 search = FALSE,
                 assay = "RNA")
```

    Warning message:
    “The following features are not present in the object: TNFRSF26, C77080, ZFP628, GM10521, MUS_MUSCULUS_NEWGENE_530, ZFP503, MUS_MUSCULUS_NEWGENE_707, ZFP335, MUS_MUSCULUS_NEWGENE_1249, 4833420G17RIK, not searching for symbol synonyms”



```R
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(brown_top100),
                 ctrl = length(brown_top100),
                 name = 'Brown_Module',
                 search = FALSE,
                 assay = "RNA")
```

    Warning message:
    “The following features are not present in the object: MAN1A, CLEC4N, MUS_MUSCULUS_NEWGENE_1705, MUS_MUSCULUS_NEWGENE_3963, MUS_MUSCULUS_NEWGENE_3605, AHSA2, G6PDX, not searching for symbol synonyms”



```R
library(scCustomize)
```

    [1m[22mscCustomize v1.1.3
    If you find the scCustomize useful please cite.
    See 'samuel-marsh.github.io/scCustomize/articles/FAQ.html' for citation info.
    



```R
DefaultAssay(immune.combined)  <- 'RNA'
```


```R
pal <- viridis::viridis(n = 10, option = "D")
```


```R
FeaturePlot_scCustom(immune.combined,features = 'Yellow_Module1')
```

    Warning message:
    “[1m[22mSome of the plotted features are from meta.data slot.
    [36m•[39m Please check that `na_cutoff` param is being set appropriately for those features.”
    [1m[22m
    NOTE: [32mFeaturePlot_scCustom[39m uses a specified `na_cutoff` when plotting to
    color cells with no expression as background color separate from color scale.
    Please ensure `na_cutoff` value is appropriate for feature being plotted.
    Default setting is appropriate for use when plotting from 'RNA' assay.
    When `na_cutoff` not appropriate (e.g., module scores) set to NULL to
    plot all cells in gradient color palette.
    
    -----This message will be shown once per session.-----



    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_13_1.png)
    



```R
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               #axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = plot.margin )+ stat_summary(fun = mean, geom='point', size = 10, colour = "black", shape = 45)
       return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(angle = 90,hjust = 1), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)
}

```


```R
levels(immune.combined@meta.data$celltype_Brief)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Mono/Mø'</li><li>'SMC'</li><li>'EC'</li><li>'T'</li><li>'NK'</li><li>'B'</li><li>'PC'</li><li>'MC'</li><li>'PLF'</li></ol>




```R
immune.combined@meta.data$celltype_Brief  <- factor(immune.combined@meta.data$celltype_Brief, 
                                                    levels = c('Mono/Mø','SMC','EC','T','NK','B','PC','MC','PLF'))
```


```R
#原来的cols
color <- c('#F6313E','#4FA8DA','#FFA300','#46A040','#89774A','#A65AC2','#0DB2AA','#00441B','#9A9A9A')
```


```R
Idents(immune.combined)  <- 'celltype_Brief'
```


```R
StackedVlnPlot(immune.combined, c('Yellow_Module1', 'Brown_Module1'), pt.size=0,cols=color) + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 45)
```


    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_19_0.png)
    



```R
p <- StackedVlnPlot(immune.combined, c('Yellow_Module1', 'Brown_Module1'), pt.size=0,cols=color) + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 45)
ggsave(p, file = './Two_Module_VlnPlot.pdf')
```

    [1m[22mSaving 6.67 x 6.67 in image



```R
DefaultAssay(immune.combined)  <- 'RNA'
FeaturePlot_scCustom(immune.combined, features = 'ALDOA') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
FeaturePlot_scCustom(immune.combined, features = 'CREG1') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
FeaturePlot_scCustom(immune.combined, features = 'LGMN') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
FeaturePlot_scCustom(immune.combined, features = 'PKM') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
```


    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_21_0.png)
    



    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_21_1.png)
    



    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_21_2.png)
    



    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_21_3.png)
    



```R
DefaultAssay(immune.combined)  <- 'RNA'
p1  <- FeaturePlot_scCustom(immune.combined, features = 'ALDOA') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
p2  <- FeaturePlot_scCustom(immune.combined, features = 'CREG1') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
p3  <- FeaturePlot_scCustom(immune.combined, features = 'LGMN') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
p4  <- FeaturePlot_scCustom(immune.combined, features = 'PKM') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
```


```R
p0 <- (p1|p2)/(p3|p4) 
p0
```


    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_23_0.png)
    



```R
ggsave(p0, file = './Selected_Gene_FeaturePlot.pdf')
```

    [1m[22mSaving 6.67 x 6.67 in image



```R
DefaultAssay(immune.combined)  <- 'RNA'
p1  <- FeaturePlot_scCustom(immune.combined, features = 'ALDOA') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none')
p2  <- FeaturePlot_scCustom(immune.combined, features = 'CREG1') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none')
p3  <- FeaturePlot_scCustom(immune.combined, features = 'LGMN') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none')
p4  <- FeaturePlot_scCustom(immune.combined, features = 'PKM') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none')
```


```R
p0.1 <- (p1|p2)/(p3|p4) 
p0.1
```


    
![png](Step10.5_Human_Plaque_new_files/Step10.5_Human_Plaque_new_26_0.png)
    



```R
ggsave(p0.1, file = './Selected_Gene_FeaturePlot.tiff')
```

    [1m[22mSaving 6.67 x 6.67 in image
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “Unable to calculate text width/height (using zero)”

