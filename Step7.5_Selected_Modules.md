```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(viridis)
source("~/钱教授数据分析v2.0/tools/scTools.R")
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
    Loading required package: viridisLite
    



```R
immune.combined <- readRDS(file = '~/AS/AS_Mouse/AS_Mouse2/Step4_Cluster_Annotatoin.Rds')
```


```R
Idents(immune.combined) <-'Annotation_Formal'
DimPlot(immune.combined, label = T) + NoLegend()
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_3_0.png)
    



```R
levels(immune.combined@meta.data$Annotation_Formal)
levels(immune.combined@meta.data$Annotation_Brief)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Resident like Macrophages'</li><li>'Inflammatory Macrophages'</li><li>'Foam cells'</li><li>'IFN induced Macrophages'</li><li>'Dendritic cells_1'</li><li>'T cells'</li><li>'Monocytes'</li><li>'Dendritic cells_2'</li><li>'Fibroblasts'</li><li>'Endothelial cells'</li><li>'Hematopoietic stem cells'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Res Mφs'</li><li>'Inflam Mφs'</li><li>'Foam cells'</li><li>'IFN Mφs'</li><li>'DCs_1'</li><li>'TCs'</li><li>'Monocytes'</li><li>'DCs_2'</li><li>'FBs'</li><li>'ECs'</li><li>'HSCs'</li></ol>




```R
library(RColorBrewer)
brewer.pal(8,"Accent")
brewer.pal(8,"Pastel2")
brewer.pal(9,"Greys")
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#7FC97F'</li><li>'#BEAED4'</li><li>'#FDC086'</li><li>'#FFFF99'</li><li>'#386CB0'</li><li>'#F0027F'</li><li>'#BF5B17'</li><li>'#666666'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#B3E2CD'</li><li>'#FDCDAC'</li><li>'#CBD5E8'</li><li>'#F4CAE4'</li><li>'#E6F5C9'</li><li>'#FFF2AE'</li><li>'#F1E2CC'</li><li>'#CCCCCC'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#FFFFFF'</li><li>'#F0F0F0'</li><li>'#D9D9D9'</li><li>'#BDBDBD'</li><li>'#969696'</li><li>'#737373'</li><li>'#525252'</li><li>'#252525'</li><li>'#000000'</li></ol>




```R
#原来的cols
colors <- c("#325A9B","#FA0087","#DAA520","#006400",'#00BFFF',
           "#EE82EE","#16FF32","#782AB6","#7B68EE","#3283FE",
           "#FF00FF")
p <- DimPlot(immune.combined, cols = colors, pt.size = 0.1,group.by = 'Annotation_Formal') + NoLegend() 
p
ggsave("Seurat_cluster.pdf", p, width = 12, height = 12)
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_6_0.png)
    



```R
FeaturePlot(immune.combined, features= 'Glud1')
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_7_0.png)
    



```R
Idents(immune.combined) = 'Annotation_Formal'
VlnPlot(immune.combined, features = 'Glud1',split.by = 'group')
```

    The default behaviour of split.by has changed.
    Separate violin plots are now plotted side-by-side.
    To restore the old behaviour of a single split violin,
    set split.plot = TRUE.
          
    This message will be shown once per session.
    
    Warning message:
    “[1m[22mGroups with fewer than two data points have been dropped.”
    Warning message:
    “[1m[22mGroups with fewer than two data points have been dropped.”



    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_8_1.png)
    



```R
data.frame(celltype =levels(Idents(immune.combined)), colors)
```


<table class="dataframe">
<caption>A data.frame: 11 × 2</caption>
<thead>
	<tr><th scope=col>celltype</th><th scope=col>colors</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Resident like Macrophages</td><td>#325A9B</td></tr>
	<tr><td>Inflammatory Macrophages </td><td>#FA0087</td></tr>
	<tr><td>Foam cells               </td><td>#DAA520</td></tr>
	<tr><td>IFN induced Macrophages  </td><td>#006400</td></tr>
	<tr><td>Dendritic cells_1        </td><td>#00BFFF</td></tr>
	<tr><td>T cells                  </td><td>#EE82EE</td></tr>
	<tr><td>Monocytes                </td><td>#16FF32</td></tr>
	<tr><td>Dendritic cells_2        </td><td>#782AB6</td></tr>
	<tr><td>Fibroblasts              </td><td>#7B68EE</td></tr>
	<tr><td>Endothelial cells        </td><td>#3283FE</td></tr>
	<tr><td>Hematopoietic stem cells </td><td>#FF00FF</td></tr>
</tbody>
</table>




```R
p <- DimPlot(immune.combined,group.by = 'group')  + theme_void()+ NoLegend()  + 
scale_color_manual(values = c('C' = '#0C5F8B','3W' ='#B99831','12W' ='#BD302A'))
p
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_10_0.png)
    



```R
#原来的cols
p <- DimPlot(immune.combined,group.by = 'group')  + theme_void()+ NoLegend()  + 
scale_color_manual(values = c('C' = '#0C5F8B','3W' ='#B99831','12W' ='#BD302A'))
p
ggsave("Seurat_group.pdf", p, width = 12, height = 12)
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_11_0.png)
    



```R
modules_ME_yellow  <- data.table::fread('~/AS_HG/08.module_result/yellow-module-gene.txt')
modules_ME_brown  <- data.table::fread('~/AS_HG/08.module_result/brown-module-gene.txt')
```


```R
yellow_top100  <- modules_ME_yellow  %>% 
  arrange(-connectivity)  %>% 
   .[1:100,] %>% 
    pull(gene) 
brown_top100  <- modules_ME_brown  %>% 
  arrange(-connectivity)  %>% 
   .[1:100,] %>% 
   pull(gene)
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
    “The following features are not present in the object: Mus_musculus_newGene_530, Mus_musculus_newGene_707, Mus_musculus_newGene_1249, not searching for symbol synonyms”



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
    “The following features are not present in the object: Otulinl, Mus_musculus_newGene_1705, Mus_musculus_newGene_3963, Mus_musculus_newGene_3605, not searching for symbol synonyms”



```R
library(scCustomize)
```

    [1m[22mscCustomize v1.1.3
    If you find the scCustomize useful please cite.
    See 'samuel-marsh.github.io/scCustomize/articles/FAQ.html' for citation info.
    



```R
FeaturePlot(immune.combined, features = 'Yellow_Module1')
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_17_0.png)
    



```R
FeaturePlot(immune.combined, features = 'Brown_Module1')
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_18_0.png)
    



```R
DefaultAssay(immune.combined)  <- 'RNA'
```


```R
pal <- viridis(n = 10, option = "D")
```


```R
FeaturePlot_scCustom(immune.combined,features = 'Yellow_Module1',cols = viridis_magma_dark_high)
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



    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_21_1.png)
    



```R
VlnPlot(immune.combined, features = 'Yellow_Module1')
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_22_0.png)
    



```R
VlnPlot(immune.combined, features = 'Brown_Module1')
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_23_0.png)
    



```R
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
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
#原来的cols
colors <- c("#325A9B","#FA0087","#DAA520","#006400",'#00BFFF',
           "#EE82EE","#16FF32","#782AB6","#7B68EE","#3283FE",
           "#FF00FF")
```


```R
Idents(immune.combined)  <- 'Annotation_Brief'
```


```R
StackedVlnPlot(immune.combined, c('Yellow_Module1', 'Brown_Module1'), pt.size=0,cols=colors) + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 45)
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_27_0.png)
    



```R
p <- StackedVlnPlot(immune.combined, c('Yellow_Module1', 'Brown_Module1'), pt.size=0,cols=colors) + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 45)
ggsave(p, file = './Two_Module_VlnPlot.pdf')
```

    [1m[22mSaving 6.67 x 6.67 in image
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Res Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'Inflam Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <cf>”
    Warning message in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    “conversion failure on 'IFN Mφs' in 'mbcsToSbcs': dot substituted for <86>”



```R
getwd()
```


'/home/jfckkiu/AS_HG/Final_Results/Figure4_New_Mouse_Cd45'



```R
DefaultAssay(immune.combined)  <- 'RNA'
p1  <- FeaturePlot_scCustom(immune.combined, features = 'Aldoa') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
p2  <- FeaturePlot_scCustom(immune.combined, features = 'Creg1') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
p3  <- FeaturePlot_scCustom(immune.combined, features = 'Lgmn') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
p4  <- FeaturePlot_scCustom(immune.combined, features = 'Pkm') + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 8))
```


```R
p0 <- (p1|p2)/(p3|p4) 
p0
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_31_0.png)
    



```R
ggsave(p0, file = './Selected_Gene_FeaturePlot.pdf')
```

    [1m[22mSaving 6.67 x 6.67 in image



```R
DefaultAssay(immune.combined)  <- 'RNA'
p1  <- FeaturePlot_scCustom(immune.combined, features = 'Aldoa',pt.size = 0.1) + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none')
p2  <- FeaturePlot_scCustom(immune.combined, features = 'Creg1',pt.size = 0.1) + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none')
p3  <- FeaturePlot_scCustom(immune.combined, features = 'Lgmn',pt.size = 0.1) + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none')
p4  <- FeaturePlot_scCustom(immune.combined, features = 'Pkm',pt.size = 0.1) + theme_void() +theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5, size = 0),legend.position = 'none') 
```


```R
p0.1 <- (p1|p2|p3|p4) 
p0.1
```


    
![png](Step7.5_Selected_Modules_files/Step7.5_Selected_Modules_34_0.png)
    



```R
ggsave(p0.1, file = './Selected_Gene_FeaturePlot.tiff', height = 5, width = 19, units = 'cm')
```

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

