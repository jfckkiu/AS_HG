```R
############select yellow and brown module
rm(list = ls())
options(stringsAsFactors = F)
```


```R
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggpubr)
library(viridis)
library(patchwork)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library(ReactomePA)
source('~/钱教授数据分析v2.0/tools/scTools.R')
```

    ReactomePA v1.34.0  For help: https://guangchuangyu.github.io/ReactomePA
    
    If you use ReactomePA in published research, please cite:
    Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479
    



```R
package.version('enrichplot')
```


'1.10.2'



```R
package.version('DOSE')
```


'3.16.0'



```R
citation('enrichplot')
```


    
    To cite package ‘enrichplot’ in publications use:
    
      Guangchuang Yu (2021). enrichplot: Visualization of Functional
      Enrichment Result. R package version 1.10.2.
      https://yulab-smu.top/biomedical-knowledge-mining-book/
    
    A BibTeX entry for LaTeX users is
    
      @Manual{,
        title = {enrichplot: Visualization of Functional Enrichment Result},
        author = {Guangchuang Yu},
        year = {2021},
        note = {R package version 1.10.2},
        url = {https://yulab-smu.top/biomedical-knowledge-mining-book/},
      }




```R
net <- readRDS('~/AS_HG/Net.Rds')
```


```R
yellow_module <- data.table::fread('./05.module_eigengenes.xls')%>%
as_tibble() %>%
dplyr::select(samples,MEyellow) %>%
mutate(group = rep(c('H_C',"H_O",'L_C','L_O'),c(4,4,3,4))) %>%
mutate(group = factor(group, levels = c('L_C','H_C','L_O', "H_O")))
```


```R
p1 <- ggboxplot(yellow_module, title = 'Module eigengenes of yellow module', ylab = 'score', xlab = '',
          x = 'group', y = 'MEyellow',fill  = 'group', add = 'jitter', palette = 'lancet',legend = "right") +
stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('L_O', 'H_O'),c('H_C','H_O')),label = 'p.signif') +
  theme(plot.title=element_text(hjust=0.5, size = 8),
         axis.text.x = element_text(angle=0,hjust=1,size=6),
         axis.text.y = element_text(size=6),
         axis.title.y = element_text(size=6),
          legend.title =  element_text(size = 6),
           legend.text = element_text(size = 6))
p1 
```


    
![png](Step2_Enrich_Plot_files/Step2_Enrich_Plot_7_0.png)
    



```R
brown_module <- data.table::fread('./05.module_eigengenes.xls')%>%
as_tibble() %>%
dplyr::select(samples,MEbrown) %>%
mutate(group = rep(c('H_C',"H_O",'L_C','L_O'),c(4,4,3,4))) %>%
mutate(group = factor(group, levels = c('L_C','H_C','L_O', "H_O")))
```


```R
p2 <- ggboxplot(brown_module, title = 'Module eigengenes of brown module', ylab = 'score', xlab = '',
          x = 'group', y = 'MEbrown',fill  = 'group', add = 'jitter', palette = 'lancet',legend = "right") +
stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('L_O', 'H_O'),c('H_C','H_O')),label = 'p.signif') +
  theme(plot.title=element_text(hjust=0.5, size = 8),
         axis.text.x = element_text(angle=0,hjust=1,size=6),
         axis.text.y = element_text(size=6),
         axis.title.y = element_text(size=6),
          legend.title =  element_text(size = 6),
           legend.text = element_text(size = 6))
p2
```


    
![png](Step2_Enrich_Plot_files/Step2_Enrich_Plot_9_0.png)
    



```R
p0 <- p1 + p2 + plot_layout(guides = 'collect')
```


```R
ggsave(p0, file ='Module_eigengene.pdf',  units = 'cm')
```

    [1m[22mSaving 16.9 x 16.9 cm image



```R
###select top 10% percent of genes
brown_genes <- data.table::fread('~/AS_HG/08.module_result/brown-module-gene.txt') %>%
               .[1:(round(0.1 * nrow(.))),] %>%
               mutate(gene = str_to_upper(gene)) 
```


```R
####
```


```R
 brown_entrez <- as.character(na.omit(bitr(brown_genes$gene %>% str_to_upper(), #数据集
                                            fromType="SYMBOL", #输入格式
                                            toType="ENTREZID", # 转为ENTERZID格式
                                            OrgDb="org.Hs.eg.db")[,2]))
```

    Loading required package: org.Hs.eg.db
    
    
    
    'select()' returned 1:1 mapping between keys and columns
    
    Warning message in bitr(brown_genes$gene %>% str_to_upper(), fromType = "SYMBOL", :
    “8.04% of input gene IDs are fail to map...”



```R
ego_brown <- enrichDGN(brown_entrez)
```


```R
edox_brown <- setReadable(ego_brown, 'org.Hs.eg.db', 'ENTREZID')
```


```R
genes_fc_brown <- brown_genes$connectivity
names(genes_fc_brown)<- brown_genes$gene
```


```R
p1_brown <- cnetplot(edox_brown,showCategory = 3,foldChange = genes_fc_brown,categorySize="pvalue", colorEdge = TRUE)
```


```R
ggsave(p1, file ='brown_module.pdf')
```

    [1m[22mSaving 6.67 x 6.67 in image



```R
###select top 10% percent of genes
yellow_genes <- data.table::fread('~/AS_HG/08.module_result/yellow-module-gene.txt') %>%
               .[1:(round(0.1 * nrow(.))),] %>%
               mutate(gene = str_to_upper(gene)) 
```


```R
####
```


```R
 yellow_entrez <- as.character(na.omit(bitr(yellow_genes$gene %>% str_to_upper(), #数据集
                                            fromType="SYMBOL", #输入格式
                                            toType="ENTREZID", # 转为ENTERZID格式
                                            OrgDb="org.Hs.eg.db")[,2]))
```

    'select()' returned 1:1 mapping between keys and columns
    
    Warning message in bitr(yellow_genes$gene %>% str_to_upper(), fromType = "SYMBOL", :
    “9.7% of input gene IDs are fail to map...”



```R
ego_yellow <- enrichDGN(yellow_entrez)
```


```R
edox_yellow <- setReadable(ego_yellow, 'org.Hs.eg.db', 'ENTREZID')
```


```R
genes_fc_yellow <- yellow_genes$connectivity
names(genes_fc_yellow)<- yellow_genes$gene
```


```R
p2_yellow <- cnetplot(edox_yellow,showCategory = 3,foldChange = genes_fc_yellow,categorySize="pvalue", colorEdge = TRUE)
```


```R
p0_dsg <- p1_brown + p2_yellow 
ggsave(p0_dsg, file ='./FInal_Results/Figure1/Two_Module_Disease.pdf', width = 20)
```

    [1m[22mSaving 20 x 6.67 in image



    Error in grDevices::pdf(file = filename, ..., version = version): cannot open file './FInal_Results/Figure1/Two_Module_Disease.pdf'
    Traceback:


    1. ggsave(p0_dsg, file = "./FInal_Results/Figure1/Two_Module_Disease.pdf", 
     .     width = 20)

    2. dev(filename = filename, width = dim[1], height = dim[2], bg = bg, 
     .     ...)

    3. grDevices::pdf(file = filename, ..., version = version)



```R
p0_dsg
```


    
![png](Step2_Enrich_Plot_files/Step2_Enrich_Plot_28_0.png)
    

