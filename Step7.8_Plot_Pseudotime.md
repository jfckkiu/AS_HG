```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(tradeSeq)
library(slingshot)
library(patchwork)
```


```R
immune.combined <- readRDS('~/AS/AS_Mouse/AS_Mouse2/Final_result/Figure5_Monocle_Macro/Step4_Macro_pseudotime.Rds')
```


```R
DefaultAssay(immune.combined) <- 'RNA'
immune.combined <- NormalizeData(immune.combined)
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 1000)
genes <- VariableFeatures(immune.combined)
```


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


```R
DimPlot(immune.combined)
```


    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_6_0.png)
    



```R
counts_cds <- immune.combined@assays$RNA@counts
```


```R
set.seed(3)
pdf("RNA_evaluateK.pdf")
icMat <- evaluateK(counts = counts_cds, sds = ss, k=3:10, nGenes = 200,parallel = TRUE,
                   verbose=TRUE)
dev.off()
```

      |======================================================================| 100%
    
      |======================================================================| 100%
    
      |======================================================================| 100%
    
      |======================================================================| 100%
    
      |======================================================================| 100%
    
      |======================================================================| 100%
    
      |======================================================================| 100%
    
      |======================================================================| 100%
    



<strong>png:</strong> 2



```R
group <- immune.combined@meta.data$group
group <- as.character(group)
group[group == c("C")] =1
group[group == c("3W")] =2
group[group == c("12W")] =2
group <- as.factor(group)
```


```R
sce   <- fitGAM(counts = counts_cds[genes,],sds = ss,
                  conditions = group,nknots = 7,parallel = TRUE)
```

    Fitting lineages with multiple conditions. This method has been tested on a couple of datasets, but is still in an experimental phase.
    


      |======================================================================| 100%
    

saveRDS(sce,file = 'RNA_tradeseq.Rds')

```R
startRes <- startVsEndTest(sce)
```


```R
startRes %>%
 arrange(-waldStat) %>%
filter(logFClineage1 > 0) %>%
rownames_to_column(var = 'gene') %>%
mutate(rank = 1:nrow(.)) 
```


<table class="dataframe">
<caption>A data.frame: 438 × 6</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>waldStat</th><th scope=col>df</th><th scope=col>pvalue</th><th scope=col>logFClineage1</th><th scope=col>rank</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Cd63   </td><td>105.91596</td><td>1</td><td>0.000000e+00</td><td>4.933959</td><td> 1</td></tr>
	<tr><td>Ctsb   </td><td>105.65869</td><td>1</td><td>0.000000e+00</td><td>3.341890</td><td> 2</td></tr>
	<tr><td>Ftl1   </td><td> 83.89654</td><td>1</td><td>0.000000e+00</td><td>2.187480</td><td> 3</td></tr>
	<tr><td>Cd68   </td><td> 81.52365</td><td>1</td><td>0.000000e+00</td><td>2.968081</td><td> 4</td></tr>
	<tr><td>Fth1   </td><td> 78.61124</td><td>1</td><td>0.000000e+00</td><td>2.279504</td><td> 5</td></tr>
	<tr><td>Lgmn   </td><td> 78.14355</td><td>1</td><td>0.000000e+00</td><td>3.145638</td><td> 6</td></tr>
	<tr><td>Ms4a7  </td><td> 66.78184</td><td>1</td><td>3.330669e-16</td><td>4.737644</td><td> 7</td></tr>
	<tr><td>Ctsd   </td><td> 60.43556</td><td>1</td><td>7.549517e-15</td><td>3.963134</td><td> 8</td></tr>
	<tr><td>Cstb   </td><td> 48.00360</td><td>1</td><td>4.254375e-12</td><td>3.578606</td><td> 9</td></tr>
	<tr><td>Ctss   </td><td> 46.61052</td><td>1</td><td>8.659407e-12</td><td>2.186612</td><td>10</td></tr>
	<tr><td>Ctsz   </td><td> 46.24197</td><td>1</td><td>1.045142e-11</td><td>2.243448</td><td>11</td></tr>
	<tr><td>Prdx1  </td><td> 45.84483</td><td>1</td><td>1.279998e-11</td><td>2.590644</td><td>12</td></tr>
	<tr><td>Trem2  </td><td> 45.58800</td><td>1</td><td>1.459333e-11</td><td>4.316834</td><td>13</td></tr>
	<tr><td>Lgals3 </td><td> 41.12575</td><td>1</td><td>1.427430e-10</td><td>3.598529</td><td>14</td></tr>
	<tr><td>Cd72   </td><td> 35.59832</td><td>1</td><td>2.424962e-09</td><td>5.513323</td><td>15</td></tr>
	<tr><td>Pld3   </td><td> 34.32934</td><td>1</td><td>4.653126e-09</td><td>5.498683</td><td>16</td></tr>
	<tr><td>Lpl    </td><td> 33.37772</td><td>1</td><td>7.588839e-09</td><td>5.559441</td><td>17</td></tr>
	<tr><td>Vcam1  </td><td> 32.12882</td><td>1</td><td>1.442808e-08</td><td>7.788103</td><td>18</td></tr>
	<tr><td>Lipa   </td><td> 31.09428</td><td>1</td><td>2.457952e-08</td><td>3.912834</td><td>19</td></tr>
	<tr><td>Gclm   </td><td> 30.37973</td><td>1</td><td>3.552202e-08</td><td>2.701063</td><td>20</td></tr>
	<tr><td>Hexb   </td><td> 29.73215</td><td>1</td><td>4.960517e-08</td><td>3.058178</td><td>21</td></tr>
	<tr><td>Vim    </td><td> 27.35488</td><td>1</td><td>1.693365e-07</td><td>3.016714</td><td>22</td></tr>
	<tr><td>Capg   </td><td> 26.91395</td><td>1</td><td>2.127181e-07</td><td>3.282377</td><td>23</td></tr>
	<tr><td>Adcy3  </td><td> 22.44847</td><td>1</td><td>2.158576e-06</td><td>5.651157</td><td>24</td></tr>
	<tr><td>Cd9    </td><td> 22.42987</td><td>1</td><td>2.179575e-06</td><td>3.307456</td><td>25</td></tr>
	<tr><td>Id2    </td><td> 22.19216</td><td>1</td><td>2.466795e-06</td><td>3.491650</td><td>26</td></tr>
	<tr><td>Psap   </td><td> 21.85862</td><td>1</td><td>2.934955e-06</td><td>1.343514</td><td>27</td></tr>
	<tr><td>Gpr137b</td><td> 21.42901</td><td>1</td><td>3.671736e-06</td><td>4.054308</td><td>28</td></tr>
	<tr><td>Mmp12  </td><td> 21.03829</td><td>1</td><td>4.501963e-06</td><td>7.632419</td><td>29</td></tr>
	<tr><td>Sdf2l1 </td><td> 20.95009</td><td>1</td><td>4.714047e-06</td><td>2.178872</td><td>30</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>Lipg     </td><td>8.560094e-05</td><td>1</td><td>0.9926180</td><td>2.062402e-01</td><td>409</td></tr>
	<tr><td>Phlda3   </td><td>5.392512e-05</td><td>1</td><td>0.9941409</td><td>1.956846e-02</td><td>410</td></tr>
	<tr><td>Scgb1a1  </td><td>1.738339e-05</td><td>1</td><td>0.9966734</td><td>5.846166e+04</td><td>411</td></tr>
	<tr><td>Gm14461  </td><td>2.999702e-09</td><td>1</td><td>0.9999563</td><td>7.309868e+02</td><td>412</td></tr>
	<tr><td>Gzmb     </td><td>2.236464e-11</td><td>1</td><td>0.9999962</td><td>6.311774e+01</td><td>413</td></tr>
	<tr><td>Mmp27    </td><td>2.396970e-12</td><td>1</td><td>0.9999988</td><td>2.066340e+01</td><td>414</td></tr>
	<tr><td>Tmem178  </td><td>1.966053e-12</td><td>1</td><td>0.9999989</td><td>1.871407e+01</td><td>415</td></tr>
	<tr><td>Adgrg3   </td><td>1.645827e-12</td><td>1</td><td>0.9999990</td><td>1.712232e+01</td><td>416</td></tr>
	<tr><td>Tspan15  </td><td>1.209521e-12</td><td>1</td><td>0.9999991</td><td>7.857944e+00</td><td>417</td></tr>
	<tr><td>Klra3    </td><td>1.180062e-12</td><td>1</td><td>0.9999991</td><td>1.449850e+01</td><td>418</td></tr>
	<tr><td>Serpina3g</td><td>1.043277e-12</td><td>1</td><td>0.9999992</td><td>1.363234e+01</td><td>419</td></tr>
	<tr><td>Tox      </td><td>9.806597e-13</td><td>1</td><td>0.9999992</td><td>1.321690e+01</td><td>420</td></tr>
	<tr><td>Dmkn     </td><td>8.267122e-13</td><td>1</td><td>0.9999993</td><td>1.213522e+01</td><td>421</td></tr>
	<tr><td>Bdh2     </td><td>6.257864e-13</td><td>1</td><td>0.9999994</td><td>1.055805e+01</td><td>422</td></tr>
	<tr><td>Npy      </td><td>3.538437e-13</td><td>1</td><td>0.9999995</td><td>7.939191e+00</td><td>423</td></tr>
	<tr><td>Nkg7     </td><td>3.456652e-13</td><td>1</td><td>0.9999995</td><td>7.846904e+00</td><td>424</td></tr>
	<tr><td>Htr2b    </td><td>3.400116e-13</td><td>1</td><td>0.9999995</td><td>7.782469e+00</td><td>425</td></tr>
	<tr><td>Tenm4    </td><td>3.240470e-13</td><td>1</td><td>0.9999995</td><td>7.597566e+00</td><td>426</td></tr>
	<tr><td>Shisa9   </td><td>2.414452e-13</td><td>1</td><td>0.9999996</td><td>6.558127e+00</td><td>427</td></tr>
	<tr><td>Scn3b    </td><td>2.265837e-13</td><td>1</td><td>0.9999996</td><td>6.353087e+00</td><td>428</td></tr>
	<tr><td>Il12a    </td><td>2.061214e-13</td><td>1</td><td>0.9999996</td><td>6.059433e+00</td><td>429</td></tr>
	<tr><td>Rhov     </td><td>1.965373e-13</td><td>1</td><td>0.9999996</td><td>5.916883e+00</td><td>430</td></tr>
	<tr><td>Tarm1    </td><td>1.273269e-13</td><td>1</td><td>0.9999997</td><td>4.762452e+00</td><td>431</td></tr>
	<tr><td>Kcnj13   </td><td>1.195567e-13</td><td>1</td><td>0.9999997</td><td>4.614849e+00</td><td>432</td></tr>
	<tr><td>Gm44205  </td><td>1.080024e-13</td><td>1</td><td>0.9999997</td><td>4.386189e+00</td><td>433</td></tr>
	<tr><td>Ly6d     </td><td>9.066867e-14</td><td>1</td><td>0.9999998</td><td>4.018825e+00</td><td>434</td></tr>
	<tr><td>Arg2     </td><td>6.494564e-14</td><td>1</td><td>0.9999998</td><td>3.401304e+00</td><td>435</td></tr>
	<tr><td>Cyp1b1   </td><td>3.442253e-14</td><td>1</td><td>0.9999999</td><td>2.476235e+00</td><td>436</td></tr>
	<tr><td>Tacstd2  </td><td>9.513914e-16</td><td>1</td><td>1.0000000</td><td>4.116708e-01</td><td>437</td></tr>
	<tr><td>Arg1     </td><td>4.739125e-16</td><td>1</td><td>1.0000000</td><td>2.905491e-01</td><td>438</td></tr>
</tbody>
</table>




```R
selected_genes_up <- startRes %>%
 arrange(-waldStat) %>%
filter(logFClineage1 > 0) %>%
rownames_to_column(var = 'gene') %>%
mutate(rank = 1:nrow(.)) %>%
pull(gene) %>%
.[c(1:6,11,22)]
```


```R
selected_genes_down <- startRes %>%
 arrange(-waldStat) %>%
filter(logFClineage1 < 0) %>%
rownames_to_column(var = 'gene') %>%
mutate(rank = 1:nrow(.)) %>%
pull(gene) %>%
.[c(1:4,6:7,18,19)]
```


```R
immune.combined@meta.data$Annotation_Brief %>% unique()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>Inflam Mφs</li><li>Foam cells</li><li>Monocytes</li><li>IFN Mφs</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Res Mφs'</li><li>'Inflam Mφs'</li><li>'Foam cells'</li><li>'IFN Mφs'</li><li>'DCs_1'</li><li>'TCs'</li><li>'Monocytes'</li><li>'DCs_2'</li><li>'FBs'</li><li>'ECs'</li><li>'HSCs'</li></ol>
</details>



```R
p0 <- tibble(celltype = immune.combined@meta.data$Annotation_Brief, pseudotime = immune.combined@meta.data$curve1_pct)  %>% 
 mutate(celltype = fct_recode(celltype, 'Inflam Mφs' ='IFN Mφs')) %>% 
 mutate(cell = 'Macrophage_Lineage') %>% 
  ggplot(aes(cell,pseudotime)) + geom_jitter(aes(colour = celltype),width = 0.2,size = 0.1) + theme_void() + theme(legend.position = 'none') + 
   scale_color_manual(values = c('Monocytes'= "#16FF32", 'Inflam Mφs' = "#FA0087", 'Foam cells' = "#DAA520")) + coord_flip()
p0
```


    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_17_0.png)
    



```R
p1  <- tibble(expression = immune.combined@assays$RNA@data['Aldoa',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime', title = 'Aldoa') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
p1
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_18_1.png)
    



```R
p0 / p1 
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_19_1.png)
    



```R
p2  <- tibble(expression = immune.combined@assays$RNA@data['Creg1',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime', title = 'Creg1') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
p2
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_20_1.png)
    



```R
p3  <- tibble(expression = immune.combined@assays$RNA@data['Lgmn',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime', title = 'Lgmn') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
p3
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_21_1.png)
    



```R
p4  <- tibble(expression = immune.combined@assays$RNA@data['Pkm',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime', title = 'Pkm') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
p4
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_22_1.png)
    



```R
p <- (p0 / p1 ) | (p0 / p2 ) | (p0 / p3 ) | (p0 / p4 ) 
ggsave(p, file = './Selected_Gene_Pseudotime.pdf', height = 8, width = 18, units = 'cm')
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



```R
tibble(expression = immune.combined@assays$RNA@data['Lpl',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
```

    Warning message:
    “[1m[22mUsing `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    [36mℹ[39m Please use `linewidth` instead.”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_24_1.png)
    



```R
tibble(expression = immune.combined@assays$RNA@data['Aldoa',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_25_1.png)
    



```R
tibble(expression = immune.combined@assays$RNA@data['Creg1',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_26_1.png)
    



```R
tibble(expression = immune.combined@assays$RNA@data['Lgmn',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_27_1.png)
    



```R
tibble(expression = immune.combined@assays$RNA@data['Pkm',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step7.8_Plot_Pseudotime_files/Step7.8_Plot_Pseudotime_28_1.png)
    



```R
tibble(expression = immune.combined@assays$RNA@data['Aldoa',], pseudotime = immune.combined@meta.data$curve1_pct) %>%
ggplot(aes(x = pseudotime, y = expression)) + geom_smooth(se = TRUE, method = 'loess',size =0.2) + 
      theme_classic() + labs( y = 'Gene expression', x = 'Scaled_Pseudotime') +
      theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(hjust=1),axis.title.x = element_text(), axis.title.y = element_text(), 
         axis.text.y = element_text(),legend.title = element_text(),
           legend.text = element_text())
```


```R
immune.combined@assays$RNA@data[union(selected_genes_up,selected_genes_down),] %>%
as.data.frame() %>%
t() %>%
bind_cols(pseudotime = immune.combined@meta.data$curve1_pct) %>%
pivot_longer(-pseudotime, names_to = 'gene', values_to = 'counts') %>%
mutate(gene = factor(gene,levels = c(selected_genes_up,selected_genes_down))) %>%
ggplot()  + geom_smooth(mapping = aes(x = pseudotime, y = counts),se = TRUE, method = 'loess',size =0.2) +
ggpubr::theme_classic2() + facet_wrap(~gene,scales = "free") + 
    labs(y = 'Gene expression', x = 'Scaled_Pseudotime') +
      theme(plot.title=element_text(hjust=0.5,size = 4),axis.text.x = element_text(hjust=1,size = 0),axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), 
         axis.text.y = element_text(size = 0),legend.title = element_text(),legend.position = 'none',
           legend.text = element_text(),strip.background = element_blank(),strip.text = element_text(size=6)) + ggthemes::scale_color_hc()
```
