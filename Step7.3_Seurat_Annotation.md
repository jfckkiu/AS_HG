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
work_path <- paste0(getwd(),'/')
```


```R
immune.combined <- readRDS(file = paste0(work_path, 'Step4_Seurat_Reduce_Dim.Rds'))
```


```R
Idents(immune.combined) <- 'integrated_snn_res.1.1'
```


```R
#Findmarker
#DefaultAssay(immune.combined) <- 'RNA'
#Idents(immune.combined) <- 'integrated_snn_res.1.1'
#immune.combined.markers <- FindAllMarkers(immune.combined,min.pct = 0.05,logfc.threshold = 0.25, only.posT=T,
#                                          min.diff.pct=0.1, only.pos = T, random.seed = 666)
#write.table(immune.combined.markers,file=paste(work_path,'Findmarkers_res_11_2.csv'),sep=',',quote=F)
```


```R
DefaultAssay(immune.combined) <- 'RNA'
DimPlot(immune.combined,group.by = 'integrated_snn_res.1.1',label = T) + NoLegend()
DimPlot(immune.combined,group.by = 'group') 
```


    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_6_0.png)
    



    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_6_1.png)
    



```R
FeaturePlot(immune.combined, features = c('Cd9'))
```


    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_7_0.png)
    



```R
markers <- FindMarkers(immune.combined,ident.1 =15, ident.2 = 16,min.pct = 0.05,logfc.threshold = 0.25)
markers %>%
arrange(avg_log2FC)
```


<table class="dataframe">
<caption>A data.frame: 1545 × 5</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>S100a9</th><td>1.864518e-03</td><td>-4.756436</td><td>0.041</td><td>0.131</td><td>1.000000e+00</td></tr>
	<tr><th scope=row>S100a8</th><td>3.728914e-02</td><td>-4.730024</td><td>0.088</td><td>0.153</td><td>1.000000e+00</td></tr>
	<tr><th scope=row>Plac8</th><td>7.070759e-43</td><td>-4.635609</td><td>0.041</td><td>0.759</td><td>1.172403e-38</td></tr>
	<tr><th scope=row>Il1b</th><td>1.790189e-17</td><td>-3.204045</td><td>0.268</td><td>0.664</td><td>2.968313e-13</td></tr>
	<tr><th scope=row>Slpi</th><td>4.823473e-08</td><td>-2.909235</td><td>0.036</td><td>0.226</td><td>7.997800e-04</td></tr>
	<tr><th scope=row>Chil3</th><td>2.666928e-18</td><td>-2.763455</td><td>0.010</td><td>0.365</td><td>4.422034e-14</td></tr>
	<tr><th scope=row>Msrb1</th><td>6.761738e-43</td><td>-2.660776</td><td>0.366</td><td>0.949</td><td>1.121164e-38</td></tr>
	<tr><th scope=row>Gngt2</th><td>8.562900e-34</td><td>-2.618410</td><td>0.247</td><td>0.847</td><td>1.419814e-29</td></tr>
	<tr><th scope=row>Ly6c2</th><td>5.923191e-21</td><td>-2.458466</td><td>0.010</td><td>0.409</td><td>9.821243e-17</td></tr>
	<tr><th scope=row>Thbs1</th><td>1.887662e-11</td><td>-2.235120</td><td>0.067</td><td>0.343</td><td>3.129933e-07</td></tr>
	<tr><th scope=row>Ms4a4c</th><td>1.244351e-28</td><td>-2.159840</td><td>0.093</td><td>0.650</td><td>2.063259e-24</td></tr>
	<tr><th scope=row>Lst1</th><td>5.044431e-31</td><td>-2.091527</td><td>0.521</td><td>0.927</td><td>8.364172e-27</td></tr>
	<tr><th scope=row>Hp</th><td>4.715458e-29</td><td>-1.883201</td><td>0.010</td><td>0.540</td><td>7.818701e-25</td></tr>
	<tr><th scope=row>Clec4e</th><td>1.391847e-20</td><td>-1.806705</td><td>0.088</td><td>0.540</td><td>2.307822e-16</td></tr>
	<tr><th scope=row>Gsr</th><td>1.896849e-28</td><td>-1.734991</td><td>0.082</td><td>0.628</td><td>3.145166e-24</td></tr>
	<tr><th scope=row>Ccrl2</th><td>4.921039e-04</td><td>-1.589357</td><td>0.366</td><td>0.496</td><td>1.000000e+00</td></tr>
	<tr><th scope=row>Ms4a6c</th><td>1.221448e-29</td><td>-1.586063</td><td>0.536</td><td>0.949</td><td>2.025282e-25</td></tr>
	<tr><th scope=row>Il1rn</th><td>1.228838e-09</td><td>-1.547629</td><td>0.175</td><td>0.445</td><td>2.037535e-05</td></tr>
	<tr><th scope=row>Ifitm3</th><td>1.508346e-23</td><td>-1.511220</td><td>0.928</td><td>0.985</td><td>2.500989e-19</td></tr>
	<tr><th scope=row>Hes1</th><td>5.591353e-09</td><td>-1.470905</td><td>0.098</td><td>0.350</td><td>9.271022e-05</td></tr>
	<tr><th scope=row>Osm</th><td>9.147274e-12</td><td>-1.469938</td><td>0.139</td><td>0.467</td><td>1.516710e-07</td></tr>
	<tr><th scope=row>Srgn</th><td>1.872576e-18</td><td>-1.463636</td><td>0.830</td><td>0.978</td><td>3.104918e-14</td></tr>
	<tr><th scope=row>Ptgs2</th><td>3.433402e-08</td><td>-1.452053</td><td>0.057</td><td>0.270</td><td>5.692923e-04</td></tr>
	<tr><th scope=row>Tnfaip2</th><td>1.208355e-19</td><td>-1.447259</td><td>0.098</td><td>0.533</td><td>2.003574e-15</td></tr>
	<tr><th scope=row>B2m</th><td>2.613877e-32</td><td>-1.438790</td><td>0.974</td><td>1.000</td><td>4.334070e-28</td></tr>
	<tr><th scope=row>Prdx5</th><td>1.301159e-24</td><td>-1.437454</td><td>0.670</td><td>0.956</td><td>2.157452e-20</td></tr>
	<tr><th scope=row>Cebpb</th><td>2.569267e-21</td><td>-1.419042</td><td>0.696</td><td>0.949</td><td>4.260101e-17</td></tr>
	<tr><th scope=row>Tgm2</th><td>2.343789e-19</td><td>-1.415635</td><td>0.067</td><td>0.489</td><td>3.886236e-15</td></tr>
	<tr><th scope=row>Nr4a1</th><td>1.581781e-15</td><td>-1.403538</td><td>0.289</td><td>0.686</td><td>2.622751e-11</td></tr>
	<tr><th scope=row>Itgal</th><td>4.890400e-18</td><td>-1.391714</td><td>0.062</td><td>0.453</td><td>8.108772e-14</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>mt-Nd3</th><td>6.956424e-12</td><td>1.022150</td><td>0.938</td><td>0.883</td><td>1.153445e-07</td></tr>
	<tr><th scope=row>Batf3</th><td>2.072078e-11</td><td>1.061276</td><td>0.526</td><td>0.204</td><td>3.435713e-07</td></tr>
	<tr><th scope=row>mt-Co3</th><td>3.922912e-10</td><td>1.083734</td><td>0.995</td><td>1.000</td><td>6.504580e-06</td></tr>
	<tr><th scope=row>Socs6</th><td>9.811051e-10</td><td>1.091445</td><td>0.433</td><td>0.139</td><td>1.626770e-05</td></tr>
	<tr><th scope=row>mt-Nd4</th><td>8.412090e-10</td><td>1.101656</td><td>0.990</td><td>1.000</td><td>1.394809e-05</td></tr>
	<tr><th scope=row>mt-Co2</th><td>2.771530e-11</td><td>1.106304</td><td>1.000</td><td>0.993</td><td>4.595474e-07</td></tr>
	<tr><th scope=row>Dok2</th><td>7.604786e-14</td><td>1.112472</td><td>0.670</td><td>0.336</td><td>1.260950e-09</td></tr>
	<tr><th scope=row>Hspa1b</th><td>2.988819e-05</td><td>1.138198</td><td>0.696</td><td>0.547</td><td>4.955760e-01</td></tr>
	<tr><th scope=row>mt-Nd2</th><td>8.433553e-13</td><td>1.155645</td><td>0.985</td><td>0.956</td><td>1.398367e-08</td></tr>
	<tr><th scope=row>Malat1</th><td>9.148361e-17</td><td>1.156781</td><td>1.000</td><td>1.000</td><td>1.516890e-12</td></tr>
	<tr><th scope=row>Fn1</th><td>3.511531e-17</td><td>1.171349</td><td>0.851</td><td>0.526</td><td>5.822470e-13</td></tr>
	<tr><th scope=row>Olfm1</th><td>1.387567e-15</td><td>1.172815</td><td>0.593</td><td>0.190</td><td>2.300725e-11</td></tr>
	<tr><th scope=row>Jund</th><td>7.222115e-13</td><td>1.186197</td><td>0.954</td><td>0.832</td><td>1.197499e-08</td></tr>
	<tr><th scope=row>Plxnd1</th><td>2.204314e-16</td><td>1.188227</td><td>0.649</td><td>0.219</td><td>3.654973e-12</td></tr>
	<tr><th scope=row>H2-Ab1</th><td>3.393491e-22</td><td>1.240769</td><td>1.000</td><td>0.985</td><td>5.626748e-18</td></tr>
	<tr><th scope=row>H2-Eb1</th><td>5.447739e-21</td><td>1.264822</td><td>1.000</td><td>0.905</td><td>9.032896e-17</td></tr>
	<tr><th scope=row>Pltp</th><td>8.925093e-22</td><td>1.297255</td><td>0.938</td><td>0.693</td><td>1.479870e-17</td></tr>
	<tr><th scope=row>Gm26917</th><td>8.668553e-06</td><td>1.333668</td><td>0.680</td><td>0.533</td><td>1.437333e-01</td></tr>
	<tr><th scope=row>Clec4b1</th><td>1.588642e-16</td><td>1.341989</td><td>0.567</td><td>0.139</td><td>2.634127e-12</td></tr>
	<tr><th scope=row>mt-Co1</th><td>1.968073e-15</td><td>1.360639</td><td>1.000</td><td>1.000</td><td>3.263262e-11</td></tr>
	<tr><th scope=row>Gm42418</th><td>6.146713e-15</td><td>1.387062</td><td>0.954</td><td>0.708</td><td>1.019186e-10</td></tr>
	<tr><th scope=row>Cd81</th><td>1.589229e-22</td><td>1.390493</td><td>0.876</td><td>0.555</td><td>2.635100e-18</td></tr>
	<tr><th scope=row>mt-Nd1</th><td>5.803402e-14</td><td>1.402600</td><td>0.995</td><td>0.993</td><td>9.622621e-10</td></tr>
	<tr><th scope=row>Ahnak</th><td>8.446160e-22</td><td>1.465263</td><td>0.974</td><td>0.766</td><td>1.400458e-17</td></tr>
	<tr><th scope=row>Mgl2</th><td>1.259125e-15</td><td>1.677039</td><td>0.722</td><td>0.358</td><td>2.087755e-11</td></tr>
	<tr><th scope=row>Furin</th><td>1.392064e-12</td><td>1.918956</td><td>0.773</td><td>0.453</td><td>2.308181e-08</td></tr>
	<tr><th scope=row>Fcrls</th><td>1.590678e-33</td><td>2.365155</td><td>0.866</td><td>0.314</td><td>2.637504e-29</td></tr>
	<tr><th scope=row>Lpl</th><td>1.930156e-31</td><td>3.401876</td><td>0.789</td><td>0.190</td><td>3.200392e-27</td></tr>
	<tr><th scope=row>Lyz1</th><td>8.315889e-23</td><td>3.439074</td><td>0.722</td><td>0.226</td><td>1.378858e-18</td></tr>
	<tr><th scope=row>Retnla</th><td>2.319146e-29</td><td>4.867743</td><td>0.866</td><td>0.489</td><td>3.845376e-25</td></tr>
</tbody>
</table>




```R
markers %>%
arrange(desc(avg_log2FC))
```


<table class="dataframe">
<caption>A data.frame: 1545 × 5</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Retnla</th><td>2.319146e-29</td><td>4.867743</td><td>0.866</td><td>0.489</td><td>3.845376e-25</td></tr>
	<tr><th scope=row>Lyz1</th><td>8.315889e-23</td><td>3.439074</td><td>0.722</td><td>0.226</td><td>1.378858e-18</td></tr>
	<tr><th scope=row>Lpl</th><td>1.930156e-31</td><td>3.401876</td><td>0.789</td><td>0.190</td><td>3.200392e-27</td></tr>
	<tr><th scope=row>Fcrls</th><td>1.590678e-33</td><td>2.365155</td><td>0.866</td><td>0.314</td><td>2.637504e-29</td></tr>
	<tr><th scope=row>Furin</th><td>1.392064e-12</td><td>1.918956</td><td>0.773</td><td>0.453</td><td>2.308181e-08</td></tr>
	<tr><th scope=row>Mgl2</th><td>1.259125e-15</td><td>1.677039</td><td>0.722</td><td>0.358</td><td>2.087755e-11</td></tr>
	<tr><th scope=row>Ahnak</th><td>8.446160e-22</td><td>1.465263</td><td>0.974</td><td>0.766</td><td>1.400458e-17</td></tr>
	<tr><th scope=row>mt-Nd1</th><td>5.803402e-14</td><td>1.402600</td><td>0.995</td><td>0.993</td><td>9.622621e-10</td></tr>
	<tr><th scope=row>Cd81</th><td>1.589229e-22</td><td>1.390493</td><td>0.876</td><td>0.555</td><td>2.635100e-18</td></tr>
	<tr><th scope=row>Gm42418</th><td>6.146713e-15</td><td>1.387062</td><td>0.954</td><td>0.708</td><td>1.019186e-10</td></tr>
	<tr><th scope=row>mt-Co1</th><td>1.968073e-15</td><td>1.360639</td><td>1.000</td><td>1.000</td><td>3.263262e-11</td></tr>
	<tr><th scope=row>Clec4b1</th><td>1.588642e-16</td><td>1.341989</td><td>0.567</td><td>0.139</td><td>2.634127e-12</td></tr>
	<tr><th scope=row>Gm26917</th><td>8.668553e-06</td><td>1.333668</td><td>0.680</td><td>0.533</td><td>1.437333e-01</td></tr>
	<tr><th scope=row>Pltp</th><td>8.925093e-22</td><td>1.297255</td><td>0.938</td><td>0.693</td><td>1.479870e-17</td></tr>
	<tr><th scope=row>H2-Eb1</th><td>5.447739e-21</td><td>1.264822</td><td>1.000</td><td>0.905</td><td>9.032896e-17</td></tr>
	<tr><th scope=row>H2-Ab1</th><td>3.393491e-22</td><td>1.240769</td><td>1.000</td><td>0.985</td><td>5.626748e-18</td></tr>
	<tr><th scope=row>Plxnd1</th><td>2.204314e-16</td><td>1.188227</td><td>0.649</td><td>0.219</td><td>3.654973e-12</td></tr>
	<tr><th scope=row>Jund</th><td>7.222115e-13</td><td>1.186197</td><td>0.954</td><td>0.832</td><td>1.197499e-08</td></tr>
	<tr><th scope=row>Olfm1</th><td>1.387567e-15</td><td>1.172815</td><td>0.593</td><td>0.190</td><td>2.300725e-11</td></tr>
	<tr><th scope=row>Fn1</th><td>3.511531e-17</td><td>1.171349</td><td>0.851</td><td>0.526</td><td>5.822470e-13</td></tr>
	<tr><th scope=row>Malat1</th><td>9.148361e-17</td><td>1.156781</td><td>1.000</td><td>1.000</td><td>1.516890e-12</td></tr>
	<tr><th scope=row>mt-Nd2</th><td>8.433553e-13</td><td>1.155645</td><td>0.985</td><td>0.956</td><td>1.398367e-08</td></tr>
	<tr><th scope=row>Hspa1b</th><td>2.988819e-05</td><td>1.138198</td><td>0.696</td><td>0.547</td><td>4.955760e-01</td></tr>
	<tr><th scope=row>Dok2</th><td>7.604786e-14</td><td>1.112472</td><td>0.670</td><td>0.336</td><td>1.260950e-09</td></tr>
	<tr><th scope=row>mt-Co2</th><td>2.771530e-11</td><td>1.106304</td><td>1.000</td><td>0.993</td><td>4.595474e-07</td></tr>
	<tr><th scope=row>mt-Nd4</th><td>8.412090e-10</td><td>1.101656</td><td>0.990</td><td>1.000</td><td>1.394809e-05</td></tr>
	<tr><th scope=row>Socs6</th><td>9.811051e-10</td><td>1.091445</td><td>0.433</td><td>0.139</td><td>1.626770e-05</td></tr>
	<tr><th scope=row>mt-Co3</th><td>3.922912e-10</td><td>1.083734</td><td>0.995</td><td>1.000</td><td>6.504580e-06</td></tr>
	<tr><th scope=row>Batf3</th><td>2.072078e-11</td><td>1.061276</td><td>0.526</td><td>0.204</td><td>3.435713e-07</td></tr>
	<tr><th scope=row>mt-Nd3</th><td>6.956424e-12</td><td>1.022150</td><td>0.938</td><td>0.883</td><td>1.153445e-07</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>Itgal</th><td>4.890400e-18</td><td>-1.391714</td><td>0.062</td><td>0.453</td><td>8.108772e-14</td></tr>
	<tr><th scope=row>Nr4a1</th><td>1.581781e-15</td><td>-1.403538</td><td>0.289</td><td>0.686</td><td>2.622751e-11</td></tr>
	<tr><th scope=row>Tgm2</th><td>2.343789e-19</td><td>-1.415635</td><td>0.067</td><td>0.489</td><td>3.886236e-15</td></tr>
	<tr><th scope=row>Cebpb</th><td>2.569267e-21</td><td>-1.419042</td><td>0.696</td><td>0.949</td><td>4.260101e-17</td></tr>
	<tr><th scope=row>Prdx5</th><td>1.301159e-24</td><td>-1.437454</td><td>0.670</td><td>0.956</td><td>2.157452e-20</td></tr>
	<tr><th scope=row>B2m</th><td>2.613877e-32</td><td>-1.438790</td><td>0.974</td><td>1.000</td><td>4.334070e-28</td></tr>
	<tr><th scope=row>Tnfaip2</th><td>1.208355e-19</td><td>-1.447259</td><td>0.098</td><td>0.533</td><td>2.003574e-15</td></tr>
	<tr><th scope=row>Ptgs2</th><td>3.433402e-08</td><td>-1.452053</td><td>0.057</td><td>0.270</td><td>5.692923e-04</td></tr>
	<tr><th scope=row>Srgn</th><td>1.872576e-18</td><td>-1.463636</td><td>0.830</td><td>0.978</td><td>3.104918e-14</td></tr>
	<tr><th scope=row>Osm</th><td>9.147274e-12</td><td>-1.469938</td><td>0.139</td><td>0.467</td><td>1.516710e-07</td></tr>
	<tr><th scope=row>Hes1</th><td>5.591353e-09</td><td>-1.470905</td><td>0.098</td><td>0.350</td><td>9.271022e-05</td></tr>
	<tr><th scope=row>Ifitm3</th><td>1.508346e-23</td><td>-1.511220</td><td>0.928</td><td>0.985</td><td>2.500989e-19</td></tr>
	<tr><th scope=row>Il1rn</th><td>1.228838e-09</td><td>-1.547629</td><td>0.175</td><td>0.445</td><td>2.037535e-05</td></tr>
	<tr><th scope=row>Ms4a6c</th><td>1.221448e-29</td><td>-1.586063</td><td>0.536</td><td>0.949</td><td>2.025282e-25</td></tr>
	<tr><th scope=row>Ccrl2</th><td>4.921039e-04</td><td>-1.589357</td><td>0.366</td><td>0.496</td><td>1.000000e+00</td></tr>
	<tr><th scope=row>Gsr</th><td>1.896849e-28</td><td>-1.734991</td><td>0.082</td><td>0.628</td><td>3.145166e-24</td></tr>
	<tr><th scope=row>Clec4e</th><td>1.391847e-20</td><td>-1.806705</td><td>0.088</td><td>0.540</td><td>2.307822e-16</td></tr>
	<tr><th scope=row>Hp</th><td>4.715458e-29</td><td>-1.883201</td><td>0.010</td><td>0.540</td><td>7.818701e-25</td></tr>
	<tr><th scope=row>Lst1</th><td>5.044431e-31</td><td>-2.091527</td><td>0.521</td><td>0.927</td><td>8.364172e-27</td></tr>
	<tr><th scope=row>Ms4a4c</th><td>1.244351e-28</td><td>-2.159840</td><td>0.093</td><td>0.650</td><td>2.063259e-24</td></tr>
	<tr><th scope=row>Thbs1</th><td>1.887662e-11</td><td>-2.235120</td><td>0.067</td><td>0.343</td><td>3.129933e-07</td></tr>
	<tr><th scope=row>Ly6c2</th><td>5.923191e-21</td><td>-2.458466</td><td>0.010</td><td>0.409</td><td>9.821243e-17</td></tr>
	<tr><th scope=row>Gngt2</th><td>8.562900e-34</td><td>-2.618410</td><td>0.247</td><td>0.847</td><td>1.419814e-29</td></tr>
	<tr><th scope=row>Msrb1</th><td>6.761738e-43</td><td>-2.660776</td><td>0.366</td><td>0.949</td><td>1.121164e-38</td></tr>
	<tr><th scope=row>Chil3</th><td>2.666928e-18</td><td>-2.763455</td><td>0.010</td><td>0.365</td><td>4.422034e-14</td></tr>
	<tr><th scope=row>Slpi</th><td>4.823473e-08</td><td>-2.909235</td><td>0.036</td><td>0.226</td><td>7.997800e-04</td></tr>
	<tr><th scope=row>Il1b</th><td>1.790189e-17</td><td>-3.204045</td><td>0.268</td><td>0.664</td><td>2.968313e-13</td></tr>
	<tr><th scope=row>Plac8</th><td>7.070759e-43</td><td>-4.635609</td><td>0.041</td><td>0.759</td><td>1.172403e-38</td></tr>
	<tr><th scope=row>S100a8</th><td>3.728914e-02</td><td>-4.730024</td><td>0.088</td><td>0.153</td><td>1.000000e+00</td></tr>
	<tr><th scope=row>S100a9</th><td>1.864518e-03</td><td>-4.756436</td><td>0.041</td><td>0.131</td><td>1.000000e+00</td></tr>
</tbody>
</table>




```R
DimPlot(immune.combined, group.by = 'integrated_snn_res.1.1',label = T) +NoLegend()
DimPlot(immune.combined, group.by = 'group') 
```


    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_10_0.png)
    



    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_10_1.png)
    



```R
Idents(immune.combined) <- 'integrated_snn_res.1.1'
new.cluster.ids<-c("Resident like Macrophages",'Resident like Macrophages','Inflammatory Macrophages','Inflammatory Macrophages','Resident like Macrophages',
                  'Inflammatory Macrophages','Resident like Macrophages','Resident like Macrophages','Resident like Macrophages','Foam cells',
                  'Resident like Macrophages','Inflammatory Macrophages','IFN induced Macrophages','Dendritic cells_1','T cells',
                   'Dendritic cells_1','Monocytes','Dendritic cells_2','Fibroblasts','Endothelial cells',
                  'Hematopoietic stem cells','Fibroblasts','Fibroblasts')
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined,reduction="umap",label=T)+NoLegend()
```


    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_11_0.png)
    



```R
immune.combined[["Annotation_Formal"]] <- Idents(immune.combined)
```


```R
DimPlot(immune.combined, group.by = 'Annotation_Formal',label = T) +NoLegend()
```


    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_13_0.png)
    



```R
Idents(immune.combined) <- 'integrated_snn_res.1.1'
new.cluster.ids<-c("Res Mφs",'Res Mφs','Inflam Mφs','Inflam Mφs','Res Mφs',
                  'Inflam Mφs','Res Mφs','Res Mφs','Res Mφs','Foam cells',
                  'Res Mφs','Inflam Mφs','IFN Mφs','DCs_1','TCs',
                   'DCs_1','Monocytes','DCs_2','FBs','ECs',
                  'HSCs','FBs','FBs')
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined,reduction="umap",label=T)+NoLegend()
```


    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_14_0.png)
    



```R
immune.combined[["Annotation_Brief"]] <- Idents(immune.combined)
```


```R
DimPlot(immune.combined, group.by = 'Annotation_Brief',label = T) +NoLegend()
```


    
![png](Step7.3_Seurat_Annotation_files/Step7.3_Seurat_Annotation_16_0.png)
    



```R
saveRDS(immune.combined, file = 'Step4_Cluster_Annotatoin.Rds')
```
