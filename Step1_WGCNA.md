```R
###
rm(list = ls())
options(stringsAsFactors = F)
```


```R
library(WGCNA) 
library(DESeq2)
library(tidyverse)
```

    Loading required package: dynamicTreeCut
    
    Loading required package: fastcluster
    
    
    Attaching package: ‘fastcluster’
    
    
    The following object is masked from ‘package:stats’:
    
        hclust
    
    
    
    
    
    Attaching package: ‘WGCNA’
    
    
    The following object is masked from ‘package:stats’:
    
        cor
    
    
    Loading required package: S4Vectors
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    Loading required package: parallel
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:parallel’:
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min
    
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following object is masked from ‘package:base’:
    
        expand.grid
    
    
    Loading required package: IRanges
    
    Loading required package: GenomicRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: ‘Biobase’
    
    
    The following object is masked from ‘package:MatrixGenerics’:
    
        rowMedians
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        anyMissing, rowMedians
    
    
    Warning message in system("timedatectl", intern = TRUE):
    “running command 'timedatectl' had status 1”
    ── [1mAttaching packages[22m ─────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
    
    [32m✔[39m [34mggplot2[39m 3.4.0     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.5     [32m✔[39m [34mdplyr  [39m 1.0.7
    [32m✔[39m [34mtidyr  [39m 1.1.4     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 2.0.2     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mcollapse()[39m   masks [34mIRanges[39m::collapse()
    [31m✖[39m [34mdplyr[39m::[32mcombine()[39m    masks [34mBiobase[39m::combine(), [34mBiocGenerics[39m::combine()
    [31m✖[39m [34mdplyr[39m::[32mcount()[39m      masks [34mmatrixStats[39m::count()
    [31m✖[39m [34mdplyr[39m::[32mdesc()[39m       masks [34mIRanges[39m::desc()
    [31m✖[39m [34mtidyr[39m::[32mexpand()[39m     masks [34mS4Vectors[39m::expand()
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m     masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mfirst()[39m      masks [34mS4Vectors[39m::first()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m        masks [34mstats[39m::lag()
    [31m✖[39m [34mggplot2[39m::[32mPosition()[39m masks [34mBiocGenerics[39m::Position(), [34mbase[39m::Position()
    [31m✖[39m [34mpurrr[39m::[32mreduce()[39m     masks [34mGenomicRanges[39m::reduce(), [34mIRanges[39m::reduce()
    [31m✖[39m [34mdplyr[39m::[32mrename()[39m     masks [34mS4Vectors[39m::rename()
    [31m✖[39m [34mdplyr[39m::[32mslice()[39m      masks [34mIRanges[39m::slice()
    



```R
package.version('WGCNA')
```


'1.71'



```R
exp_dat <-  read_csv('./RNA_counts.csv') %>%
 column_to_rownames(var ='gene')
```

    [1mRows: [22m[34m15871[39m [1mColumns: [22m[34m17[39m
    
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m  (1): gene
    [32mdbl[39m (16): H_C_1, H_C_2, H_C_3, H_C_4, H_O_1, H_O_2, H_O_3, H_O_4, L_C_1, L_C...
    
    
    [36mℹ[39m Use [30m[47m[30m[47m`spec()`[47m[30m[49m[39m to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set [30m[47m[30m[47m`show_col_types = FALSE`[47m[30m[49m[39m to quiet this message.
    



```R
high_glucose <- c(rep(1,8),rep(0,8))
oxldl <- c(rep(0,4),rep(1,4),rep(0,4),rep(1,4))
high_glucose
oxldl
metadata <- data.frame(high_glucose,oxldl)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0</li><li>0</li><li>0</li><li>0</li><li>1</li><li>1</li><li>1</li><li>1</li><li>0</li><li>0</li><li>0</li><li>0</li><li>1</li><li>1</li><li>1</li><li>1</li></ol>




```R
metadata$high_glucose <- as.factor(metadata$high_glucose)
metadata$oxldl <- as.factor(metadata$oxldl)
```


```R
# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = exp_dat, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)
```

    converting counts to integer mode
    



```R
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)
```


```R
# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data
```


```R
normalized_counts
```


<table class="dataframe">
<caption>A matrix: 16 × 15871 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Gpnmb</th><th scope=col>Psap</th><th scope=col>Fth1</th><th scope=col>Ctsd</th><th scope=col>Ctsb</th><th scope=col>mt-Co1</th><th scope=col>Lyz2</th><th scope=col>Ctss</th><th scope=col>Mpeg1</th><th scope=col>Ftl1</th><th scope=col>⋯</th><th scope=col>Epo</th><th scope=col>Kif28</th><th scope=col>AL607142.1</th><th scope=col>Hist1h2ah</th><th scope=col>Gm20538</th><th scope=col>Gm11084</th><th scope=col>Gm20662</th><th scope=col>Gm49342</th><th scope=col>Exosc6</th><th scope=col>Gm28053</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>19.84392</td><td>19.66701</td><td>19.32427</td><td>19.60235</td><td>19.31819</td><td>18.19933</td><td>19.02063</td><td>18.78647</td><td>18.62897</td><td>18.28188</td><td>⋯</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>6.964300</td><td>6.964300</td><td>7.099573</td><td>7.099573</td></tr>
	<tr><th scope=row>H_C_2</th><td>19.82143</td><td>19.68751</td><td>19.26940</td><td>19.57475</td><td>19.32144</td><td>18.37896</td><td>18.99820</td><td>18.75755</td><td>18.62014</td><td>18.15879</td><td>⋯</td><td>7.096656</td><td>7.096656</td><td>6.964300</td><td>6.964300</td><td>7.096656</td><td>7.096656</td><td>6.964300</td><td>7.096656</td><td>7.096656</td><td>6.964300</td></tr>
	<tr><th scope=row>H_C_3</th><td>19.79431</td><td>19.64258</td><td>19.25138</td><td>19.53206</td><td>19.30767</td><td>18.29319</td><td>19.04634</td><td>18.81283</td><td>18.62392</td><td>18.13978</td><td>⋯</td><td>7.088083</td><td>7.088083</td><td>7.088083</td><td>6.964300</td><td>7.088083</td><td>6.964300</td><td>6.964300</td><td>7.088083</td><td>7.088083</td><td>6.964300</td></tr>
	<tr><th scope=row>H_C_4</th><td>19.80935</td><td>19.65189</td><td>19.25727</td><td>19.54011</td><td>19.29872</td><td>18.39850</td><td>19.00260</td><td>18.77303</td><td>18.61877</td><td>18.09527</td><td>⋯</td><td>6.964300</td><td>7.096213</td><td>7.096213</td><td>6.964300</td><td>7.096213</td><td>7.096213</td><td>7.096213</td><td>6.964300</td><td>7.096213</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_1</th><td>19.40186</td><td>19.48225</td><td>19.36663</td><td>19.45605</td><td>19.10526</td><td>18.35323</td><td>19.09379</td><td>18.23198</td><td>18.15746</td><td>18.34140</td><td>⋯</td><td>7.089577</td><td>7.089577</td><td>7.089577</td><td>7.089577</td><td>7.089577</td><td>6.964300</td><td>6.964300</td><td>7.089577</td><td>7.089577</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_2</th><td>19.30420</td><td>19.44577</td><td>19.33320</td><td>19.40336</td><td>19.09562</td><td>18.32884</td><td>19.19804</td><td>18.48864</td><td>18.22323</td><td>18.20287</td><td>⋯</td><td>7.089424</td><td>7.089424</td><td>7.089424</td><td>6.964300</td><td>6.964300</td><td>7.089424</td><td>7.089424</td><td>7.089424</td><td>7.089424</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_3</th><td>19.38136</td><td>19.55646</td><td>19.42910</td><td>19.42724</td><td>19.13191</td><td>18.36379</td><td>19.09936</td><td>18.34826</td><td>18.14207</td><td>18.46552</td><td>⋯</td><td>7.091864</td><td>7.091864</td><td>7.091864</td><td>7.091864</td><td>7.091864</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.091864</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_4</th><td>19.41995</td><td>19.50978</td><td>19.35323</td><td>19.41914</td><td>19.14318</td><td>18.42633</td><td>19.07696</td><td>18.24008</td><td>18.10835</td><td>18.38856</td><td>⋯</td><td>7.085405</td><td>7.085405</td><td>6.964300</td><td>7.085405</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.085405</td><td>7.085405</td><td>7.085405</td></tr>
	<tr><th scope=row>L_C_1</th><td>19.76070</td><td>19.75923</td><td>19.25629</td><td>19.58424</td><td>19.30039</td><td>18.13375</td><td>18.97041</td><td>18.72967</td><td>18.44761</td><td>18.10005</td><td>⋯</td><td>7.099935</td><td>7.099935</td><td>7.099935</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.099935</td><td>7.099935</td><td>7.099935</td></tr>
	<tr><th scope=row>L_C_2</th><td>19.67348</td><td>19.75032</td><td>19.24019</td><td>19.51053</td><td>19.26887</td><td>18.33763</td><td>18.96034</td><td>18.77695</td><td>18.43574</td><td>18.06239</td><td>⋯</td><td>6.964300</td><td>7.092951</td><td>7.092951</td><td>7.092951</td><td>6.964300</td><td>7.092951</td><td>6.964300</td><td>6.964300</td><td>7.092951</td><td>6.964300</td></tr>
	<tr><th scope=row>L_C_3</th><td>19.69776</td><td>19.72907</td><td>19.22201</td><td>19.54604</td><td>19.26165</td><td>18.23735</td><td>18.86380</td><td>18.63988</td><td>18.33978</td><td>18.10834</td><td>⋯</td><td>7.097444</td><td>6.964300</td><td>7.097444</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.097444</td><td>6.964300</td></tr>
	<tr><th scope=row>L_C_4</th><td>20.10351</td><td>19.71374</td><td>19.63206</td><td>19.04233</td><td>19.19639</td><td>19.30477</td><td>19.00606</td><td>18.48129</td><td>18.23226</td><td>17.48928</td><td>⋯</td><td>7.088116</td><td>7.088116</td><td>7.088116</td><td>6.964300</td><td>6.964300</td><td>7.088116</td><td>6.964300</td><td>7.088116</td><td>7.088116</td><td>6.964300</td></tr>
	<tr><th scope=row>L_O_1</th><td>19.30021</td><td>19.64602</td><td>19.38676</td><td>19.47236</td><td>19.09137</td><td>18.34954</td><td>19.02381</td><td>18.12189</td><td>18.00268</td><td>18.25962</td><td>⋯</td><td>7.098663</td><td>6.964300</td><td>7.098663</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.098663</td><td>7.098663</td><td>7.098663</td><td>7.098663</td></tr>
	<tr><th scope=row>L_O_2</th><td>19.34710</td><td>19.62101</td><td>19.38559</td><td>19.48709</td><td>19.11748</td><td>18.26185</td><td>18.97049</td><td>18.13506</td><td>18.01441</td><td>18.36307</td><td>⋯</td><td>7.089137</td><td>7.089137</td><td>7.089137</td><td>6.964300</td><td>6.964300</td><td>7.089137</td><td>7.089137</td><td>7.089137</td><td>7.089137</td><td>7.089137</td></tr>
	<tr><th scope=row>L_O_3</th><td>19.31497</td><td>19.61604</td><td>19.35164</td><td>19.46266</td><td>19.09135</td><td>18.30043</td><td>18.97454</td><td>18.08274</td><td>18.01496</td><td>18.33809</td><td>⋯</td><td>6.964300</td><td>7.095753</td><td>6.964300</td><td>7.095753</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.095753</td><td>7.095753</td><td>7.095753</td></tr>
	<tr><th scope=row>L_O_4</th><td>19.30582</td><td>19.60616</td><td>19.36783</td><td>19.44114</td><td>19.10754</td><td>18.33136</td><td>18.96661</td><td>18.05261</td><td>18.01056</td><td>18.31798</td><td>⋯</td><td>7.091786</td><td>6.964300</td><td>7.091786</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.091786</td><td>7.091786</td><td>6.964300</td></tr>
</tbody>
</table>




```R
exp_mt <-  as.data.frame(normalized_counts)
```


```R
##(2)判断数据质量--缺失值
gsg = goodSamplesGenes(exp_mt)
gsg$allOK
# [1] TRUE
```

     Flagging genes and samples with too many missing values...
      ..step 1



TRUE



```R
##(3)判断数据质量--离群点样本
sampleTree = stats::hclust(dist(exp_mt), method = "average")
```


```R
pdf("01.samples_cluster_tree.pdf",width=25,height=8)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="")
cutHeight <- 30
abline(h=cutHeight,col="red")
dev.off()
```


<strong>png:</strong> 2



```R
##如下图所示，存在一个显著离群点；剔除掉
clust = cutreeStatic(sampleTree, cutHeight = 30, minSize = 2)
table(clust)
# clust
# 0   1 
# 1 134 
keepSamples = (clust!=0)
exp_mt_f = exp_mt[keepSamples, ]
```


    clust
    0 1 2 
    1 8 7 



```R
exp_mt_f
```


<table class="dataframe">
<caption>A data.frame: 15 × 15871</caption>
<thead>
	<tr><th></th><th scope=col>Gpnmb</th><th scope=col>Psap</th><th scope=col>Fth1</th><th scope=col>Ctsd</th><th scope=col>Ctsb</th><th scope=col>mt-Co1</th><th scope=col>Lyz2</th><th scope=col>Ctss</th><th scope=col>Mpeg1</th><th scope=col>Ftl1</th><th scope=col>⋯</th><th scope=col>Epo</th><th scope=col>Kif28</th><th scope=col>AL607142.1</th><th scope=col>Hist1h2ah</th><th scope=col>Gm20538</th><th scope=col>Gm11084</th><th scope=col>Gm20662</th><th scope=col>Gm49342</th><th scope=col>Exosc6</th><th scope=col>Gm28053</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>19.84392</td><td>19.66701</td><td>19.32427</td><td>19.60235</td><td>19.31819</td><td>18.19933</td><td>19.02063</td><td>18.78647</td><td>18.62897</td><td>18.28188</td><td>⋯</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>7.099573</td><td>6.964300</td><td>6.964300</td><td>7.099573</td><td>7.099573</td></tr>
	<tr><th scope=row>H_C_2</th><td>19.82143</td><td>19.68751</td><td>19.26940</td><td>19.57475</td><td>19.32144</td><td>18.37896</td><td>18.99820</td><td>18.75755</td><td>18.62014</td><td>18.15879</td><td>⋯</td><td>7.096656</td><td>7.096656</td><td>6.964300</td><td>6.964300</td><td>7.096656</td><td>7.096656</td><td>6.964300</td><td>7.096656</td><td>7.096656</td><td>6.964300</td></tr>
	<tr><th scope=row>H_C_3</th><td>19.79431</td><td>19.64258</td><td>19.25138</td><td>19.53206</td><td>19.30767</td><td>18.29319</td><td>19.04634</td><td>18.81283</td><td>18.62392</td><td>18.13978</td><td>⋯</td><td>7.088083</td><td>7.088083</td><td>7.088083</td><td>6.964300</td><td>7.088083</td><td>6.964300</td><td>6.964300</td><td>7.088083</td><td>7.088083</td><td>6.964300</td></tr>
	<tr><th scope=row>H_C_4</th><td>19.80935</td><td>19.65189</td><td>19.25727</td><td>19.54011</td><td>19.29872</td><td>18.39850</td><td>19.00260</td><td>18.77303</td><td>18.61877</td><td>18.09527</td><td>⋯</td><td>6.964300</td><td>7.096213</td><td>7.096213</td><td>6.964300</td><td>7.096213</td><td>7.096213</td><td>7.096213</td><td>6.964300</td><td>7.096213</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_1</th><td>19.40186</td><td>19.48225</td><td>19.36663</td><td>19.45605</td><td>19.10526</td><td>18.35323</td><td>19.09379</td><td>18.23198</td><td>18.15746</td><td>18.34140</td><td>⋯</td><td>7.089577</td><td>7.089577</td><td>7.089577</td><td>7.089577</td><td>7.089577</td><td>6.964300</td><td>6.964300</td><td>7.089577</td><td>7.089577</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_2</th><td>19.30420</td><td>19.44577</td><td>19.33320</td><td>19.40336</td><td>19.09562</td><td>18.32884</td><td>19.19804</td><td>18.48864</td><td>18.22323</td><td>18.20287</td><td>⋯</td><td>7.089424</td><td>7.089424</td><td>7.089424</td><td>6.964300</td><td>6.964300</td><td>7.089424</td><td>7.089424</td><td>7.089424</td><td>7.089424</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_3</th><td>19.38136</td><td>19.55646</td><td>19.42910</td><td>19.42724</td><td>19.13191</td><td>18.36379</td><td>19.09936</td><td>18.34826</td><td>18.14207</td><td>18.46552</td><td>⋯</td><td>7.091864</td><td>7.091864</td><td>7.091864</td><td>7.091864</td><td>7.091864</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.091864</td><td>6.964300</td></tr>
	<tr><th scope=row>H_O_4</th><td>19.41995</td><td>19.50978</td><td>19.35323</td><td>19.41914</td><td>19.14318</td><td>18.42633</td><td>19.07696</td><td>18.24008</td><td>18.10835</td><td>18.38856</td><td>⋯</td><td>7.085405</td><td>7.085405</td><td>6.964300</td><td>7.085405</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.085405</td><td>7.085405</td><td>7.085405</td></tr>
	<tr><th scope=row>L_C_1</th><td>19.76070</td><td>19.75923</td><td>19.25629</td><td>19.58424</td><td>19.30039</td><td>18.13375</td><td>18.97041</td><td>18.72967</td><td>18.44761</td><td>18.10005</td><td>⋯</td><td>7.099935</td><td>7.099935</td><td>7.099935</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.099935</td><td>7.099935</td><td>7.099935</td></tr>
	<tr><th scope=row>L_C_2</th><td>19.67348</td><td>19.75032</td><td>19.24019</td><td>19.51053</td><td>19.26887</td><td>18.33763</td><td>18.96034</td><td>18.77695</td><td>18.43574</td><td>18.06239</td><td>⋯</td><td>6.964300</td><td>7.092951</td><td>7.092951</td><td>7.092951</td><td>6.964300</td><td>7.092951</td><td>6.964300</td><td>6.964300</td><td>7.092951</td><td>6.964300</td></tr>
	<tr><th scope=row>L_C_3</th><td>19.69776</td><td>19.72907</td><td>19.22201</td><td>19.54604</td><td>19.26165</td><td>18.23735</td><td>18.86380</td><td>18.63988</td><td>18.33978</td><td>18.10834</td><td>⋯</td><td>7.097444</td><td>6.964300</td><td>7.097444</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.097444</td><td>6.964300</td></tr>
	<tr><th scope=row>L_O_1</th><td>19.30021</td><td>19.64602</td><td>19.38676</td><td>19.47236</td><td>19.09137</td><td>18.34954</td><td>19.02381</td><td>18.12189</td><td>18.00268</td><td>18.25962</td><td>⋯</td><td>7.098663</td><td>6.964300</td><td>7.098663</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.098663</td><td>7.098663</td><td>7.098663</td><td>7.098663</td></tr>
	<tr><th scope=row>L_O_2</th><td>19.34710</td><td>19.62101</td><td>19.38559</td><td>19.48709</td><td>19.11748</td><td>18.26185</td><td>18.97049</td><td>18.13506</td><td>18.01441</td><td>18.36307</td><td>⋯</td><td>7.089137</td><td>7.089137</td><td>7.089137</td><td>6.964300</td><td>6.964300</td><td>7.089137</td><td>7.089137</td><td>7.089137</td><td>7.089137</td><td>7.089137</td></tr>
	<tr><th scope=row>L_O_3</th><td>19.31497</td><td>19.61604</td><td>19.35164</td><td>19.46266</td><td>19.09135</td><td>18.30043</td><td>18.97454</td><td>18.08274</td><td>18.01496</td><td>18.33809</td><td>⋯</td><td>6.964300</td><td>7.095753</td><td>6.964300</td><td>7.095753</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.095753</td><td>7.095753</td><td>7.095753</td></tr>
	<tr><th scope=row>L_O_4</th><td>19.30582</td><td>19.60616</td><td>19.36783</td><td>19.44114</td><td>19.10754</td><td>18.33136</td><td>18.96661</td><td>18.05261</td><td>18.01056</td><td>18.31798</td><td>⋯</td><td>7.091786</td><td>6.964300</td><td>7.091786</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>6.964300</td><td>7.091786</td><td>7.091786</td><td>6.964300</td></tr>
</tbody>
</table>




```R
rownames(exp_mt_f) 
high_glucose <- c(rep(1,8),rep(0,7))
oxldl <- c(rep(0,4),rep(1,4),rep(0,3),rep(1,4))
high_glucose
oxldl
trait_dat_f <- data.frame(high_glucose,oxldl)
rownames(trait_dat_f) <- rownames(exp_mt_f) 
trait_dat_f
identical(rownames(exp_mt_f), rownames(trait_dat_f))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'H_C_1'</li><li>'H_C_2'</li><li>'H_C_3'</li><li>'H_C_4'</li><li>'H_O_1'</li><li>'H_O_2'</li><li>'H_O_3'</li><li>'H_O_4'</li><li>'L_C_1'</li><li>'L_C_2'</li><li>'L_C_3'</li><li>'L_O_1'</li><li>'L_O_2'</li><li>'L_O_3'</li><li>'L_O_4'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0</li><li>0</li><li>0</li><li>0</li><li>1</li><li>1</li><li>1</li><li>1</li><li>0</li><li>0</li><li>0</li><li>1</li><li>1</li><li>1</li><li>1</li></ol>




<table class="dataframe">
<caption>A data.frame: 15 × 2</caption>
<thead>
	<tr><th></th><th scope=col>high_glucose</th><th scope=col>oxldl</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_2</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_3</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_4</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_O_1</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>H_O_2</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>H_O_3</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>H_O_4</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>L_C_1</th><td>0</td><td>0</td></tr>
	<tr><th scope=row>L_C_2</th><td>0</td><td>0</td></tr>
	<tr><th scope=row>L_C_3</th><td>0</td><td>0</td></tr>
	<tr><th scope=row>L_O_1</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>L_O_2</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>L_O_3</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>L_O_4</th><td>0</td><td>1</td></tr>
</tbody>
</table>




TRUE



```R
##汇总最终数据
exp_dat = exp_mt_f
dim(exp_dat)
# [1]  134 3600
exp_dat[1:4,1:4]

trait_dat = trait_dat_f
dim(trait_dat)
# [1] 134  25
trait_dat[1:4,1:2]
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>15</li><li>15871</li></ol>




<table class="dataframe">
<caption>A data.frame: 4 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Gpnmb</th><th scope=col>Psap</th><th scope=col>Fth1</th><th scope=col>Ctsd</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>19.84392</td><td>19.66701</td><td>19.32427</td><td>19.60235</td></tr>
	<tr><th scope=row>H_C_2</th><td>19.82143</td><td>19.68751</td><td>19.26940</td><td>19.57475</td></tr>
	<tr><th scope=row>H_C_3</th><td>19.79431</td><td>19.64258</td><td>19.25138</td><td>19.53206</td></tr>
	<tr><th scope=row>H_C_4</th><td>19.80935</td><td>19.65189</td><td>19.25727</td><td>19.54011</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>15</li><li>2</li></ol>




<table class="dataframe">
<caption>A data.frame: 4 × 2</caption>
<thead>
	<tr><th></th><th scope=col>high_glucose</th><th scope=col>oxldl</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_2</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_3</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_4</th><td>1</td><td>0</td></tr>
</tbody>
</table>




```R
sampleTree2 <- hclust(dist(exp_dat),method="average")
traitColors <- numbers2colors(trait_dat,colors=blueWhiteRed(50))
pdf("01.cluster_trait.pdf",width=25,height=10)
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels=names(trait_dat),
                    main="Sample dendrogram and trait heatmap")
dev.off()
```


<strong>png:</strong> 2



```R
powers <- c(c(1:10),seq(from=12,to=30,by=2))
networkType <- "signed"
sft <- pickSoftThreshold(exp_dat,
                         powerVector=powers,
                         networkType="signed",
                         verbose=5,
                         dataIsExpr = TRUE,
                         corFnc = cor)
pdf(file="02.soft_threshold.pdf",width=8,height=4.5)
par(mfrow=c(1,2))
cex1 <- 0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft threshold (power)",
     ylab="Scale free topology model fit, signed R^2",
     type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")

plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab="Soft threshold (power)",ylab="Mean connectivity",
     type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")
dev.off()
```

    pickSoftThreshold: will use block size 2818.
     pickSoftThreshold: calculating connectivity for given powers...
       ..working on genes 1 through 2818 of 15871


    Warning message:
    “executing %dopar% sequentially: no parallel backend registered”


       ..working on genes 2819 through 5636 of 15871
       ..working on genes 5637 through 8454 of 15871
       ..working on genes 8455 through 11272 of 15871
       ..working on genes 11273 through 14090 of 15871
       ..working on genes 14091 through 15871 of 15871
       Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    1      1  0.00414  9.06          0.938  7940.0   7940.00   8080
    2      2  0.37800 -6.59          0.652  4450.0   4400.00   5030
    3      3  0.81100 -3.87          0.878  2700.0   2620.00   3560
    4      4  0.90400 -2.74          0.941  1750.0   1660.00   2710
    5      5  0.93200 -2.23          0.958  1190.0   1090.00   2170
    6      6  0.94800 -1.95          0.968   847.0    742.00   1810
    7      7  0.95600 -1.78          0.972   623.0    519.00   1550
    8      8  0.96300 -1.67          0.975   473.0    372.00   1350
    9      9  0.96800 -1.60          0.979   368.0    271.00   1200
    10    10  0.97000 -1.54          0.980   293.0    200.00   1080
    11    12  0.96700 -1.47          0.975   196.0    114.00    896
    12    14  0.96700 -1.42          0.977   139.0     68.50    765
    13    16  0.96700 -1.38          0.978   103.0     42.30    667
    14    18  0.96800 -1.35          0.981    79.3     27.10    590
    15    20  0.96900 -1.32          0.982    62.8     17.70    529
    16    22  0.96700 -1.30          0.983    50.8     11.90    478
    17    24  0.96400 -1.29          0.982    42.0      8.09    435
    18    26  0.96600 -1.28          0.985    35.2      5.63    399
    19    28  0.96800 -1.26          0.987    29.9      4.04    367
    20    30  0.96700 -1.25          0.988    25.7      2.91    340



<strong>png:</strong> 2



```R
#histogram of k
softPower <- 3
adjacency = adjacency(exp_dat, power = softPower);
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average");
k <- softConnectivity(datE=exp_dat,power=softPower) 
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
```

    ..connectivity..
    ..matrix multiplication (system BLAS)..
    ..normalization..
    problem: 32844930 TOM entries are larger than 1.
    ..done.
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 



<table class="dataframe">
<caption>A data.frame: 1 × 2</caption>
<thead>
	<tr><th scope=col>scaleFreeRsquared</th><th scope=col>slope</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>0.86</td><td>11.4</td></tr>
</tbody>
</table>




```R
pdf(file="histogram_of_k.pdf",width=8,height=4.5)
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()
```


<table class="dataframe">
<caption>A data.frame: 1 × 2</caption>
<thead>
	<tr><th scope=col>scaleFreeRsquared</th><th scope=col>slope</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>0.86</td><td>11.4</td></tr>
</tbody>
</table>




<strong>png:</strong> 2



```R
networkType
```


'signed'



```R
saveRDS(exp_dat,file = 'WGCNA_Input.Rds')
```


```R
#SYSTEM ERROR, RUN ON RStudio,one step-wgcna
#net <- blockwiseModules(exp_dat,power=softPower,TOMType = "signed", corType = "pearson",
#                        networkType= "signed",numericLabels=TRUE,nThreads = 4,
#                        mergeCutHeight=0.25,minModuleSize=30, 
#                        maxBlockSize=30000,saveTOMs=TRUE,
#                        saveTOMFileBase="WGCNA_TOM",verbose=5)
```


```R
####WGCNA MODULE
net <- readRDS('./Net.Rds')
load('blockwiseTOM-block-block.1.RData')
```


```R
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
```


```R
unique(net$colors)
table(net$colors)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'turquoise'</li><li>'yellow'</li><li>'brown'</li><li>'green'</li><li>'black'</li><li>'blue'</li><li>'red'</li><li>'grey'</li><li>'magenta'</li><li>'purple'</li><li>'pink'</li></ol>




    
        black      blue     brown     green      grey   magenta      pink    purple 
          669      2310      2241      1336      2603       377       586       133 
          red turquoise    yellow 
         1048      2543      2025 



```R
pdf(file="03.auto_modules_color.pdf",width=12,height=4)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels=FALSE,addGuide=TRUE)
dev.off()
```


<strong>png:</strong> 2



```R
#(2) 134个样本对于18个模块的特征值
dim(net$MEs)
net$MEs[1:4,1:4]
net$MEsOK
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>15</li><li>11</li></ol>




<table class="dataframe">
<caption>A data.frame: 4 × 4</caption>
<thead>
	<tr><th></th><th scope=col>MEturquoise</th><th scope=col>MEyellow</th><th scope=col>MEpurple</th><th scope=col>MEmagenta</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>0.3479056</td><td>0.06483648</td><td>-0.42962929</td><td>-0.16612619</td></tr>
	<tr><th scope=row>H_C_2</th><td>0.3291747</td><td>0.11738321</td><td>-0.05140988</td><td>-0.08359337</td></tr>
	<tr><th scope=row>H_C_3</th><td>0.3057041</td><td>0.08705117</td><td>-0.10647700</td><td>-0.17777469</td></tr>
	<tr><th scope=row>H_C_4</th><td>0.3302650</td><td>0.08013889</td><td>-0.07690690</td><td>-0.18388134</td></tr>
</tbody>
</table>




TRUE



```R
MEs0 <- moduleEigengenes(exp_dat[,net$goodGenes],
                         moduleColors[net$blockGenes[[1]]])$eigengenes
MEs <- orderMEs(MEs0)
rownames(MEs) <- rownames(exp_dat[,net$goodGenes])
text <- cbind(rownames(MEs),MEs)
colnames(text)[1] <- "samples"
write.table(text,file="05.module_eigengenes.xls",
            quote=F,sep="\t",row.names=F)
names(MEs) <- substring(names(MEs),3)
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss),method="average")
```


```R
pdf(file="05.modules_cluster_tree.pdf",width=7,height=5)
plot(METree,main="Clustering of module eigengenes",xlab="",sub="")
dev.off()
```


<strong>png:</strong> 2



```R
moduleCor <- corAndPvalue(MEs,use="p")
rowLabels <- paste("ME",names(MEs),sep="")
textMatrix <- paste(signif(moduleCor$cor,2),
                    "\n(",signif(moduleCor$p,1),")",sep="")
dim(textMatrix) <- dim(moduleCor$cor)
pdf(file="05.modules_relationships.pdf",12,8)
par(mar=c(10,10,1,2))
labeledHeatmap(Matrix=moduleCor$cor,
               textMatrix=textMatrix,
               xLabels=rowLabels,yLabels=rowLabels,
               xSymbols=names(MEs),ySymbols=names(MEs),
               colorLabels=TRUE,colors=blueWhiteRed(50),
               setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()

```


<strong>png:</strong> 2



```R
text <- paste("cor=",round(moduleCor$cor,4),
              ";p-value=",round(moduleCor$p,4),sep="")
dim(text) <- dim(moduleCor$cor)
rownames(text) <- rowLabels
colnames(text) <- rowLabels
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="05.modules_relationships.xls",
            quote=F,sep="\t",row.names=F)
```


```R
dir.create("05.expression_ME")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  dir <- "05.expression_ME/"
  pdf(file=paste(dir,"05.expression_ME_",which.module,".pdf",sep=""),
      25,10)
  ME <- MEs[,which.module]
  ME <- t(as.matrix(MEs[,which.module]))
  colnames(ME) <- rownames(exp_dat[,net$goodGenes])
  layout(matrix(c(1,2)),heights=c(1.5,3))
  par(mar=c(0.3,9,3,5))
  plotMat(t(scale(exp_dat[,net$goodGenes][,moduleColors[net$blockGenes[[1]]]==which.module])),
          nrgcols=30,rlabels=F,rcols=which.module,
          main=paste(which.module),cex.main=1)
  par(mar=c(5,4,0,1))
  barplot(ME,col=which.module,main="",cex.names=1,cex.axis=1,
          ylab="module eigengene",las=3)
  dev.off()
}
```

    Warning message in dir.create("05.expression_ME"):
    “'05.expression_ME' already exists”



```R
MEs
```


<table class="dataframe">
<caption>A data.frame: 15 × 11</caption>
<thead>
	<tr><th></th><th scope=col>greenyellow</th><th scope=col>purple</th><th scope=col>pink</th><th scope=col>magenta</th><th scope=col>red</th><th scope=col>black</th><th scope=col>turquoise</th><th scope=col>yellow</th><th scope=col>green</th><th scope=col>blue</th><th scope=col>brown</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td> 0.06483648</td><td> 0.3479056</td><td>-0.42962929</td><td>-0.3006342</td><td>-0.16612619</td><td>-0.166691429</td><td> 0.37217886</td><td> 0.2473024</td><td>-0.35536647</td><td>-0.33439396</td><td>-0.1362663</td></tr>
	<tr><th scope=row>H_C_2</th><td> 0.11738321</td><td> 0.3291747</td><td>-0.05140988</td><td>-0.2598914</td><td>-0.08359337</td><td>-0.187085595</td><td> 0.31017143</td><td> 0.1460524</td><td> 0.43691470</td><td>-0.33282755</td><td>-0.1837052</td></tr>
	<tr><th scope=row>H_C_3</th><td> 0.08705117</td><td> 0.3057041</td><td>-0.10647700</td><td>-0.3271898</td><td>-0.17777469</td><td> 0.748707277</td><td> 0.32178586</td><td> 0.1495088</td><td> 0.04792924</td><td>-0.34544858</td><td>-0.1390583</td></tr>
	<tr><th scope=row>H_C_4</th><td> 0.08013889</td><td> 0.3302650</td><td>-0.07690690</td><td>-0.3006632</td><td>-0.18388134</td><td>-0.082296854</td><td> 0.35814598</td><td> 0.2390226</td><td> 0.43194359</td><td>-0.36639903</td><td>-0.1591204</td></tr>
	<tr><th scope=row>H_O_1</th><td>-0.32375795</td><td>-0.1728430</td><td>-0.07725040</td><td>-0.1575922</td><td>-0.31176078</td><td> 0.106549938</td><td> 0.09684268</td><td> 0.2796818</td><td>-0.04663179</td><td> 0.13726849</td><td> 0.3059360</td></tr>
	<tr><th scope=row>H_O_2</th><td>-0.28556597</td><td>-0.1953911</td><td> 0.07180182</td><td>-0.2286414</td><td>-0.31812605</td><td> 0.146008639</td><td> 0.12127753</td><td> 0.2754940</td><td>-0.21768303</td><td> 0.08105576</td><td> 0.3312856</td></tr>
	<tr><th scope=row>H_O_3</th><td>-0.29375465</td><td>-0.1713330</td><td>-0.33123680</td><td>-0.1189650</td><td>-0.25225192</td><td>-0.152315643</td><td> 0.08421899</td><td> 0.2541347</td><td>-0.23095087</td><td> 0.15896513</td><td> 0.3054608</td></tr>
	<tr><th scope=row>H_O_4</th><td>-0.33524897</td><td>-0.1821726</td><td> 0.27850605</td><td>-0.1624323</td><td>-0.30173478</td><td> 0.366939851</td><td> 0.03108302</td><td> 0.2468759</td><td>-0.10249877</td><td> 0.11566289</td><td> 0.3090116</td></tr>
	<tr><th scope=row>L_C_1</th><td> 0.39220735</td><td> 0.2206508</td><td>-0.06211244</td><td> 0.1350396</td><td> 0.32829828</td><td> 0.009161127</td><td>-0.02914290</td><td>-0.3163401</td><td>-0.30142982</td><td>-0.17215541</td><td>-0.3712201</td></tr>
	<tr><th scope=row>L_C_2</th><td> 0.44991562</td><td> 0.1923974</td><td> 0.70858572</td><td> 0.2068644</td><td> 0.35663888</td><td> 0.045930998</td><td>-0.13909755</td><td>-0.3753166</td><td>-0.16432159</td><td>-0.13864362</td><td>-0.3829476</td></tr>
	<tr><th scope=row>L_C_3</th><td> 0.43087720</td><td> 0.1447375</td><td> 0.24565880</td><td> 0.2831441</td><td> 0.44764239</td><td>-0.053428272</td><td>-0.21063130</td><td>-0.4264354</td><td>-0.09351377</td><td>-0.07473560</td><td>-0.3977671</td></tr>
	<tr><th scope=row>L_O_1</th><td>-0.04109841</td><td>-0.3012707</td><td>-0.04510761</td><td> 0.3610419</td><td> 0.18186829</td><td>-0.318271212</td><td>-0.36486752</td><td>-0.2392347</td><td> 0.07197072</td><td> 0.34471080</td><td> 0.1096165</td></tr>
	<tr><th scope=row>L_O_2</th><td>-0.11154906</td><td>-0.2710016</td><td> 0.08076816</td><td> 0.2615140</td><td> 0.15224979</td><td>-0.179699080</td><td>-0.31042009</td><td>-0.1780448</td><td>-0.14245337</td><td> 0.30769576</td><td> 0.1456241</td></tr>
	<tr><th scope=row>L_O_3</th><td>-0.12429937</td><td>-0.2769946</td><td>-0.07802412</td><td> 0.2729275</td><td> 0.13838859</td><td>-0.158123928</td><td>-0.29362230</td><td>-0.1307765</td><td> 0.27549061</td><td> 0.30315123</td><td> 0.1451345</td></tr>
	<tr><th scope=row>L_O_4</th><td>-0.10713554</td><td>-0.2998286</td><td>-0.12716612</td><td> 0.3354780</td><td> 0.19016291</td><td>-0.125385818</td><td>-0.34792267</td><td>-0.1719244</td><td> 0.39060061</td><td> 0.31609367</td><td> 0.1180161</td></tr>
</tbody>
</table>




```R
# 06. Relating modules to external information (samples/traits)
sample_cor <- cor(t(exp_dat[,net$goodGenes]),t(exp_dat[,net$goodGenes]),
                  use='pairwise.complete.obs')
moduleSampleCor <- cor(MEs,sample_cor,use="p")
nSamples <- nrow(exp_dat[,net$goodGenes])
moduleSamplePvalue <- corPvalueStudent(moduleSampleCor,nSamples)
textMatrix <- paste(signif(moduleSampleCor,2),
                    "\n(",signif(moduleSamplePvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleSampleCor)
rowLabels <- paste("ME",names(MEs),sep="")
```


```R
pdf(file="06.modules_samples_relationships.pdf",22,6)
par(mar=c(5,12,1,1))
labeledHeatmap(Matrix=moduleSampleCor,
               xLabels=colnames(sample_cor),
               yLabels=rowLabels,ySymbols=names(MEs),
               colorLabels=TRUE,colors=blueWhiteRed(50),
               setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()

```


<strong>png:</strong> 2



```R
text <- paste("cor=",round(moduleSampleCor,4),
              ";p-value=",round(moduleSamplePvalue,4),sep="")
dim(text) <- dim(moduleSampleCor)
rownames(text) <- rownames(moduleSampleCor)
colnames(text) <- colnames(moduleSampleCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="06.modules_samples_relationships.xls",
            quote=F,sep="\t",row.names=F)
```


```R
moduleTraitCor <- cor(MEs,trait_dat,use="p")
nSamples <- nrow(exp_dat[,net$goodGenes])
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste(signif(moduleTraitCor,2),
                    "\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix) <- dim(moduleTraitCor)
rowLabels <- paste("ME",names(MEs),sep="")
```


```R
pdf(file="06.modules_traits_relationships.pdf",15,8)
par(mar=c(7,12,1,2))
labeledHeatmap(Matrix=moduleTraitCor,
               textMatrix=textMatrix,
               xLabels=colnames(trait_dat),
               yLabels=rowLabels,ySymbols=names(MEs),
               colorLabels=TRUE,colors=blueWhiteRed(50),
               setStdMargins=FALSE,xLabelsAngle=90,zlim=c(-1,1))
dev.off()
```


<strong>png:</strong> 2



```R
text <- paste("cor=",round(moduleTraitCor,4),
              ";p-value=",round(moduleTraitPvalue,4),sep="")
dim(text) <- dim(moduleTraitCor)
rownames(text) <- rownames(moduleTraitCor)
colnames(text) <- colnames(moduleTraitCor)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="06.modules_traits_relationships.xls",
            quote=F,sep="\t",row.names=F)
```


```R
# 07. Relating nodes to external information (samples/traits)
modNames <- names(MEs)
geneModuleMembership <- cor(exp_dat[,net$goodGenes],MEs,use="p")
nSamples <- nrow(exp_dat[,net$goodGenes])
MMPvalue <- corPvalueStudent(geneModuleMembership,nSamples)
colnames(geneModuleMembership) <- paste("MM",modNames,sep="")
colnames(MMPvalue) <- paste("p.MM",modNames,sep="")
text <- paste("cor=",round(geneModuleMembership,4),
              ";p-value=",round(MMPvalue,4),sep="")
dim(text) <- dim(geneModuleMembership)
rownames(text) <- rownames(geneModuleMembership)
colnames(text) <- colnames(geneModuleMembership)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "modules"
write.table(text,file="07.genes_module_membership.xls",
            quote=F,sep="\t",row.names=F)
```


```R
trait_dat
```


<table class="dataframe">
<caption>A data.frame: 15 × 2</caption>
<thead>
	<tr><th></th><th scope=col>high_glucose</th><th scope=col>oxldl</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_2</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_3</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_C_4</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>H_O_1</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>H_O_2</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>H_O_3</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>H_O_4</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>L_C_1</th><td>0</td><td>0</td></tr>
	<tr><th scope=row>L_C_2</th><td>0</td><td>0</td></tr>
	<tr><th scope=row>L_C_3</th><td>0</td><td>0</td></tr>
	<tr><th scope=row>L_O_1</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>L_O_2</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>L_O_3</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>L_O_4</th><td>0</td><td>1</td></tr>
</tbody>
</table>




```R
net$MEs %>%
rownames_to_column(var ='Samples') %>%
write_csv('Samples_ME.csv')
```


```R
trait_dat$high_glucose
trait_dat$oxldl
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>1</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li><li>0</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0</li><li>0</li><li>0</li><li>0</li><li>1</li><li>1</li><li>1</li><li>1</li><li>0</li><li>0</li><li>0</li><li>1</li><li>1</li><li>1</li><li>1</li></ol>




```R
as.data.frame(trait_dat$high_glucose,trait_dat$oxldl)
```

    Warning message in as.data.frame.numeric(trait_dat$high_glucose, trait_dat$oxldl):
    “'row.names' is not a character vector of length 15 -- omitting it. Will be an error!”



<table class="dataframe">
<caption>A data.frame: 15 × 1</caption>
<thead>
	<tr><th scope=col>trait_dat$high_glucose</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1</td></tr>
	<tr><td>1</td></tr>
	<tr><td>1</td></tr>
	<tr><td>1</td></tr>
	<tr><td>1</td></tr>
	<tr><td>1</td></tr>
	<tr><td>1</td></tr>
	<tr><td>1</td></tr>
	<tr><td>0</td></tr>
	<tr><td>0</td></tr>
	<tr><td>0</td></tr>
	<tr><td>0</td></tr>
	<tr><td>0</td></tr>
	<tr><td>0</td></tr>
	<tr><td>0</td></tr>
</tbody>
</table>




```R
modNames <- names(MEs)
trigly <- data.frame(glucose = trait_dat$high_glucose,oxldl = trait_dat$oxldl)
```


```R
modNames <- names(MEs)
trigly <- data.frame(glucose = trait_dat$high_glucose,oxldl = trait_dat$oxldl)
geneTraitSignificance <- as.data.frame(cor(exp_dat[,net$goodGenes],trigly,use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste("GS.",names(trigly),sep="")
names(GSPvalue) <- paste("p.GS.",names(trigly),sep="")
head(geneTraitSignificance)
head(GSPvalue)
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>GS.glucose</th><th scope=col>GS.oxldl</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Gpnmb</th><td> 0.25451282</td><td>-0.9707996</td></tr>
	<tr><th scope=row>Psap</th><td>-0.52419170</td><td>-0.7603028</td></tr>
	<tr><th scope=row>Fth1</th><td> 0.05826131</td><td> 0.8905958</td></tr>
	<tr><th scope=row>Ctsd</th><td>-0.05024026</td><td>-0.8884569</td></tr>
	<tr><th scope=row>Ctsb</th><td> 0.20119783</td><td>-0.9784023</td></tr>
	<tr><th scope=row>mt-Co1</th><td> 0.42323048</td><td> 0.3757233</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>p.GS.glucose</th><th scope=col>p.GS.oxldl</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Gpnmb</th><td>0.35997540</td><td>1.940657e-09</td></tr>
	<tr><th scope=row>Psap</th><td>0.04487236</td><td>1.001176e-03</td></tr>
	<tr><th scope=row>Fth1</th><td>0.83660087</td><td>8.536343e-06</td></tr>
	<tr><th scope=row>Ctsd</th><td>0.85887222</td><td>9.629722e-06</td></tr>
	<tr><th scope=row>Ctsb</th><td>0.47211711</td><td>2.783012e-10</td></tr>
	<tr><th scope=row>mt-Co1</th><td>0.11597338</td><td>1.675473e-01</td></tr>
</tbody>
</table>




```R
dir.create("07.MM_vs_trigly")
for(i in 1:(ncol(MEs)-1)) {
  which.module <- labels2colors(i)
  column <- match(which.module,modNames)
  moduleGenes <- moduleColors[net$blockGenes[[1]]]==which.module
  dir <- "07.MM_vs_trigly/"
  pdf(file=paste(dir,"07.",which.module,"_MM_vs_trigly.pdf",sep=""),6,6)
  verboseScatterplot(geneModuleMembership[moduleGenes,column],
                     geneTraitSignificance[moduleGenes,1],
                     xlab=paste("Module membership (MM) in",which.module,"module"),
                     ylab="Gene significance for trigly",
                     main=paste("Module membership vs. gene significance\n"),
                     col=which.module)
  dev.off()
}

```

    Warning message in dir.create("07.MM_vs_trigly"):
    “'07.MM_vs_trigly' already exists”



```R
text <- cbind(geneTraitSignificance,GSPvalue)
text <- cbind(rownames(text),text)
colnames(text)[1] <- "genes"
write.table(text,file="07.genes_trait_significance.xls",
            quote=F,sep="\t",row.names=F)
```


```R
# 08. Exporting network
ATOM <- as.matrix(TOM)
TOM1 <- ATOM[1:round((nrow(ATOM)/2)),1:round((nrow(ATOM)/2))]
TOM2 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),1:round((nrow(ATOM)/2))]
TOM3 <- ATOM[1:round((nrow(ATOM)/2)),(round(nrow(ATOM)/2)+1):nrow(ATOM)]
TOM4 <- ATOM[(round(nrow(ATOM)/2)+1):nrow(ATOM),(round(nrow(ATOM)/2)+1):nrow(ATOM)]
dir.create("08.module_result")
setwd("08.module_result")
```

    Warning message in dir.create("08.module_result"):
    “'08.module_result' already exists”



```R
for(i in 1:(ncol(MEs)-1)) {
  module <- labels2colors(i)
  inModule <- moduleColors[net$blockGenes[[1]]]==module
  genename <- colnames(exp_dat[,net$goodGenes])
  modGenes <- genename[inModule]
  modTOM1 <- TOM1[inModule[1:round((nrow(ATOM)/2))],
                  inModule[1:round((nrow(ATOM)/2))]]
  modTOM2 <- TOM2[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)],
                  inModule[1:round((nrow(ATOM)/2))]]
  modTOM3 <- TOM3[inModule[1:round((nrow(ATOM)/2))],
                  inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]
  modTOM4 <- TOM4[inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)],
                  inModule[(round(nrow(ATOM)/2)+1):nrow(ATOM)]]
  modTOM <- rbind(cbind(modTOM1,modTOM3),cbind(modTOM2,modTOM4))
  IMConn <- softConnectivity(exp_dat[,net$goodGenes][,modGenes])
  
  cyt1 <- exportNetworkToCytoscape(modTOM,
                                   edgeFile=paste("CytoscapeInput-edges-",paste(module,collapse="-"),".txt",sep=""),
                                   nodeFile=paste("CytoscapeInput-nodes-",paste(module,collapse="-"),".txt",sep=""),
                                   weighted=TRUE,threshold=0.02,
                                   nodeNames=modGenes,altNodeNames=modGenes,
                                   nodeAttr=moduleColors[net$blockGenes[[1]]][inModule])
  
  out <- cbind(modGenes,IMConn)
  colnames(out) <- c("gene","connectivity")
  out <- out[order(as.numeric(out[,2]),decreasing=T),]
  write.table(out,paste(module,"-module-gene.txt",sep=""),
              sep="\t",quote=F,row.names=F)
  
  nTop <- 0.05*length(modGenes)
  top <- (rank(-IMConn) <= nTop)
  out <- cbind(modGenes[top],IMConn[top])
  colnames(out) <- c("gene","connectivity")
  out <- out[order(as.numeric(out[,2]),decreasing=T),]
  write.table(out,paste(module,"-5%hubgene.txt",sep=""),
              sep="\t",quote=F,row.names=F)
  
  nTop <- 0.1*length(modGenes)
  top <- (rank(-IMConn) <= nTop)
  out <- cbind(modGenes[top],IMConn[top])
  colnames(out) <- c("gene","connectivity")
  out <- out[order(as.numeric(out[,2]),decreasing=T),]
  write.table(out,paste(module,"-10%hubgene.txt",sep=""),
              sep="\t",quote=F,row.names=F)
  
  nTop <- 25
  top <- (rank(-IMConn) <= nTop)
  out <- cbind(modGenes[top],IMConn[top])
  colnames(out) <- c("gene","connectivity")
  out <- out[order(as.numeric(out[,2]),decreasing=T),]
  write.table(out,paste(module,"-top25_hubgene.txt",sep=""),
              sep="\t",quote=F,row.names=F)
}
```

     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 



```R
genename <- colnames(exp_dat[,net$goodGenes])
modGenes <- genename
IMConn <- softConnectivity(exp_dat[,net$goodGenes][,modGenes])
cyt1 <- exportNetworkToCytoscape(ATOM,
                                 edgeFile=paste("CytoscapeInput-edges-overall.txt",sep=""),
                                 nodeFile=paste("CytoscapeInput-nodes-overall.txt",sep=""),
                                 weighted=TRUE,threshold=0.02,
                                 nodeNames=modGenes,altNodeNames=modGenes,
                                 nodeAttr=moduleColors[net$blockGenes[[1]]])
out <- cbind(modGenes,IMConn)
colnames(out) <- c("gene","connectivity")
out <- out[order(as.numeric(out[,2]),decreasing=T),]
write.table(out,paste("overall_gene_connectivity.txt",sep=""),
            sep="\t",quote=F,row.names=F)
IntraMConn <- intramodularConnectivity.fromExpr(exp_dat[,net$goodGenes],
                                                moduleColors[net$blockGenes[[1]]],
                                                networkType=networkType,
                                                power=softPower)
merged_IntraMConn <- cbind(colnames(exp_dat[,net$goodGenes]),
                           moduleColors[net$blockGenes[[1]]],IntraMConn)
colnames(merged_IntraMConn)[1:2] <- c("gene","module")
head(merged_IntraMConn)
write.table(merged_IntraMConn,file="intraModular_connectivity.xls",
            quote=F,sep="\t",row.names=FALSE)
```

     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 
     softConnectivity: FYI: connecitivty of genes with less than 5 valid samples will be returned as NA.
     ..calculating connectivities.. 



<table class="dataframe">
<caption>A data.frame: 6 × 6</caption>
<thead>
	<tr><th></th><th scope=col>gene</th><th scope=col>module</th><th scope=col>kTotal</th><th scope=col>kWithin</th><th scope=col>kOut</th><th scope=col>kDiff</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>Gpnmb </td><td>purple     </td><td>7291.728</td><td>2542.000</td><td>4749.728</td><td>-2207.728</td></tr>
	<tr><th scope=row>2</th><td>Psap  </td><td>greenyellow</td><td>6925.774</td><td>2021.394</td><td>4904.380</td><td>-2882.986</td></tr>
	<tr><th scope=row>3</th><td>Fth1  </td><td>brown      </td><td>7159.509</td><td>2226.606</td><td>4932.903</td><td>-2706.297</td></tr>
	<tr><th scope=row>4</th><td>Ctsd  </td><td>purple     </td><td>6974.573</td><td>2517.519</td><td>4457.054</td><td>-1939.536</td></tr>
	<tr><th scope=row>5</th><td>Ctsb  </td><td>purple     </td><td>7293.711</td><td>2541.922</td><td>4751.789</td><td>-2209.868</td></tr>
	<tr><th scope=row>6</th><td>mt-Co1</td><td>brown      </td><td>6986.586</td><td>1966.883</td><td>5019.703</td><td>-3052.820</td></tr>
</tbody>
</table>


