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


```R
immune.combined <- readRDS('~/AS/Step3_AddDoublet_Seurat.Rds')
```


```R
#cell cycle score
load("~/钱教授数据分析v2.0/human_cycle.rda")
DefaultAssay(immune.combined) <- 'RNA'
immune.combined <- CellCycleScoring(immune.combined,g2m.features = g2m_genes,s.features = s_genes)
```


```R
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
```


```R
#stress score
stree.genes <- read_csv('~/Resources/stress.genes.csv') %>% pull (gene)
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(stree.genes),
                 ctrl = length(stree.genes),
                 name = 'Stress.Score',
                 assay = "integrated")
```

    [1mRows: [22m[34m32[39m [1mColumns: [22m[34m1[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m (1): gene
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
    Warning message:
    “The following features are not present in the object: FOSB, JUNB, UBC, ZFP36, HSPA8, MT1, IER2, DNAJA1, JUND, PPP1R15A, HSPE1, DUSP1, HSP90AB1, not searching for symbol synonyms”



```R
###ifn genes 
library(msigdbr)
mdb <- msigdbr(species = "Homo sapiens", category = "C2")
IFN_pathway <- mdb[str_detect(mdb$gs_name,str_to_upper("BROWNE_INTERFERON_RESPONSIVE_GENES")),]
geneSets <- lapply(unique(IFN_pathway$gs_name), function(x){print(x);IFN_pathway$gene_symbol[IFN_pathway$gs_name == x]})
names(geneSets) <- "BROWNE_INTERFERON_RESPONSIVE_GENES"
ifn_gene <- geneSets$BROWNE_INTERFERON_RESPONSIVE_GENES[geneSets$BROWNE_INTERFERON_RESPONSIVE_GENES %in% rownames(immune.combined)]
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(ifn_gene),
                 ctrl = length(ifn_gene),
                 name = 'IFN.Score',
                 assay = "RNA")
```

    [1] "BROWNE_INTERFERON_RESPONSIVE_GENES"



```R
features <- VariableFeatures(immune.combined)
length(features)
```


3000



```R
#Find cell cycle related features
source("~/钱教授数据分析v2.0/tools/scTools.R")
cycle_related <- extract_cellcycle(immune.combined, features, cores = 1, cutoff = 0.09, assay = "integrated")
length(cycle_related)
cycle_related_position <- match(cycle_related, features)
features_subtr_1 <- features[-c(cycle_related_position)]
```


322



```R
#scaleData
immune.combined<- ScaleData(immune.combined, verbose = FALSE)
```


```R
immune.combined <- RunPCA(immune.combined, npcs = 50, do.print = TRUE, verbose = FALSE, pcs.print = 1:5,genes.print = 5,set.seed=666, features =  features)
```


```R
VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca")
DimPlot(immune.combined, reduction = "pca")
DimHeatmap(immune.combined, dims = 1:30, cells = 500, balanced = TRUE)
#immune.combined <- JackStraw(immune.combined,dims=50,num.replicate = 50)
#immune.combined <- ScoreJackStraw(immune.combined, dims = 1:50)
#JackStrawPlot(immune.combined, dims = 1:50)
ElbowPlot(immune.combined,ndims=50)
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_11_0.png)
    



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_11_1.png)
    



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_11_2.png)
    



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_11_3.png)
    



```R
#PCA SET 25
immune.combined <- FindNeighbors(immune.combined, dims = 1:25)
res.use <- c(seq(0.1,3,by=0.2),3,5)
res.use
  immune.combined <- FindClusters(object =immune.combined,resolution = res.use)
clust.tree.out <- clustree(immune.combined)+
  theme(legend.position = "bottom")+
  scale_color_brewer(palette = "Set1")
clust.tree.out
```

    Computing nearest neighbor graph
    
    Computing SNN
    



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0.1</li><li>0.3</li><li>0.5</li><li>0.7</li><li>0.9</li><li>1.1</li><li>1.3</li><li>1.5</li><li>1.7</li><li>1.9</li><li>2.1</li><li>2.3</li><li>2.5</li><li>2.7</li><li>2.9</li><li>3</li><li>5</li></ol>



    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9757
    Number of communities: 9
    Elapsed time: 11 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9470
    Number of communities: 18
    Elapsed time: 10 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9298
    Number of communities: 22
    Elapsed time: 10 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9175
    Number of communities: 24
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9054
    Number of communities: 27
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8945
    Number of communities: 29
    Elapsed time: 10 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8846
    Number of communities: 34
    Elapsed time: 7 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8756
    Number of communities: 39
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8684
    Number of communities: 44
    Elapsed time: 9 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8619
    Number of communities: 45
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8551
    Number of communities: 48
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8490
    Number of communities: 52
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8442
    Number of communities: 52
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8393
    Number of communities: 58
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8342
    Number of communities: 56
    Elapsed time: 8 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8318
    Number of communities: 58
    Elapsed time: 9 seconds
    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 48624
    Number of edges: 1823291
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7937
    Number of communities: 80
    Elapsed time: 8 seconds


    Warning message:
    “The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
    Please use the `.add` argument instead.
    [90mThis warning is displayed once every 8 hours.[39m
    [90mCall `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.[39m”
    Warning message in RColorBrewer::brewer.pal(n, pal):
    “n too large, allowed maximum for palette Set1 is 9
    Returning the palette you asked for with that many colors
    ”
    Warning message:
    “Removed 449 rows containing missing values (geom_point).”



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_12_4.png)
    



```R
immune.combined <- RunUMAP(immune.combined,reduction.key = "UMAP_", seed.use=666,dims  = 1:25)
```

    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    17:50:25 UMAP embedding parameters a = 0.9922 b = 1.112
    
    17:50:25 Read 48624 rows and found 25 numeric columns
    
    17:50:25 Using Annoy for neighbor search, n_neighbors = 30
    
    17:50:25 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    |
    
    17:50:28 Writing NN index file to temp file /tmp/RtmpM30DIS/file8a4601b1b46
    
    17:50:28 Searching Annoy index using 1 thread, search_k = 3000
    
    17:50:40 Annoy recall = 100%
    
    17:50:40 Commencing smooth kNN distance calibration using 1 thread
    
    17:50:42 Initializing from normalized Laplacian + noise
    
    17:50:46 Commencing optimization for 200 epochs, with 2132656 positive edges
    
    17:51:00 Optimization finished
    



```R
immune.combined@meta.data$Stress.Score1
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>-0.0666544878311134</li><li>-0.0738759474219872</li><li>0.048382518859524</li><li>-0.026150277504219</li><li>-0.0335963814859453</li><li>-0.0251245470933316</li><li>0.0686942338186471</li><li>-0.0323065403823726</li><li>-0.00443855487678985</li><li>-0.0405693652983748</li><li>0.025399709460054</li><li>-0.0110218854755907</li><li>-0.0367808544553066</li><li>-0.0490148967041726</li><li>-0.0801974252256246</li><li>-0.064427720598147</li><li>-0.00181403346021931</li><li>0.101076644820056</li><li>-0.0554329795994628</li><li>-0.0854497477354915</li><li>-0.0504671393081662</li><li>0.00355652547421549</li><li>-0.0292693844326188</li><li>0.0389814735351569</li><li>-0.0382552726531642</li><li>-0.0487506436074566</li><li>-0.00441827656790839</li><li>-0.078061236926398</li><li>-0.0150115376624363</li><li>-0.0600521361559936</li><li>-0.01239041980364</li><li>-0.118952668594297</li><li>0.0242350924100708</li><li>-0.0086432898252658</li><li>0.0741424738087741</li><li>0.0113013829518421</li><li>0.0664389722885068</li><li>-0.0507079612902981</li><li>-0.0700955421715783</li><li>-0.0533273395901273</li><li>-0.0355823855390491</li><li>-0.0848109749122798</li><li>-0.10029154265142</li><li>-0.014011658664747</li><li>-0.0521666493573535</li><li>-0.0487611376683682</li><li>0.0850558341721352</li><li>-0.0845255599993933</li><li>-0.0632214467079285</li><li>0.0320553122899987</li><li>-0.0645462495538601</li><li>-0.0599163666743508</li><li>-0.00385516810881444</li><li>-0.0316487984873311</li><li>-0.0766091598013287</li><li>-0.0354493673450289</li><li>-0.0507280965131566</li><li>-0.0568389665922452</li><li>0.0674015047992613</li><li>-0.0228917027482111</li><li>0.0100838770567095</li><li>-0.00696371939585896</li><li>-0.0404253909175305</li><li>-0.0289188067907679</li><li>-0.0905193913776353</li><li>-0.00656707295772754</li><li>-0.0873189327783411</li><li>-0.0677998279796289</li><li>-0.0508217818207955</li><li>-0.0283917870549486</li><li>-0.037868791303237</li><li>0.0590572879123816</li><li>-0.0547672938378433</li><li>-0.0703737159400897</li><li>-0.0094185429473326</li><li>-0.0885114865869745</li><li>-0.020130442382622</li><li>-0.0479226765405389</li><li>-0.0411169813353815</li><li>-0.0257227018914267</li><li>-0.0357362155412946</li><li>-0.0415007700047448</li><li>-0.031758787849458</li><li>-0.0809051576243655</li><li>-0.000174499623141214</li><li>-0.0896970714880252</li><li>-0.0956800814877981</li><li>-0.049922575756838</li><li>-0.0511992632676018</li><li>-0.0522915990469706</li><li>-0.0505187406982241</li><li>-0.046637193521948</li><li>-0.0528797949190759</li><li>-0.0182136854404137</li><li>0.157609761891044</li><li>-0.0671114350197743</li><li>-0.00718717997525509</li><li>-0.0216455466137774</li><li>-0.0217433521383519</li><li>0.0361366367828464</li><li>-0.0747716778344549</li><li>-0.127266461744111</li><li>-0.0276837706715387</li><li>-0.0294843763498984</li><li>-0.0109243278211807</li><li>-0.0483884698217913</li><li>0.00113917506739833</li><li>-0.0365485412382201</li><li>-0.0875803418606639</li><li>-0.00132703893339861</li><li>0.0318177502457398</li><li>-0.0470453832183324</li><li>-0.0402219408799275</li><li>-0.0624516625656205</li><li>-0.0122277811761164</li><li>0.0881371472205261</li><li>0.0230141702395705</li><li>-0.024731478323333</li><li>0.292681053263492</li><li>0.0123682291751936</li><li>-0.0566412860494845</li><li>-0.0934109330981504</li><li>0.0252047028775644</li><li>-0.0428688386287718</li><li>-0.0968347205091152</li><li>-0.038345948476124</li><li>-0.0146036424860554</li><li>-0.0793115394176611</li><li>-0.00275700639834697</li><li>-0.0552037431620277</li><li>-0.0150306670307385</li><li>-0.0498395589939198</li><li>-0.0404458511368812</li><li>-0.0255202070119733</li><li>-0.0631091816626018</li><li>0.0574733156996335</li><li>-0.0840955524968244</li><li>-0.0136315952362863</li><li>-0.00145844327053614</li><li>-0.0464988166530143</li><li>0.0010424974105662</li><li>-0.0723136347299311</li><li>-0.0649179818235524</li><li>0.0387491370180653</li><li>0.00470178725868814</li><li>-0.0731354270519458</li><li>-0.0712832661010472</li><li>-0.0359677567424611</li><li>0.00478933672877162</li><li>-0.0697716128161064</li><li>-0.0508189068570787</li><li>-0.0043210912736631</li><li>-0.0102088037089051</li><li>0.0155448497820454</li><li>-0.068922781463891</li><li>-0.0839646455911917</li><li>-0.0882132610004257</li><li>-0.0488666578202333</li><li>-0.0722719312330746</li><li>-0.0404963907587539</li><li>-0.0645670119969883</li><li>0.00591833090353648</li><li>0.0115584019769019</li><li>-0.0450695642477612</li><li>-0.0252512074564838</li><li>-0.0117045771940897</li><li>-0.0680257094475833</li><li>-0.0389043771368703</li><li>-0.0425466475957363</li><li>-0.0524451514072619</li><li>0.00139996434410586</li><li>-0.0709768959004456</li><li>-0.0312683265843474</li><li>-0.0602916476311805</li><li>1.26784398508036</li><li>0.0200470388172101</li><li>0.0938706686239756</li><li>0.00207862655388651</li><li>-0.0220343893096195</li><li>-0.0460473661188392</li><li>-0.100920316489935</li><li>0.00832959719248076</li><li>-0.0635002235791133</li><li>-0.0750248921649951</li><li>-0.0089652647246613</li><li>0.0567389356848807</li><li>-0.0583511767831093</li><li>-0.0400554617177359</li><li>-0.0576167840352331</li><li>-0.0555127787861053</li><li>-0.0130757042449474</li><li>-0.055653764255653</li><li>-0.0778222749566976</li><li>-0.0606055419380724</li><li>0.0341707374780343</li><li>-0.0200976051321338</li><li>-0.0383006933379721</li><li>0.00956405015444861</li><li>0.00770771910033266</li><li>0.010014829904458</li><li>⋯</li><li>-0.0827381398077077</li><li>-0.0492936608581075</li><li>-0.00912362567733421</li><li>-0.0509832538179482</li><li>-0.0204744470423982</li><li>-0.0607833868843636</li><li>-0.00426270136219767</li><li>-0.011980970012206</li><li>-0.094248882402605</li><li>-0.103169366950419</li><li>-0.0485562100296471</li><li>-0.0752119378637581</li><li>-0.058328903697548</li><li>-0.0847265483148911</li><li>-0.0627394594807959</li><li>0.0194293431584241</li><li>-0.0590933058761773</li><li>-0.0453477423628039</li><li>-0.03088452537158</li><li>-0.0536905250845395</li><li>0.0346605351375204</li><li>-0.0781163201465906</li><li>-0.0343681031826622</li><li>-0.0780023208511081</li><li>-0.0943251207372913</li><li>-0.07547066782901</li><li>-0.0540109032794144</li><li>-0.0912886154541005</li><li>-0.0485171170892278</li><li>0.00409960042188934</li><li>-0.0227752743908602</li><li>-0.0110042277108551</li><li>-0.0439549982619603</li><li>-0.0644624070045592</li><li>-0.105802096342425</li><li>0.0271352577380996</li><li>-0.0496762981875047</li><li>-0.0311090049882159</li><li>0.0164267482484459</li><li>0.00386297850777684</li><li>-0.00318060988094271</li><li>0.0611175262294193</li><li>-0.0883495133352492</li><li>0.011439577139783</li><li>-0.110555482784879</li><li>-0.0590494784507824</li><li>-0.0306310971903084</li><li>-0.114896319975343</li><li>-0.0613652971658448</li><li>-0.0165007137365437</li><li>-0.029102398741561</li><li>0.0621589799885186</li><li>-0.0549719884727759</li><li>-0.0887711752713314</li><li>0.0215519222857806</li><li>0.0297054988200376</li><li>-0.0877263244220818</li><li>-0.0167138773143734</li><li>-0.0312340928970071</li><li>0.0067462193768997</li><li>0.00036404606612464</li><li>-0.0141684706549819</li><li>-0.0388643618826261</li><li>-0.0663636173859404</li><li>-0.0748274054648072</li><li>-0.0198747396015279</li><li>0.125641120188502</li><li>-0.0783798539913566</li><li>-0.0180631002207099</li><li>0.0173096232097381</li><li>-0.0272716244145186</li><li>-0.00540098693563076</li><li>-0.0452425645423871</li><li>-0.12153844395165</li><li>-0.125780500578325</li><li>-0.0388880596111086</li><li>-0.0835938381428192</li><li>-0.0166931365171339</li><li>-0.0977817198901355</li><li>-0.0821314132980847</li><li>0.0314439517139982</li><li>-0.109379200731962</li><li>-0.0735318111419373</li><li>-0.0373866958043579</li><li>-0.0585123912796256</li><li>-0.0173730625916449</li><li>-0.0549520439068674</li><li>-0.0916790203237346</li><li>-0.058067567759879</li><li>-0.0606950291877954</li><li>-0.0909099539552469</li><li>0.0222408001085521</li><li>-0.0225587474416229</li><li>-0.0293132019293449</li><li>-0.0359418110812979</li><li>-0.0735726408281381</li><li>-0.0283421208564512</li><li>0.00530168027629557</li><li>-0.129000076610006</li><li>-0.0799703895875452</li><li>0.00891024907635141</li><li>-0.0855078859958481</li><li>-0.0489900978077556</li><li>-0.0239491773167113</li><li>-0.0282248450061204</li><li>-0.0513024926833721</li><li>-0.0414122557120904</li><li>-0.0527386298263484</li><li>-0.0248925272417032</li><li>-0.0384727326371573</li><li>-0.072240799085061</li><li>-0.0402042578777849</li><li>-0.0204516346441819</li><li>-0.000452966590118398</li><li>-0.0606514205168711</li><li>-0.041559568482949</li><li>-0.00256663687692418</li><li>-0.0319851685300218</li><li>-0.0431529668766238</li><li>-0.0371173636411831</li><li>-0.0544515309813299</li><li>-0.0308289388525085</li><li>-0.0487993368863372</li><li>-0.0441804863785352</li><li>-0.00245637682103177</li><li>-0.098226246108477</li><li>-0.022341002163371</li><li>-0.0790982943066555</li><li>-0.0370800456053287</li><li>-0.063413813505307</li><li>-0.0315985033096963</li><li>-0.0671472937830964</li><li>-0.129987343088148</li><li>-0.0158654549518997</li><li>-0.0355682966671307</li><li>-0.0590264268495247</li><li>-0.0547874027200237</li><li>-0.0656475365314802</li><li>0.0219067646336422</li><li>-0.0532070161988431</li><li>0.0377288961379842</li><li>-0.0015900739888374</li><li>-0.0699269724634762</li><li>-0.0307076261890804</li><li>-0.0527760725134383</li><li>-0.0215179696009459</li><li>-0.0324506080920378</li><li>-0.0338011389370909</li><li>-0.0268231782097489</li><li>-0.0336429554066844</li><li>-0.0185704489894381</li><li>-0.0813732836054161</li><li>-0.0396800298638067</li><li>-0.0439268403979788</li><li>0.0216814653859645</li><li>-0.0401544773585033</li><li>-0.00544040822665585</li><li>-0.00446663535268343</li><li>-0.0423077200851032</li><li>-0.0165101442381501</li><li>-0.0791645137465845</li><li>-0.0422448064891901</li><li>-0.0554789484465841</li><li>-0.107299572408655</li><li>-0.0638770319244123</li><li>-0.0579442904925633</li><li>-0.0587703166672314</li><li>-0.0330463369874313</li><li>-0.111373584179021</li><li>-0.100805745325989</li><li>-0.044551064327668</li><li>-0.073755736664203</li><li>0.00607178360273859</li><li>-0.104383672028076</li><li>-0.0223134465434726</li><li>-0.0295119297656342</li><li>-0.0979378590607274</li><li>-0.0681212564670517</li><li>-0.0410871903164212</li><li>-0.0725092264051528</li><li>-0.0312134778084214</li><li>0.00068072772712964</li><li>-0.0427407684906622</li><li>-0.0921229522750205</li><li>-0.0249890322884921</li><li>-0.0673073471621926</li><li>0.0140836552180424</li><li>0.0753273913915434</li><li>-0.0552182151424342</li><li>-0.0841623090837957</li><li>-0.0576696580977904</li><li>-0.022951247796796</li><li>-0.0619624077953798</li><li>-0.0599119792292472</li><li>-0.0921044103723688</li><li>-0.0490832964190956</li><li>-0.0517524079474129</li><li>-0.064051510330649</li><li>-0.0233579489615152</li><li>0.0686303057824926</li></ol>




```R
FeaturePlot(immune.combined, features = 'IFN.Score1')
FeaturePlot(immune.combined, features = 'Stress.Score1')
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_15_0.png)
    



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_15_1.png)
    



```R
FeaturePlot(immune.combined, features = 'G2M.Score')
FeaturePlot(immune.combined, features = 'S.Score')
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_16_0.png)
    



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_16_1.png)
    



```R
DefaultAssay(immune.combined) <- 'RNA'
FeaturePlot(immune.combined, features = 'CD3E')
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_17_0.png)
    



```R
DimPlot(immune.combined, group.by = 'integrated_snn_res.0.1',label = T) +NoLegend()
DimPlot(immune.combined, group.by = 'integrated_snn_res.1.1',label = T) +NoLegend()
DimPlot(immune.combined, group.by = 'DoubletFinder',label = T) +NoLegend()
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_18_0.png)
    



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_18_1.png)
    



    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_18_2.png)
    



```R
#Findmarker
DefaultAssay(immune.combined) <- 'RNA'
Idents(immune.combined) <- 'integrated_snn_res.0.1'
immune.combined.markers <- FindAllMarkers(immune.combined,min.pct = 0.05,logfc.threshold = 0.25, only.posT=T, assay = 'RNA',
                                          min.diff.pct=0.1, only.pos = T, random.seed = 666)

```

    Calculating cluster 0
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Calculating cluster 8
    



    Error in paste(path, "Findmarkers_res_0.1.csv"): object 'path' not found
    Traceback:


    1. write.table(immune.combined.markers, file = paste(path, "Findmarkers_res_0.1.csv"), 
     .     sep = ",", quote = F)

    2. paste(path, "Findmarkers_res_0.1.csv")



```R
write.table(immune.combined.markers,file='Findmarkers_res_0.1.csv',sep=',',quote=F)
```


```R
############doublet finder
tibble(nFeature = immune.combined@meta.data$nFeature_RNA, UMI = immune.combined@meta.data$nCount_RNA, cluster = immune.combined@meta.data$`integrated_snn_res.0.1`) %>%
group_by(cluster) %>%
summarise(nFeature_Mean = mean(nFeature), UMI_mean = mean(UMI)) 
```


<table class="dataframe">
<caption>A tibble: 9 × 3</caption>
<thead>
	<tr><th scope=col>cluster</th><th scope=col>nFeature_Mean</th><th scope=col>UMI_mean</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>0</td><td>1104.880</td><td>3156.674</td></tr>
	<tr><td>1</td><td>2008.048</td><td>6554.905</td></tr>
	<tr><td>2</td><td>2184.103</td><td>8900.105</td></tr>
	<tr><td>3</td><td>2180.614</td><td>7327.183</td></tr>
	<tr><td>4</td><td>1124.093</td><td>3687.182</td></tr>
	<tr><td>5</td><td>1189.853</td><td>3119.824</td></tr>
	<tr><td>6</td><td>1345.636</td><td>4222.729</td></tr>
	<tr><td>7</td><td>2156.234</td><td>6847.303</td></tr>
	<tr><td>8</td><td>1581.244</td><td>8329.566</td></tr>
</tbody>
</table>




```R
tibble(DoubletFinder = immune.combined@meta.data$DoubletFinder, DoubletScore = immune.combined@meta.data$DoubletScore, cluster = immune.combined@meta.data$`integrated_snn_res.0.1`) %>%
group_by(cluster) %>%
mutate(Doublet = ifelse(DoubletFinder == 'Doublet', '1','0')) %>%
mutate(Doublet = as.numeric(Doublet)) %>%
summarise(DoubletScore = mean(DoubletScore), cellcount = n(), Doublet = sum(Doublet))
```


<table class="dataframe">
<caption>A tibble: 9 × 4</caption>
<thead>
	<tr><th scope=col>cluster</th><th scope=col>DoubletScore</th><th scope=col>cellcount</th><th scope=col>Doublet</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>0</td><td>0.1183079</td><td>17286</td><td>149</td></tr>
	<tr><td>1</td><td>0.2488774</td><td>10289</td><td>630</td></tr>
	<tr><td>2</td><td>0.2445581</td><td>10001</td><td>674</td></tr>
	<tr><td>3</td><td>0.2376731</td><td> 6312</td><td>790</td></tr>
	<tr><td>4</td><td>0.1275146</td><td> 1884</td><td>103</td></tr>
	<tr><td>5</td><td>0.1639672</td><td> 1606</td><td> 18</td></tr>
	<tr><td>6</td><td>0.2137556</td><td>  546</td><td> 38</td></tr>
	<tr><td>7</td><td>0.2335350</td><td>  389</td><td>  5</td></tr>
	<tr><td>8</td><td>0.2128419</td><td>  311</td><td> 25</td></tr>
</tbody>
</table>




```R
Idents(immune.combined) <- 'integrated_snn_res.0.1'
new.cluster.ids<-c("T cells","Smooth muscle cells","Monocytes/Macrophages", "Endothelial cells","B cells", 
                   "Natural killer cells", "Mast cells", "Proliferating cells", "Plasma cells")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined,reduction="umap",label=T)+NoLegend()
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_23_0.png)
    



```R
immune.combined[['celltype']] <- Idents(immune.combined)
```


```R
immune.combined@meta.data <- immune.combined@meta.data %>%
                                 mutate(celltype_Brief = fct_recode(celltype,
                                        'Mono/Mø' = 'Monocytes/Macrophages',
                                        'T' = 'T cells',
                                        'B' = 'B cells',       
                                       'SMC'= 'Smooth muscle cells',
                                        'EC'= 'Endothelial cells',
                                       'NK'= 'Natural killer cells',
                                    'PLF' = 'Proliferating cells',
                                      'MC'  = 'Mast cells',
                                    'PC'    = 'Plasma cells'                              
                                        ))
immune.combined@meta.data$celltype_Brief <- factor(immune.combined@meta.data$celltype_Brief,
                                                  levels = c('Mono/Mø','SMC','EC',
                                                             'T','NK','PLF',
                                                             'B','PC','MC'))
```


```R
DimPlot(immune.combined,group.by ='celltype_Brief',label=T)+NoLegend()
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_26_0.png)
    



```R
saveRDS(immune.combined, file = 'Step1_AddDoublet_Seurat.Rds')
```


```R
levels(Idents(immune.combined))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'T cells'</li><li>'Smooth muscle cells'</li><li>'Monocytes/Macrophages'</li><li>'Endothelial cells'</li><li>'B cells'</li><li>'Natural killer cells'</li><li>'Mast cells'</li><li>'Proliferating cells'</li><li>'Plasma cells'</li></ol>




```R
BuenColors::tf1_color_maps[c(1:12)] %>% unname()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#F6313E'</li><li>'#46A040'</li><li>'#0081C9'</li><li>'#A65AC2'</li><li>'#FFA300'</li><li>'#FFFF32'</li><li>'#89774A'</li><li>'#FF6A80'</li><li>'#999999'</li><li>'#0DB2AA'</li><li>'#001588'</li><li>'#00441B'</li></ol>




```R
color <- c('#46A040','#0081C9','#F6313E','#FFA300','#A65AC2','#89774A','#00441B','#999999','#0DB2AA')
names(color) <- levels(Idents(immune.combined))
```


```R
DimPlot(immune.combined, group.by ='celltype')  + scale_color_manual(values = color)
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_31_0.png)
    



```R
p1 <- DimPlot(immune.combined, group.by ='celltype')  + theme_void() + 
NoLegend() + scale_color_manual(values = color) + labs(title = '')
p1
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_32_0.png)
    



```R
immune.combined@meta.data$group
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>'AC'</li><li>⋯</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li><li>'PA'</li></ol>




```R
p2 <- DimPlot(immune.combined, group.by ='group')  + theme_void() + 
NoLegend() + scale_color_viridis(begin = 0.4, end = 1,discrete = TRUE) + labs(title = '')
p2
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_34_0.png)
    



```R
p3 <- DimPlot(immune.combined, group.by ='sample')  + theme_void() + 
NoLegend()  + labs(title = '')
p3
ggsave(p3 ,file ='Seurat_group.pdf', height = 10, width = 10)
```


    
![png](Step10.2_Set_Idents_files/Step10.2_Set_Idents_35_0.png)
    



```R
library(patchwork)
p0 <- p1 + p2 
ggsave(p0,file ='Seurat_celltype.pdf', width = 20, height = 10)
```
