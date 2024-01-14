```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(doParallel)
library(RColorBrewer)
library(DoubletFinder)
```


```R
GSM4837523 <- Read10X(data.dir = '~/AS/GSE159677/GSM4837523/')
```


```R
GSM4837524 <- Read10X(data.dir = '~/AS/GSE159677/GSM4837524/')
```


```R
GSM4837525 <- Read10X(data.dir = '~/AS/GSE159677/GSM4837525/')
```


```R
GSM4837526 <- Read10X(data.dir = '~/AS/GSE159677/GSM4837526/')
```


```R
GSM4837527 <- Read10X(data.dir = '~/AS/GSE159677/GSM4837527/')
```


```R
GSM4837528 <- Read10X(data.dir = '~/AS/GSE159677/GSM4837528/')
```


```R
Count.list <- list(GSM4837523, GSM4837524, GSM4837525, GSM4837526, GSM4837527, GSM4837528)
```


```R
colnames(Count.list[[1]]) <- str_c(colnames(Count.list[[1]]), "_AC1")
colnames(Count.list[[2]]) <- str_c(colnames(Count.list[[2]]), "_PA1")
colnames(Count.list[[3]]) <- str_c(colnames(Count.list[[3]]), "_AC2")
colnames(Count.list[[4]]) <- str_c(colnames(Count.list[[4]]), "_PA2")
colnames(Count.list[[5]]) <- str_c(colnames(Count.list[[5]]), "_AC3")
colnames(Count.list[[6]]) <- str_c(colnames(Count.list[[6]]), "_PA3")
```


```R
#create_metadata
sam.info.list <- vector("list",6) 
```


```R
group_type <- rep(c("AC","PA"),3)
```


```R
group_type
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'AC'</li><li>'PA'</li><li>'AC'</li><li>'PA'</li><li>'AC'</li><li>'PA'</li></ol>




```R
for (i in 1:6){
sam.info.list[[i]] <- as.data.frame(matrix(data=NA, nrow=ncol(Count.list[[i]]),ncol=2))
rownames(sam.info.list[[i]]) <- colnames(Count.list[[i]])
colnames(sam.info.list[[i]]) <- c("sample","group")
sam.info.list[[i]][,1] <- rep(paste0("patient_", (i+1) %/% 2,group_type[[i]]), ncol(Count.list[[i]]))
sam.info.list[[i]][,2] <- rep(group_type[[i]], ncol(Count.list[[i]]))
}
```


```R
Seurat.list <- vector("list",6)
```


```R
for (i in 1:6) {
  Seurat.list[[i]] <- CreateSeuratObject(
  counts = Count.list[[i]], 
  meta.data=sam.info.list[[i]],
  project = "Human_carotid_plaque", 
  min.cells = 3, min.features = 50,
  names.field = 2,names.delim = "_")
 }
```


```R
immune.combined.list <- lapply(X = Seurat.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
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
    
    	Found 8024 anchors
    
    Filtering anchors
    
    	Retained 4611 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 19713 anchors
    
    Filtering anchors
    
    	Retained 8764 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 8593 anchors
    
    Filtering anchors
    
    	Retained 3545 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 11612 anchors
    
    Filtering anchors
    
    	Retained 6430 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 6758 anchors
    
    Filtering anchors
    
    	Retained 4421 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 16906 anchors
    
    Filtering anchors
    
    	Retained 11457 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 19515 anchors
    
    Filtering anchors
    
    	Retained 9489 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 8597 anchors
    
    Filtering anchors
    
    	Retained 3711 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 25292 anchors
    
    Filtering anchors
    
    	Retained 11711 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 13624 anchors
    
    Filtering anchors
    
    	Retained 5940 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 8268 anchors
    
    Filtering anchors
    
    	Retained 5597 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 7068 anchors
    
    Filtering anchors
    
    	Retained 5245 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 9522 anchors
    
    Filtering anchors
    
    	Retained 6963 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 6945 anchors
    
    Filtering anchors
    
    	Retained 5437 anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 9625 anchors
    
    Filtering anchors
    
    	Retained 6726 anchors
    



```R
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
```

    Merging dataset 4 into 3
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    
    Merging dataset 6 into 5
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    
    Merging dataset 2 into 1
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    
    Merging dataset 5 6 into 3 4
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    
    Merging dataset 1 2 into 3 4 5 6
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    



```R
DefaultAssay(immune.combined) <- 'RNA'
```


```R
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern = "^MT-")
```


```R
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_22_0.png)
    



```R
quantile(immune.combined@meta.data$nCount_RNA,  c(.25, .50,  .75, .90, .95,.99))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>25%</dt><dd>2526</dd><dt>50%</dt><dd>4473</dd><dt>75%</dt><dd>8324</dd><dt>90%</dt><dd>12710.4</dd><dt>95%</dt><dd>16337</dd><dt>99%</dt><dd>29757.2399999998</dd></dl>




```R
quantile(immune.combined@meta.data$nFeature_RNA,  c(.25, .50,  .75, .90, .95,.99))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>25%</dt><dd>998</dd><dt>50%</dt><dd>1501</dd><dt>75%</dt><dd>2331</dd><dt>90%</dt><dd>2981</dd><dt>95%</dt><dd>3386.2</dd><dt>99%</dt><dd>4225.24</dd></dl>




```R
plot1 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```


    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_25_0.png)
    



```R
immune.combined <- subset(immune.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA > 500 & nCount_RNA <20000
               & percent.mt < 25 )
```


```R
plot1 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(immune.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```


    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_27_0.png)
    



```R
##CCAIntegration
immune.combined.list <-SplitObject(immune.combined, split.by = "sample")
```


```R
doublet.list <- lapply(X = immune.combined.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x)
    x <- RunPCA(x)
    x <- RunUMAP(x, dims = 1:20)
})
```

    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  TYROBP, FCER1G, AIF1, COTL1, LYZ, CTSS, PLEK, CD14, CD68, C5AR1 
    	   TYMP, CD83, GPR183, PLAUR, PLIN2, FCGR3A, MS4A7, BCL2A1, RGCC, LST1 
    	   ZNF331, GLUL, RNASET2, C1orf162, SAT1, MNDA, HLA-DRA, IFI30, HLA-DPA1, CTSH 
    Negative:  BGN, IGFBP7, MGP, CALD1, ADIRF, CLU, CRIP2, CAV1, NNMT, SELENOM 
    	   MYL9, OGN, IGFBP2, TPM2, EFEMP1, CTGF, C1R, TAGLN, CAVIN1, SPARCL1 
    	   TPM1, COX7A1, DSTN, EMP2, ID3, CRYAB, SOD3, MFGE8, CD151, CAVIN3 
    PC_ 2 
    Positive:  CST3, PSAP, FTL, AIF1, CD68, CTSZ, CD63, CD14, C5AR1, CTSB 
    	   MARCKS, TYROBP, FCER1G, LYZ, CEBPD, SAT1, PLAUR, MS4A7, BLVRB, GLUL 
    	   HLA-DRB5, CTSS, TIMP1, IER3, HLA-DRA, CTSL, ASAH1, MAFB, TYMP, CFD 
    Negative:  IL32, CD2, GZMA, TRBC2, TRAC, CD3E, CD69, CCL5, CD3D, CD3G 
    	   CST7, CTSW, GZMK, NKG7, TRBC1, GZMM, APOBEC3G, GIMAP7, LTB, CD247 
    	   ITM2A, SPOCK2, CD7, IL7R, GZMH, DUSP2, KLRB1, CXCR3, SOCS1, HOPX 
    PC_ 3 
    Positive:  SOD3, ACTA2, FRZB, RAMP1, RGS5, FXYD1, MFAP4, TAGLN, LMOD1, ITGA8 
    	   MYH11, TPM2, COL14A1, MYL9, CRYAB, PPP1R14A, C11orf96, NOV, MAP1B, ID4 
    	   PDE5A, MYH10, MFGE8, SSPN, PLN, C2orf40, CNN1, NEXN, GEM, TPM1 
    Negative:  CLEC14A, ECSCR, VWF, RAMP2, SOX18, BMX, EDN1, CLDN5, CALCRL, PALMD 
    	   MPZL2, PECAM1, EGFL7, PTPRB, PROCR, HSPG2, NRN1, SELP, DKK2, HYAL2 
    	   NPDC1, ITLN1, ADGRL4, IGFBP4, MTUS1, PLVAP, TSPAN7, EFNA1, CDH5, AC004540.2 
    PC_ 4 
    Positive:  SPP1, SLAMF9, CD36, CSTB, MARCO, FABP5, SDC2, MMP19, PHLDA1, FABP4 
    	   SMIM25, FBP1, C15orf48, S100A10, FN1, AQP9, CCL7, GPNMB, PLIN2, ERO1A 
    	   ANPEP, TIMP1, VCAN, CYP27A1, CTSL, ENO1, P4HA1, RNASE1, CD109, MGLL 
    Negative:  HLA-DQA1, MS4A6A, C1QA, HLA-DPA1, HLA-DPB1, C1QB, C1QC, HLA-DQB1, HLA-DMB, HLA-DMA 
    	   F13A1, FGL2, FOLR2, CD74, TMEM176B, GPR34, CPVL, FCGR1B, SELENOP, TMEM176A 
    	   FCGR1A, HLA-DRA, HLA-DRB1, HLA-DRB5, P2RY13, CCL4L2, CCL3L1, CTSC, RASSF4, IGF1 
    PC_ 5 
    Positive:  LUM, COL1A1, THY1, COL6A3, COL3A1, POSTN, COL1A2, FAP, COL15A1, COL5A2 
    	   ECM1, PRSS23, PCOLCE, SERPINF1, COL5A1, COL6A1, DCN, MXRA5, PDPN, KDELR3 
    	   MMP11, COL4A1, IBSP, COL4A2, SPARC, COL6A2, LRRC17, TMEM119, ISLR, COL12A1 
    Negative:  PLN, MFAP4, TCEAL2, FRZB, ACTC1, MYH11, RAMP1, CNN1, LMOD1, CLU 
    	   ACTA2, ITGA8, LMCD1, SCX, RGS5, RERGL, MYH10, CSRP2, NEXN, SLC25A4 
    	   PPP1R14A, ITLN1, HSPB7, SYNPO2, SORBS2, CRIM1, LDOC1, HRCT1, PPP1R12B, SELP 
    
    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    20:57:15 UMAP embedding parameters a = 0.9922 b = 1.112
    
    20:57:15 Read 10446 rows and found 20 numeric columns
    
    20:57:15 Using Annoy for neighbor search, n_neighbors = 30
    
    20:57:15 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    20:57:16 Writing NN index file to temp file /tmp/Rtmpea1M1m/file4f753d4557e
    
    20:57:16 Searching Annoy index using 1 thread, search_k = 3000
    
    20:57:18 Annoy recall = 100%
    
    20:57:18 Commencing smooth kNN distance calibration using 1 thread
    
    20:57:18 Initializing from normalized Laplacian + noise
    
    20:57:19 Commencing optimization for 200 epochs, with 434918 positive edges
    
    20:57:22 Optimization finished
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  CALD1, SOD3, OGN, CRYAB, MGP, MYL9, EFEMP1, PRELP, IGFBP5, COL6A2 
    	   IGFBP7, FXYD1, BGN, COL1A2, TPM1, AEBP1, TAGLN, MFAP4, DSTN, CTGF 
    	   IGFBP6, FHL1, TPM2, SSPN, OMD, PLAC9, MFGE8, FBLN5, PTGIS, MAP1B 
    Negative:  SRGN, CD74, HLA-DPA1, HLA-DPB1, PTPRC, HLA-DRA, TYROBP, HLA-DRB1, LAPTM5, HCST 
    	   ITGB2, HLA-DQB1, FCER1G, LCP1, HLA-DRB5, VAMP8, CORO1A, PLEK, CD37, CD53 
    	   CTSS, CTSC, GPSM3, COTL1, AIF1, HLA-DMA, FYB1, CCL3, EFHD2, CD83 
    PC_ 2 
    Positive:  S100A4, CYBA, TYROBP, EMP3, FCER1G, AIF1, PLEK, ITGB2, CD68, CYBB 
    	   C5AR1, FGL2, LYZ, SPI1, GLIPR1, MNDA, PTPRC, MS4A6A, IGSF6, LST1 
    	   LGALS1, RGS10, LAPTM5, HCST, LCP1, CCL3, LSP1, SERPINA1, CLEC7A, RNASET2 
    Negative:  ACKR1, ECSCR, SPARCL1, RAMP3, IFI27, RAMP2, PLVAP, CLEC14A, VWF, IGFBP4 
    	   TM4SF1, PALMD, CCL14, ADGRL4, NRN1, CALCRL, EMCN, SOX18, STC1, SLC9A3R2 
    	   RGS16, SPRY1, HYAL2, EGFL7, OLFM1, RNASE1, C2CD4B, TNFSF10, MARCKSL1, ITGA6 
    PC_ 3 
    Positive:  ACTA2, PPP1R14A, MYH11, RAMP1, TPM2, RGS5, ITGA8, MYH10, C12orf75, MAP1B 
    	   RCAN2, FRZB, IGFBP2, FBLIM1, NEXN, CNN1, FLNA, TAGLN, NOV, NKG7 
    	   GZMA, CTSW, LMOD1, CST7, PDLIM3, CCL5, TRBC2, PDE5A, SOST, BGN 
    Negative:  CFD, SERPINF1, DCN, MGST1, SFRP2, C3, FBLN1, MFAP5, APOD, PI16 
    	   PLA2G2A, LRRN4CL, KLF4, FBLN2, CST3, SLPI, IGF1, SCARA5, CHRDL1, WISP2 
    	   GSN, ACKR3, PDGFRL, LEPR, ABCA8, CEBPD, C1S, LUM, GAS1, ABI3BP 
    PC_ 4 
    Positive:  CST3, AIF1, MS4A6A, HLA-DRA, CD14, C5AR1, CYBB, SPI1, IGSF6, CLEC7A 
    	   LYZ, MNDA, HLA-DQB1, CD83, MS4A7, FCGR2A, HLA-DRB5, SERPINA1, IL1B, HLA-DMA 
    	   CD86, LST1, LY86, LMO2, HLA-DQA1, FCGR1A, CXCL8, PLAUR, CSF1R, C1QC 
    Negative:  CTSW, GZMA, NKG7, CST7, CCL5, TRBC2, GZMM, PRF1, CD3E, CD247 
    	   GZMH, KLRD1, CD2, IL2RG, TRBC1, IL32, GNLY, TRAC, CD7, GZMB 
    	   SKAP1, LCK, RAC2, CD69, CD3G, FGFBP2, CD3D, SAMD3, HOPX, KLRB1 
    PC_ 5 
    Positive:  S100A9, S100A8, FCN1, S100A12, CSTA, AC020656.1, EREG, CFP, LILRA5, SLC11A1 
    	   APOBEC3A, CLEC12A, LGALS2, SMIM25, MCEMP1, C19orf38, LILRA2, LYZ, LILRA1, AC245128.3 
    	   G0S2, HCAR3, AREG, CSF3R, S100A6, SERPINA1, NCF2, ANPEP, RETN, FGR 
    Negative:  C1QA, C1QC, C1QB, FOLR2, MRC1, MS4A4A, VSIG4, GPR34, MSR1, VMO1 
    	   ADAP2, SLCO2B1, CSF1R, IGSF21, MAF, SELENOP, F13A1, MS4A7, CD163, CH25H 
    	   FCGBP, FCGR2A, STAB1, TREM2, CCL4L2, HPGDS, CCL8, EBI3, C3AR1, LGMN 
    
    20:57:24 UMAP embedding parameters a = 0.9922 b = 1.112
    
    20:57:24 Read 3341 rows and found 20 numeric columns
    
    20:57:24 Using Annoy for neighbor search, n_neighbors = 30
    
    20:57:24 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    20:57:24 Writing NN index file to temp file /tmp/Rtmpea1M1m/file4f71fcf1179
    
    20:57:24 Searching Annoy index using 1 thread, search_k = 3000
    
    20:57:25 Annoy recall = 100%
    
    20:57:25 Commencing smooth kNN distance calibration using 1 thread
    
    20:57:25 Initializing from normalized Laplacian + noise
    
    20:57:25 Commencing optimization for 500 epochs, with 131612 positive edges
    
    20:57:28 Optimization finished
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  SRGN, CYBA, ITGB2, GPR183, CCL5, LTB, RGS1, CCL4, GZMK, IL7R 
    	   RGCC, VAMP8, COTL1, ZNF331, NKG7, TRBC1, CD74, HLA-DPB1, UCP2, RNASET2 
    	   HLA-DQB1, TNF, CD83, BIRC3, HLA-DQA1, CCL4L2, HLA-DPA1, CCR7, CTSW, HLA-DQA2 
    Negative:  IGFBP7, CALD1, ADIRF, BGN, MGP, NNMT, CAVIN3, SPARC, SPARCL1, CTGF 
    	   CYR61, MYL9, CFH, SOD3, IFITM3, TAGLN, C1R, COL6A2, COL14A1, AEBP1 
    	   SERPING1, CAV1, TPM1, CPE, IGFBP2, CD9, TPM2, TIMP1, CNN3, COL1A2 
    PC_ 2 
    Positive:  AIF1, CD14, MS4A6A, C5AR1, CD68, FCER1G, FCGR2A, IGSF6, TYROBP, TMEM176B 
    	   C1QC, CPVL, MS4A7, C1QA, LYZ, CFD, SPI1, C1QB, F13A1, CYBB 
    	   FCGR1A, CLEC7A, CTSB, CD163, SERPINA1, MNDA, PLAUR, MAFB, CSF1R, FTL 
    Negative:  CCL5, GZMK, IL7R, LTB, TRBC1, NKG7, IGKC, KLRB1, CTSW, GZMH 
    	   CD8A, CCR7, CD7, CD8B, C12orf75, HOPX, TNF, TIGIT, PRF1, IFNG 
    	   TRGC2, KLRD1, XCL2, CXCR6, IGLC2, CD40LG, HIST1H1D, GNLY, GZMB, FKBP11 
    PC_ 3 
    Positive:  RAMP2, CLEC14A, ECSCR, PLVAP, ADGRL4, PALMD, SOX18, EMCN, PECAM1, EGFL7 
    	   HYAL2, IFI27, CALCRL, VWA1, GNG11, RBP7, RAMP3, ID1, CYYR1, CDH5 
    	   TM4SF18, VWF, HSPG2, FLT1, CXorf36, ESAM, LDB2, FAM110D, KANK3, PRCP 
    Negative:  S100A4, COL14A1, SOD3, TAGLN, TPM2, COL1A2, AEBP1, C1R, CRYAB, PCOLCE 
    	   CPE, C1S, OGN, DKK3, LUM, ACTA2, BGN, VCAN, FRZB, LTBP1 
    	   SSPN, MFGE8, PRRX1, FXYD1, THBS2, FN1, ID4, SMOC2, LGALS3BP, C2orf40 
    PC_ 4 
    Positive:  RGS5, COX4I2, CSRP2, ABCC9, MCAM, TBX2, NDUFA4L2, TBX2-AS1, RERGL, CCDC102B 
    	   HIGD1B, FAM162B, SEPT4, NOTCH3, GUCY1B1, CDH6, PGF, PTP4A3, GJA4, ACTA2 
    	   C11orf96, KCNE4, EDNRA, ISYNA1, TINAGL1, TPPP3, FILIP1L, GUCY1A2, KCNJ8, HES4 
    Negative:  CRTAC1, LUM, SFRP4, SFRP2, FAP, LTBP2, DPT, DCN, TGM2, OMD 
    	   ITGBL1, SMOC1, SLPI, SULF1, KDELR3, CTSK, VCAN, SUGCT, GAP43, CP 
    	   PCOLCE2, MMP23B, CDH11, MXRA5, HAPLN1, THBS2, COL3A1, RCN3, PRSS23, COL8A1 
    PC_ 5 
    Positive:  CD79A, BANK1, MS4A1, RALGPS2, TNFRSF13C, VPREB3, SPIB, CD24, IGHM, CD83 
    	   TNFRSF13B, HLA-DRA, FAM30A, HLA-DQA1, LINC01781, LINC01857, BASP1, CD74, IRF8, HLA-DQB1 
    	   HLA-DPA1, HLA-DQA2, HLA-DPB1, HLA-DRB5, LY86, HLA-DRB1, MARCKS, CD79B, IGHA2, MARCH1 
    Negative:  CCL5, NKG7, GZMK, RGCC, S100A4, CCL4, SRGN, S100A10, CTSW, MT2A 
    	   S100A11, ITGB2, TRBC1, GZMH, S100A6, VIM, IL7R, CD7, CD8A, KLRB1 
    	   PRF1, CD8B, FOS, XCL2, KLRD1, CSTB, GNLY, TXN, PLIN2, CITED2 
    
    20:57:35 UMAP embedding parameters a = 0.9922 b = 1.112
    
    20:57:35 Read 15074 rows and found 20 numeric columns
    
    20:57:35 Using Annoy for neighbor search, n_neighbors = 30
    
    20:57:35 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    20:57:36 Writing NN index file to temp file /tmp/Rtmpea1M1m/file4f745204a1a
    
    20:57:36 Searching Annoy index using 1 thread, search_k = 3000
    
    20:57:39 Annoy recall = 100%
    
    20:57:40 Commencing smooth kNN distance calibration using 1 thread
    
    20:57:40 Initializing from normalized Laplacian + noise
    
    20:57:41 Commencing optimization for 200 epochs, with 650156 positive edges
    
    20:57:45 Optimization finished
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  SRGN, CD69, DUSP2, TRBC2, IL32, CYBA, LTB, CCL5, ITGB2, RGCC 
    	   GZMA, COTL1, IL7R, CCL4, GPR183, GZMK, ALOX5AP, TRBC1, NKG7, RGS1 
    	   CCR7, VAMP8, CD7, CTSW, MALAT1, TNF, NFKBIA, NCF1, CD8A, CD27 
    Negative:  IGFBP7, MGP, SERPING1, IFITM3, CALD1, NNMT, GSN, SPARCL1, PLAC9, ADIRF 
    	   CST3, TIMP1, CPE, NUPR1, CD9, C1R, SOD3, CYR61, LHFPL6, CAVIN3 
    	   COL14A1, DSTN, EGR1, CTGF, TAGLN, CD151, COX7A1, SPARC, BGN, CLU 
    PC_ 2 
    Positive:  S100A4, SOD3, COL14A1, TAGLN, TPM2, FXYD1, CRYAB, BGN, ACTA2, ID4 
    	   RARRES2, FRZB, SMOC2, MFAP4, OGN, CPE, LGALS1, MYH11, AEBP1, LMOD1 
    	   TPM1, NEXN, C12orf75, MAP1B, PCOLCE, C1S, LTBP1, COL1A2, MFGE8, C1R 
    Negative:  RNASE1, ECSCR, CLEC14A, PLVAP, ACKR1, RAMP2, PECAM1, VWF, CCL14, EMCN 
    	   EGFL7, ADGRL4, AQP1, IFI27, PALMD, SOX18, RAMP3, TSPAN7, CAVIN2, HYAL2 
    	   JAM2, RBP7, MMRN1, MMRN2, CYYR1, CXorf36, NRN1, POSTN, CD93, NOSTRIN 
    PC_ 3 
    Positive:  BCAM, TINAGL1, MCAM, RERGL, CRIP2, ADIRF, CAV1, ESAM, CNN1, MYH11 
    	   CSRP2, ADAMTS1, CAV2, ID1, SORBS2, ACTA2, PPP1R14A, PLN, SNCG, SPARCL1 
    	   MALAT1, C11orf96, SELENOW, FHL5, IL32, ISYNA1, HES4, MYL9, PTP4A3, C12orf75 
    Negative:  AIF1, CD68, CD14, MS4A6A, LYZ, C5AR1, FCER1G, CLEC7A, TMEM176B, CYBB 
    	   SERPINA1, FCGR1A, IGSF6, TYROBP, MS4A7, SPI1, CD163, PLAUR, MNDA, AC020656.1 
    	   C1QC, C1QA, FCGR2A, F13A1, FGL2, C1QB, LILRB4, CTSS, PILRA, CSF1R 
    PC_ 4 
    Positive:  RERGL, CNN1, NRGN, ACTA2, PPP1R14A, MYH11, PLN, MCAM, SORBS2, FHL5 
    	   HES4, ACTG2, TPM2, CSRP2, MYL9, TINAGL1, RCAN2, TBX2-AS1, CASQ2, BCAM 
    	   WTIP, COX4I2, EFHD1, FILIP1L, ISYNA1, NDUFA4L2, LMOD1, RHOB, C11orf96, MYLK 
    Negative:  FBLN1, PI16, C3, MFAP5, PLA2G2A, CXCL14, DCN, SFRP4, OMD, SFRP2 
    	   CHRDL1, SFRP1, SCARA5, LRRN4CL, ABCA8, PDGFRL, MMP2, EFEMP1, RARRES1, SLIT2 
    	   PTGIS, ITGBL1, WISP2, GAS1, LUM, CCDC80, TNXB, SEMA3C, ABCA6, ADH1B 
    PC_ 5 
    Positive:  S100A8, S100A12, S100A9, FCN1, CSTA, SLC11A1, APOBEC3A, CFP, CLEC12A, NCF2 
    	   CSF3R, AC020656.1, FPR1, LST1, THBS1, NEAT1, FOS, CD300E, MGST1, LYZ 
    	   TREM1, CFD, PLBD1, SMIM25, MNDA, MXD1, NAMPT, FGR, SULT1A1, S100A6 
    Negative:  HLA-DQA1, CD79A, MS4A1, BANK1, HLA-DRA, CD83, HLA-DQB1, HLA-DQA2, TNFRSF13C, HLA-DPB1 
    	   HLA-DPA1, RALGPS2, SPIB, CD74, MARCKS, HLA-DRB1, VPREB3, HLA-DRB5, HLA-DMA, ADAM28 
    	   TNFRSF13B, CD24, HLA-DMB, FCRLA, BLNK, IGHM, LINC00926, IGKC, IRF8, LINC02397 
    
    20:57:48 UMAP embedding parameters a = 0.9922 b = 1.112
    
    20:57:48 Read 5190 rows and found 20 numeric columns
    
    20:57:48 Using Annoy for neighbor search, n_neighbors = 30
    
    20:57:48 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    20:57:48 Writing NN index file to temp file /tmp/Rtmpea1M1m/file4f7537a18ad
    
    20:57:48 Searching Annoy index using 1 thread, search_k = 3000
    
    20:57:49 Annoy recall = 100%
    
    20:57:50 Commencing smooth kNN distance calibration using 1 thread
    
    20:57:50 Initializing from normalized Laplacian + noise
    
    20:57:50 Commencing optimization for 500 epochs, with 210428 positive edges
    
    20:57:54 Optimization finished
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  CALD1, ADIRF, IGFBP7, SPARCL1, MGP, BGN, FN1, CAV1, CAVIN1, SELENOM 
    	   DSTN, TAGLN, MFGE8, OGN, MYL9, TPM1, TPM2, MAP1B, SOD3, FHL1 
    	   COL14A1, COL6A2, AEBP1, PLS3, CAVIN3, LHFPL6, ACTA2, CLU, CD151, SPARC 
    Negative:  TYROBP, AIF1, FCER1G, HLA-DPA1, HLA-DPB1, HLA-DRA, CYBA, CD74, MS4A6A, HLA-DRB1 
    	   CTSS, HLA-DQB1, CD14, CYBB, HLA-DRB5, CPVL, LYZ, SRGN, MNDA, IGSF6 
    	   PLEK, HLA-DMA, HLA-DQA1, HLA-DMB, FCGR2A, SPI1, LST1, CD83, MS4A7, CLEC7A 
    PC_ 2 
    Positive:  CD52, CD69, IL7R, IL32, TRBC2, TRAC, CCL5, GZMA, LTB, CD7 
    	   TRBC1, DUSP2, GZMK, CST7, CTSW, NKG7, GZMH, SOCS1, RGCC, KLRB1 
    	   ICOS, CCR7, HOPX, ITM2A, APOBEC3G, GPR171, CD40LG, CD27, KLRG1, BIRC3 
    Negative:  CST3, IFITM3, PSAP, FTL, NPC2, CTSZ, CD9, CEBPD, FCGRT, GRN 
    	   FTH1, LGALS3, GSN, MARCKS, FGL2, RHOB, GAS6, TIMP1, CD68, THBD 
    	   CD14, IER3, BLVRB, KLF4, CTSB, MAFB, CD59, HSPB1, HLA-DRB5, RNF130 
    PC_ 3 
    Positive:  VWF, RAMP2, ECSCR, MPZL2, CLEC14A, PECAM1, IFI27, PTPRB, SELP, PALMD 
    	   CALCRL, EGFL7, HYAL2, MMRN2, DKK2, PLVAP, BMX, SOX18, GJA5, RAMP3 
    	   ID1, CLDN5, RNASE1, PODXL, F5, CDH5, EDN1, ACVRL1, ASS1, FGF18 
    Negative:  ACTA2, MYH11, C11orf96, TAGLN, LMOD1, MYL9, PLN, SOD3, RGS5, ITGA8 
    	   MFAP4, TPM2, PPP1R14A, COL14A1, MYH10, FXYD1, RAMP1, MYLK, FRZB, PALLD 
    	   COL6A2, SYNPO2, LTBP1, C12orf75, FHL1, TPM1, CNN1, SMOC2, SVIL, A2M 
    PC_ 4 
    Positive:  PLN, CNN1, RAMP1, SORBS2, CSRP2, LMOD1, CPE, MCAM, MYH11, SYNPO2 
    	   MYLK, C12orf75, LDOC1, FILIP1L, ACTG2, TCEAL2, PPP1R12B, PCDH7, CARMN, MFAP4 
    	   TNS1, PRUNE2, FHL5, TIMP3, FLNA, ITGA8, SBSPON, ACTA2, ADRA2C, KCNMA1 
    Negative:  LUM, CFH, NDUFA4L2, PLAC9, COL15A1, PCOLCE2, AGT, GGT5, GAP43, PCOLCE 
    	   APOE, STEAP1, FAP, FBLN1, CCDC80, TFPI2, C7, PTN, F2R, LRRC32 
    	   COL1A2, PRSS23, THBS2, COL1A1, FBLN2, SPRY1, STEAP4, KDELR3, CCDC102B, CLEC11A 
    PC_ 5 
    Positive:  FCN1, S100A8, S100A9, RETN, SLC11A1, CSTA, CFP, S100A12, G0S2, TNFRSF1B 
    	   JAML, HCAR3, C19orf38, IL32, CD300E, IL1B, IL1RN, SERPINA1, CSF3R, ITGB2 
    	   APOBEC3A, HCAR2, AC020656.1, THBS1, FGR, CLEC4E, CCL5, AQP9, FPR1, CARD16 
    Negative:  TPSAB1, CPA3, TPSB2, MS4A2, HDC, RHEX, KIT, HPGDS, SLC18A2, RGS13 
    	   GATA2, HPGD, MLPH, VWA5A, IL1RL1, TPSD1, CALB2, AL157895.1, GCSAML, SLC45A3 
    	   STXBP6, KRT19, LIF, RAB27B, MAOB, CPM, SMYD3, PTGS1, CADPS, BTK 
    
    20:58:00 UMAP embedding parameters a = 0.9922 b = 1.112
    
    20:58:00 Read 11462 rows and found 20 numeric columns
    
    20:58:00 Using Annoy for neighbor search, n_neighbors = 30
    
    20:58:00 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    20:58:01 Writing NN index file to temp file /tmp/Rtmpea1M1m/file4f75710e42b
    
    20:58:01 Searching Annoy index using 1 thread, search_k = 3000
    
    20:58:03 Annoy recall = 100%
    
    20:58:03 Commencing smooth kNN distance calibration using 1 thread
    
    20:58:04 Initializing from normalized Laplacian + noise
    
    20:58:04 Commencing optimization for 200 epochs, with 482112 positive edges
    
    20:58:07 Optimization finished
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  SRGN, LAPTM5, PTPRC, TYROBP, ITGB2, FCER1G, CD37, LCP1, HCST, AIF1 
    	   CTSS, EVI2B, CD53, FYB1, CORO1A, VAMP8, SMAP2, GPSM3, LYZ, SPI1 
    	   CYBB, COTL1, C1orf162, PLEK, MS4A6A, C5AR1, CXCR4, MNDA, STK17B, LST1 
    Negative:  MGP, GSN, TIMP3, C1R, IGFBP7, IGFBP5, SERPING1, COL6A2, CPE, IGFBP6 
    	   C1S, COL1A2, PLAC9, CALD1, MFAP4, OGN, CCDC80, PRELP, FBLN1, DCN 
    	   EFEMP1, PRRX1, LTBP4, PCOLCE, AEBP1, SOD3, NUPR1, COL6A1, MXRA8, FBLN2 
    PC_ 2 
    Positive:  S100A4, CFD, CYBA, LGALS1, FGL2, LRP1, CD68, SERPINF1, FTL, MGST1 
    	   FBLN1, DCN, PLTP, VCAN, GPNMB, PLBD1, SFRP2, MFAP5, COL1A2, MFAP4 
    	   C1S, C3, GLIPR1, PI16, AIF1, PDGFRL, TYROBP, IGFBP6, CTSB, CEBPB 
    Negative:  RAMP3, ECSCR, ACKR1, IFI27, PLVAP, BCAM, AQP1, CLEC14A, ADGRL4, VWF 
    	   SNCG, NRN1, EMCN, PECAM1, CCL14, SLC9A3R2, TSPAN7, JAM2, CYYR1, GIMAP7 
    	   EPAS1, CAV1, MMRN2, MMRN1, PALMD, ESAM, CRIP2, PDLIM1, CAVIN2, TM4SF1 
    PC_ 3 
    Positive:  CST3, CEBPD, HLA-DRA, HLA-DRB5, HLA-DRB1, HLA-DPA1, CD14, RNASE1, CD74, HLA-DQB1 
    	   HLA-DMA, CTSZ, SELENOP, MARCKS, THBD, HLA-DPB1, CD68, FTL, PDK4, MAFB 
    	   IGFBP4, HLA-DQA1, CYBB, MS4A6A, AIF1, BLVRB, FOS, EMP1, CTSB, NAMPT 
    Negative:  CD52, CD2, CD3E, CD69, TRBC2, CD3D, CD3G, TRAC, CTSW, CCL5 
    	   IL7R, CD247, LCK, NKG7, IL2RG, IL32, GZMA, CST7, CD7, GZMM 
    	   SPOCK2, PRF1, GZMH, PRDM1, RUNX3, SKAP1, TRBC1, DUSP2, ARL4C, SAMD3 
    PC_ 4 
    Positive:  ACTA2, MYH11, TPM2, PPP1R14A, PLN, CNN1, MYLK, RCAN2, FRZB, RERGL 
    	   RGS5, LMOD1, TAGLN, RAMP1, GUCY1A1, ITGA8, ACTG2, NEXN, MYL9, NOTCH3 
    	   SYNPO2, GUCY1B1, HES4, FHL5, MRVI1, C11orf96, CAP2, ID4, TBX2, MCAM 
    Negative:  FBLN2, PLAT, SLPI, TNXB, RAMP2, MFAP5, CD34, PI16, C3, DCN 
    	   SCARA5, SFRP2, APOD, CFD, CHRDL1, TNFSF10, LRRN4CL, FBN1, FBLN1, LEPR 
    	   MGST1, PCOLCE2, GAS1, SERPINF1, CD248, PTGIS, EMP1, SEMA3C, ABCA8, PRSS23 
    PC_ 5 
    Positive:  C1QA, MRC1, C1QB, C1QC, LYVE1, FOLR2, F13A1, MS4A4A, LILRB5, P2RY14 
    	   CCL13, SLCO2B1, STAB1, VSIG4, GPR34, RNASE1, SELENOP, DAB2, LGMN, MAF 
    	   SCN9A, CD209, FCGR2B, CSF1R, MAMDC2, CCL18, PLTP, CD163, CLEC10A, MARCO 
    Negative:  FCN1, S100A8, CSTA, S100A12, BCL2A1, S100A9, SLC11A1, SERPINA1, FGR, APOBEC3A 
    	   AC020656.1, LILRA5, C19orf38, FOLR3, VCAN, FPR1, CLEC4E, LILRA2, BASP1, ASGR1 
    	   CLEC12A, CSF3R, JAML, AC245128.3, CD300E, CDA, HCAR3, FPR2, PRAM1, HCAR2 
    
    20:58:09 UMAP embedding parameters a = 0.9922 b = 1.112
    
    20:58:09 Read 3111 rows and found 20 numeric columns
    
    20:58:09 Using Annoy for neighbor search, n_neighbors = 30
    
    20:58:09 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    20:58:09 Writing NN index file to temp file /tmp/Rtmpea1M1m/file4f7756d7aa0
    
    20:58:09 Searching Annoy index using 1 thread, search_k = 3000
    
    20:58:10 Annoy recall = 100%
    
    20:58:10 Commencing smooth kNN distance calibration using 1 thread
    
    20:58:10 Initializing from normalized Laplacian + noise
    
    20:58:11 Commencing optimization for 500 epochs, with 122626 positive edges
    
    20:58:13 Optimization finished
    



```R
## pK Identification (no ground-truth) 
PK.list <- lapply(X = doublet.list, FUN = function(x) {
    x %>%
    paramSweep_v3(PCs = 1:20, sct = FALSE) %>%
    summarizeSweep(GT = FALSE) %>% 
    find.pK()
    })
PK.list
```

    Loading required package: KernSmooth
    
    KernSmooth 2.23 loaded
    Copyright M. P. Wand 1997-2009
    
    Loading required package: ROCR
    
    Loading required package: fields
    
    Loading required package: spam
    
    Spam version 2.8-0 (2022-01-05) is loaded.
    Type 'help( Spam)' or 'demo( spam)' for a short introduction 
    and overview of this package.
    Help for individual functions is also obtained by adding the
    suffix '.spam' to the function name, e.g. 'help( chol.spam)'.
    
    
    Attaching package: ‘spam’
    
    
    The following objects are masked from ‘package:base’:
    
        backsolve, forwardsolve
    
    
    Loading required package: viridis
    
    Loading required package: viridisLite
    
    See https://github.com/NCAR/Fields for
     an extensive vignette, other supplements and source code 
    


    [1] "Creating artificial doublets for pN = 5%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 10%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 15%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 20%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 25%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 30%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    NULL
    [1] "Creating artificial doublets for pN = 5%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 10%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 15%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 20%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 25%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 30%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    NULL
    [1] "Creating artificial doublets for pN = 5%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 10%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 15%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 20%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 25%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 30%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    



    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_30_37.png)
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    NULL
    [1] "Creating artificial doublets for pN = 5%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 10%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 15%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 20%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 25%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 30%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    



    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_30_50.png)
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    NULL
    [1] "Creating artificial doublets for pN = 5%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 10%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 15%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 20%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 25%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 30%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    



    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_30_63.png)
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.001..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    NULL
    [1] "Creating artificial doublets for pN = 5%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 10%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 15%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 20%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 25%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."
    [1] "Creating artificial doublets for pN = 30%"
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    



    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_30_76.png)
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Defining neighborhoods..."
    [1] "Computing pANN across all pK..."
    [1] "pK = 0.005..."
    [1] "pK = 0.01..."
    [1] "pK = 0.02..."
    [1] "pK = 0.03..."
    [1] "pK = 0.04..."
    [1] "pK = 0.05..."
    [1] "pK = 0.06..."
    [1] "pK = 0.07..."
    [1] "pK = 0.08..."
    [1] "pK = 0.09..."
    [1] "pK = 0.1..."
    [1] "pK = 0.11..."
    [1] "pK = 0.12..."
    [1] "pK = 0.13..."
    [1] "pK = 0.14..."
    [1] "pK = 0.15..."
    [1] "pK = 0.16..."
    [1] "pK = 0.17..."
    [1] "pK = 0.18..."
    [1] "pK = 0.19..."
    [1] "pK = 0.2..."
    [1] "pK = 0.21..."
    [1] "pK = 0.22..."
    [1] "pK = 0.23..."
    [1] "pK = 0.24..."
    [1] "pK = 0.25..."
    [1] "pK = 0.26..."
    [1] "pK = 0.27..."
    [1] "pK = 0.28..."
    [1] "pK = 0.29..."
    [1] "pK = 0.3..."



    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_30_78.png)
    


    NULL



<dl>
	<dt>$patient_1AC</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 32 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.001</td><td>0.8373111</td><td>0.0014745571</td><td> 567.83907</td></tr>
	<tr><td> 2</td><td>0.005</td><td>0.8594746</td><td>0.0005775851</td><td>1488.04855</td></tr>
	<tr><td> 3</td><td>0.01 </td><td>0.8421643</td><td>0.0010560462</td><td> 797.46913</td></tr>
	<tr><td> 4</td><td>0.02 </td><td>0.8365642</td><td>0.0010256731</td><td> 815.62462</td></tr>
	<tr><td> 5</td><td>0.03 </td><td>0.7956493</td><td>0.0019302868</td><td> 412.19225</td></tr>
	<tr><td> 6</td><td>0.04 </td><td>0.7698482</td><td>0.0019210824</td><td> 400.73671</td></tr>
	<tr><td> 7</td><td>0.05 </td><td>0.7392876</td><td>0.0008954368</td><td> 825.61666</td></tr>
	<tr><td> 8</td><td>0.06 </td><td>0.7231233</td><td>0.0013187637</td><td> 548.33425</td></tr>
	<tr><td> 9</td><td>0.07 </td><td>0.6892395</td><td>0.0018511581</td><td> 372.32883</td></tr>
	<tr><td>10</td><td>0.08 </td><td>0.6537851</td><td>0.0047860887</td><td> 136.60112</td></tr>
	<tr><td>11</td><td>0.09 </td><td>0.6420521</td><td>0.0094743820</td><td>  67.76718</td></tr>
	<tr><td>12</td><td>0.1  </td><td>0.6475758</td><td>0.0147940400</td><td>  43.77275</td></tr>
	<tr><td>13</td><td>0.11 </td><td>0.6428277</td><td>0.0188063690</td><td>  34.18138</td></tr>
	<tr><td>14</td><td>0.12 </td><td>0.6100023</td><td>0.0169718035</td><td>  35.94210</td></tr>
	<tr><td>15</td><td>0.13 </td><td>0.5677926</td><td>0.0112435110</td><td>  50.49958</td></tr>
	<tr><td>16</td><td>0.14 </td><td>0.5460106</td><td>0.0070212768</td><td>  77.76514</td></tr>
	<tr><td>17</td><td>0.15 </td><td>0.5502016</td><td>0.0052455252</td><td> 104.88971</td></tr>
	<tr><td>18</td><td>0.16 </td><td>0.5577076</td><td>0.0081999303</td><td>  68.01370</td></tr>
	<tr><td>19</td><td>0.17 </td><td>0.5748519</td><td>0.0123356194</td><td>  46.60098</td></tr>
	<tr><td>20</td><td>0.18 </td><td>0.6087052</td><td>0.0121070725</td><td>  50.27682</td></tr>
	<tr><td>21</td><td>0.19 </td><td>0.6413800</td><td>0.0094893084</td><td>  67.58975</td></tr>
	<tr><td>22</td><td>0.2  </td><td>0.6641036</td><td>0.0073661215</td><td>  90.15648</td></tr>
	<tr><td>23</td><td>0.21 </td><td>0.6805899</td><td>0.0058194963</td><td> 116.94996</td></tr>
	<tr><td>24</td><td>0.22 </td><td>0.6900450</td><td>0.0041728119</td><td> 165.36691</td></tr>
	<tr><td>25</td><td>0.23 </td><td>0.6945832</td><td>0.0026912871</td><td> 258.08587</td></tr>
	<tr><td>26</td><td>0.24 </td><td>0.6990678</td><td>0.0016782648</td><td> 416.54201</td></tr>
	<tr><td>27</td><td>0.25 </td><td>0.7011283</td><td>0.0010213255</td><td> 686.48855</td></tr>
	<tr><td>28</td><td>0.26 </td><td>0.6989349</td><td>0.0007351561</td><td> 950.72982</td></tr>
	<tr><td>29</td><td>0.27 </td><td>0.6936984</td><td>0.0007299806</td><td> 950.29702</td></tr>
	<tr><td>30</td><td>0.28 </td><td>0.6881199</td><td>0.0007850198</td><td> 876.56373</td></tr>
	<tr><td>31</td><td>0.29 </td><td>0.6820341</td><td>0.0007804946</td><td> 873.84858</td></tr>
	<tr><td>32</td><td>0.3  </td><td>0.6747380</td><td>0.0007386577</td><td> 913.46515</td></tr>
</tbody>
</table>
</dd>
	<dt>$patient_1PA</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 31 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.005</td><td>0.8676896</td><td>0.0007448023</td><td>1164.99313</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.8771596</td><td>0.0014311100</td><td> 612.92254</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.8720755</td><td>0.0006044080</td><td>1442.85909</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.8061417</td><td>0.0023727883</td><td> 339.74446</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.7602599</td><td>0.0029026257</td><td> 261.92145</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.7184457</td><td>0.0029356614</td><td> 244.73045</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.6870218</td><td>0.0033244479</td><td> 206.65742</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.6563817</td><td>0.0032953466</td><td> 199.18442</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.6176146</td><td>0.0036941257</td><td> 167.18830</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.6266496</td><td>0.0043035253</td><td> 145.61308</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.6079600</td><td>0.0037015115</td><td> 164.24643</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.6077280</td><td>0.0046305542</td><td> 131.24305</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.6322509</td><td>0.0016370480</td><td> 386.21402</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.6415031</td><td>0.0007086957</td><td> 905.18835</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.6020815</td><td>0.0043097361</td><td> 139.70264</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.5620379</td><td>0.0064967207</td><td>  86.51101</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.5194640</td><td>0.0065288443</td><td>  79.56446</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.4875983</td><td>0.0043013218</td><td> 113.36011</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.4803638</td><td>0.0026685813</td><td> 180.00719</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.4953827</td><td>0.0019472010</td><td> 254.40760</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.5121360</td><td>0.0010177865</td><td> 503.18613</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.5413422</td><td>0.0010056329</td><td> 538.30990</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.5679720</td><td>0.0009384114</td><td> 605.24840</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.5886728</td><td>0.0010674037</td><td> 551.49971</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.6067686</td><td>0.0013950102</td><td> 434.95640</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.6272155</td><td>0.0022060758</td><td> 284.31275</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.6454837</td><td>0.0026288002</td><td> 245.54308</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.6672401</td><td>0.0030459773</td><td> 219.05617</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.6896171</td><td>0.0029643850</td><td> 232.63411</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.7119314</td><td>0.0024939664</td><td> 285.46152</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.7304068</td><td>0.0019852882</td><td> 367.90970</td></tr>
</tbody>
</table>
</dd>
	<dt>$patient_2AC</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 32 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.001</td><td>0.7687462</td><td>0.0023102584</td><td> 332.75332</td></tr>
	<tr><td> 2</td><td>0.005</td><td>0.8150234</td><td>0.0014493398</td><td> 562.34116</td></tr>
	<tr><td> 3</td><td>0.01 </td><td>0.7895742</td><td>0.0027025941</td><td> 292.15421</td></tr>
	<tr><td> 4</td><td>0.02 </td><td>0.7004752</td><td>0.0011153519</td><td> 628.03068</td></tr>
	<tr><td> 5</td><td>0.03 </td><td>0.6234873</td><td>0.0015475881</td><td> 402.87676</td></tr>
	<tr><td> 6</td><td>0.04 </td><td>0.5844172</td><td>0.0024276084</td><td> 240.73785</td></tr>
	<tr><td> 7</td><td>0.05 </td><td>0.5385351</td><td>0.0044828039</td><td> 120.13354</td></tr>
	<tr><td> 8</td><td>0.06 </td><td>0.4986511</td><td>0.0075697693</td><td>  65.87401</td></tr>
	<tr><td> 9</td><td>0.07 </td><td>0.4606719</td><td>0.0082515014</td><td>  55.82886</td></tr>
	<tr><td>10</td><td>0.08 </td><td>0.4263205</td><td>0.0066213941</td><td>  64.38531</td></tr>
	<tr><td>11</td><td>0.09 </td><td>0.4175223</td><td>0.0059909461</td><td>  69.69222</td></tr>
	<tr><td>12</td><td>0.1  </td><td>0.3981467</td><td>0.0059017295</td><td>  67.46271</td></tr>
	<tr><td>13</td><td>0.11 </td><td>0.4166131</td><td>0.0072750899</td><td>  57.26570</td></tr>
	<tr><td>14</td><td>0.12 </td><td>0.4498435</td><td>0.0098733012</td><td>  45.56161</td></tr>
	<tr><td>15</td><td>0.13 </td><td>0.4826397</td><td>0.0097514029</td><td>  49.49438</td></tr>
	<tr><td>16</td><td>0.14 </td><td>0.5133422</td><td>0.0080734155</td><td>  63.58426</td></tr>
	<tr><td>17</td><td>0.15 </td><td>0.5353425</td><td>0.0053305907</td><td> 100.42837</td></tr>
	<tr><td>18</td><td>0.16 </td><td>0.5578465</td><td>0.0030634469</td><td> 182.09765</td></tr>
	<tr><td>19</td><td>0.17 </td><td>0.5665135</td><td>0.0008001878</td><td> 707.97564</td></tr>
	<tr><td>20</td><td>0.18 </td><td>0.5629122</td><td>0.0009907939</td><td> 568.14258</td></tr>
	<tr><td>21</td><td>0.19 </td><td>0.5523770</td><td>0.0012100125</td><td> 456.50518</td></tr>
	<tr><td>22</td><td>0.2  </td><td>0.5446792</td><td>0.0008807176</td><td> 618.44924</td></tr>
	<tr><td>23</td><td>0.21 </td><td>0.5452244</td><td>0.0005871947</td><td> 928.52401</td></tr>
	<tr><td>24</td><td>0.22 </td><td>0.5530786</td><td>0.0004822528</td><td>1146.86463</td></tr>
	<tr><td>25</td><td>0.23 </td><td>0.5642727</td><td>0.0005133157</td><td>1099.27029</td></tr>
	<tr><td>26</td><td>0.24 </td><td>0.5788079</td><td>0.0005775139</td><td>1002.24051</td></tr>
	<tr><td>27</td><td>0.25 </td><td>0.5885011</td><td>0.0006791518</td><td> 866.52355</td></tr>
	<tr><td>28</td><td>0.26 </td><td>0.5984500</td><td>0.0006635567</td><td> 901.88220</td></tr>
	<tr><td>29</td><td>0.27 </td><td>0.6023284</td><td>0.0006415999</td><td> 938.79130</td></tr>
	<tr><td>30</td><td>0.28 </td><td>0.6079386</td><td>0.0005411124</td><td>1123.49791</td></tr>
	<tr><td>31</td><td>0.29 </td><td>0.6139712</td><td>0.0005348663</td><td>1147.89653</td></tr>
	<tr><td>32</td><td>0.3  </td><td>0.6244827</td><td>0.0004035014</td><td>1547.65928</td></tr>
</tbody>
</table>
</dd>
	<dt>$patient_2PA</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 31 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.005</td><td>0.8496838</td><td>0.0018856408</td><td> 450.60744</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.8500711</td><td>0.0015255161</td><td> 557.23509</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.8099062</td><td>0.0017266672</td><td> 469.05750</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.7338759</td><td>0.0065390750</td><td> 112.22931</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.6673431</td><td>0.0104721660</td><td>  63.72541</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.6191896</td><td>0.0122040061</td><td>  50.73659</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.5644412</td><td>0.0105405042</td><td>  53.54973</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.4947129</td><td>0.0047137586</td><td> 104.95084</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.5895244</td><td>0.0060463914</td><td>  97.50020</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.6088345</td><td>0.0110631700</td><td>  55.03255</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.5696556</td><td>0.0174136006</td><td>  32.71326</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.5452143</td><td>0.0133891151</td><td>  40.72071</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.5309218</td><td>0.0072586285</td><td>  73.14354</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.5763603</td><td>0.0070229307</td><td>  82.06834</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.6060984</td><td>0.0059041604</td><td> 102.65615</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.6260666</td><td>0.0038559580</td><td> 162.36345</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.6260355</td><td>0.0044548387</td><td> 140.52932</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.6170760</td><td>0.0038510250</td><td> 160.23683</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.6124080</td><td>0.0021100374</td><td> 290.23562</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.6168539</td><td>0.0014595884</td><td> 422.62184</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.6316144</td><td>0.0028370479</td><td> 222.63086</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.6457766</td><td>0.0032871244</td><td> 196.45640</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.6715824</td><td>0.0035437765</td><td> 189.51037</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.6905317</td><td>0.0022906446</td><td> 301.45738</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.7047296</td><td>0.0012247801</td><td> 575.39276</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.7115530</td><td>0.0007860030</td><td> 905.28032</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.7125783</td><td>0.0005625598</td><td>1266.67116</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.7116028</td><td>0.0003364687</td><td>2114.91519</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.7123424</td><td>0.0002951723</td><td>2413.31013</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.7117106</td><td>0.0004653235</td><td>1529.49644</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.7136110</td><td>0.0004254221</td><td>1677.41868</td></tr>
</tbody>
</table>
</dd>
	<dt>$patient_3AC</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 32 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.001</td><td>0.8062212</td><td>0.0013839482</td><td> 582.5516</td></tr>
	<tr><td> 2</td><td>0.005</td><td>0.8888258</td><td>0.0001952709</td><td>4551.7580</td></tr>
	<tr><td> 3</td><td>0.01 </td><td>0.8997737</td><td>0.0002850822</td><td>3156.1909</td></tr>
	<tr><td> 4</td><td>0.02 </td><td>0.8772050</td><td>0.0037334347</td><td> 234.9593</td></tr>
	<tr><td> 5</td><td>0.03 </td><td>0.8471291</td><td>0.0031145639</td><td> 271.9896</td></tr>
	<tr><td> 6</td><td>0.04 </td><td>0.8288300</td><td>0.0022204460</td><td> 373.2718</td></tr>
	<tr><td> 7</td><td>0.05 </td><td>0.7663232</td><td>0.0022674333</td><td> 337.9695</td></tr>
	<tr><td> 8</td><td>0.06 </td><td>0.6820327</td><td>0.0020980173</td><td> 325.0844</td></tr>
	<tr><td> 9</td><td>0.07 </td><td>0.5963261</td><td>0.0007725093</td><td> 771.9338</td></tr>
	<tr><td>10</td><td>0.08 </td><td>0.5371119</td><td>0.0003347259</td><td>1604.6317</td></tr>
	<tr><td>11</td><td>0.09 </td><td>0.4990212</td><td>0.0002753980</td><td>1812.0002</td></tr>
	<tr><td>12</td><td>0.1  </td><td>0.4754324</td><td>0.0003510778</td><td>1354.2080</td></tr>
	<tr><td>13</td><td>0.11 </td><td>0.4751878</td><td>0.0004147130</td><td>1145.8232</td></tr>
	<tr><td>14</td><td>0.12 </td><td>0.5010039</td><td>0.0016643732</td><td> 301.0165</td></tr>
	<tr><td>15</td><td>0.13 </td><td>0.5361650</td><td>0.0039471672</td><td> 135.8354</td></tr>
	<tr><td>16</td><td>0.14 </td><td>0.5508938</td><td>0.0044695802</td><td> 123.2540</td></tr>
	<tr><td>17</td><td>0.15 </td><td>0.5501522</td><td>0.0034692377</td><td> 158.5801</td></tr>
	<tr><td>18</td><td>0.16 </td><td>0.5502657</td><td>0.0030494749</td><td> 180.4460</td></tr>
	<tr><td>19</td><td>0.17 </td><td>0.5562534</td><td>0.0023012486</td><td> 241.7181</td></tr>
	<tr><td>20</td><td>0.18 </td><td>0.5627614</td><td>0.0020957697</td><td> 268.5226</td></tr>
	<tr><td>21</td><td>0.19 </td><td>0.5693632</td><td>0.0024612920</td><td> 231.3270</td></tr>
	<tr><td>22</td><td>0.2  </td><td>0.5742787</td><td>0.0022781011</td><td> 252.0866</td></tr>
	<tr><td>23</td><td>0.21 </td><td>0.5762261</td><td>0.0020606013</td><td> 279.6398</td></tr>
	<tr><td>24</td><td>0.22 </td><td>0.5793543</td><td>0.0019509094</td><td> 296.9663</td></tr>
	<tr><td>25</td><td>0.23 </td><td>0.5839414</td><td>0.0018454374</td><td> 316.4244</td></tr>
	<tr><td>26</td><td>0.24 </td><td>0.5858197</td><td>0.0014320873</td><td> 409.0670</td></tr>
	<tr><td>27</td><td>0.25 </td><td>0.5834517</td><td>0.0010357264</td><td> 563.3261</td></tr>
	<tr><td>28</td><td>0.26 </td><td>0.5790894</td><td>0.0008950585</td><td> 646.9849</td></tr>
	<tr><td>29</td><td>0.27 </td><td>0.5772421</td><td>0.0007735784</td><td> 746.1972</td></tr>
	<tr><td>30</td><td>0.28 </td><td>0.5793554</td><td>0.0007730574</td><td> 749.4339</td></tr>
	<tr><td>31</td><td>0.29 </td><td>0.5861385</td><td>0.0008835943</td><td> 663.3570</td></tr>
	<tr><td>32</td><td>0.3  </td><td>0.5947770</td><td>0.0008398194</td><td> 708.2201</td></tr>
</tbody>
</table>
</dd>
	<dt>$patient_3PA</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 31 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.005</td><td>0.8781046</td><td>0.0008875073</td><td> 989.40554</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.8943893</td><td>0.0009736613</td><td> 918.58359</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.8768132</td><td>0.0017840588</td><td> 491.47101</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.8225369</td><td>0.0006775457</td><td>1213.99473</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.7863365</td><td>0.0012740442</td><td> 617.19723</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.7239593</td><td>0.0020146826</td><td> 359.34161</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.7036003</td><td>0.0010911986</td><td> 644.79579</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.6733706</td><td>0.0006885224</td><td> 977.99373</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.6480002</td><td>0.0041881055</td><td> 154.72393</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.6277410</td><td>0.0130450942</td><td>  48.12085</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.6654605</td><td>0.0148587515</td><td>  44.78576</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.7017789</td><td>0.0164433104</td><td>  42.67869</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.7165313</td><td>0.0137675702</td><td>  52.04486</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.7116847</td><td>0.0172815364</td><td>  41.18179</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.7134267</td><td>0.0115320723</td><td>  61.86457</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.7138329</td><td>0.0085841495</td><td>  83.15709</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.6979957</td><td>0.0060501886</td><td> 115.36760</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.6769357</td><td>0.0045684103</td><td> 148.17753</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.6356593</td><td>0.0038758402</td><td> 164.00555</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.6016074</td><td>0.0022085468</td><td> 272.39966</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.5665781</td><td>0.0020399500</td><td> 277.74116</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.5781480</td><td>0.0025094484</td><td> 230.38846</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.6133849</td><td>0.0028777780</td><td> 213.14531</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.6785944</td><td>0.0008241481</td><td> 823.38888</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.7401294</td><td>0.0014471350</td><td> 511.44462</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.7451066</td><td>0.0002083602</td><td>3576.05098</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.7266831</td><td>0.0009375986</td><td> 775.04717</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.7082304</td><td>0.0019494250</td><td> 363.30218</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.7009133</td><td>0.0027276003</td><td> 256.97069</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.6971730</td><td>0.0027547238</td><td> 253.08273</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.6978114</td><td>0.0021300094</td><td> 327.60955</td></tr>
</tbody>
</table>
</dd>
</dl>




    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_30_81.png)
    



```R
#find optimalPK
opt_pk_list <- lapply (X = PK.list, FUN = function(x) { 
    as.numeric(as.vector(x$pK[which.max(x$BCmetric)]))
})
opt_pk_list
```


<dl>
	<dt>$patient_1AC</dt>
		<dd>0.005</dd>
	<dt>$patient_1PA</dt>
		<dd>0.02</dd>
	<dt>$patient_2AC</dt>
		<dd>0.3</dd>
	<dt>$patient_2PA</dt>
		<dd>0.28</dd>
	<dt>$patient_3AC</dt>
		<dd>0.005</dd>
	<dt>$patient_3PA</dt>
		<dd>0.25</dd>
</dl>




```R
#find optimal_pk for 1:9 by artifical
i = 2
PK.list[[i]]
ggplot(PK.list[[i]], mapping = aes(x = pK, y = BCmetric)) + 
geom_col() + 
theme(axis.text.x = element_text(angle=45))
```


<table class="dataframe">
<caption>A data.frame: 31 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.005</td><td>0.8676896</td><td>0.0007448023</td><td>1164.99313</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.8771596</td><td>0.0014311100</td><td> 612.92254</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.8720755</td><td>0.0006044080</td><td>1442.85909</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.8061417</td><td>0.0023727883</td><td> 339.74446</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.7602599</td><td>0.0029026257</td><td> 261.92145</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.7184457</td><td>0.0029356614</td><td> 244.73045</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.6870218</td><td>0.0033244479</td><td> 206.65742</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.6563817</td><td>0.0032953466</td><td> 199.18442</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.6176146</td><td>0.0036941257</td><td> 167.18830</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.6266496</td><td>0.0043035253</td><td> 145.61308</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.6079600</td><td>0.0037015115</td><td> 164.24643</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.6077280</td><td>0.0046305542</td><td> 131.24305</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.6322509</td><td>0.0016370480</td><td> 386.21402</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.6415031</td><td>0.0007086957</td><td> 905.18835</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.6020815</td><td>0.0043097361</td><td> 139.70264</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.5620379</td><td>0.0064967207</td><td>  86.51101</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.5194640</td><td>0.0065288443</td><td>  79.56446</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.4875983</td><td>0.0043013218</td><td> 113.36011</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.4803638</td><td>0.0026685813</td><td> 180.00719</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.4953827</td><td>0.0019472010</td><td> 254.40760</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.5121360</td><td>0.0010177865</td><td> 503.18613</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.5413422</td><td>0.0010056329</td><td> 538.30990</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.5679720</td><td>0.0009384114</td><td> 605.24840</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.5886728</td><td>0.0010674037</td><td> 551.49971</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.6067686</td><td>0.0013950102</td><td> 434.95640</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.6272155</td><td>0.0022060758</td><td> 284.31275</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.6454837</td><td>0.0026288002</td><td> 245.54308</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.6672401</td><td>0.0030459773</td><td> 219.05617</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.6896171</td><td>0.0029643850</td><td> 232.63411</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.7119314</td><td>0.0024939664</td><td> 285.46152</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.7304068</td><td>0.0019852882</td><td> 367.90970</td></tr>
</tbody>
</table>




    
![png](Step10.1_Human_Plaques_identify_files/Step10.1_Human_Plaques_identify_32_1.png)
    



```R
opt_pk_list[[1]] <-0.05
opt_pk_list[[2]] <-0.13
opt_pk_list[[3]] <-0.02
opt_pk_list[[4]] <-0.19
opt_pk_list[[5]] <-0.09
opt_pk_list[[6]] <-0.03
```


```R
homotypic.list <- lapply(X = doublet.list, FUN = function(x) {
        x@meta.data$celltype3 %>%
        modelHomotypic()})
homotypic.list 
```


<dl>
	<dt>$patient_1AC</dt>
		<dd>0</dd>
	<dt>$patient_1PA</dt>
		<dd>0</dd>
	<dt>$patient_2AC</dt>
		<dd>0</dd>
	<dt>$patient_2PA</dt>
		<dd>0</dd>
	<dt>$patient_3AC</dt>
		<dd>0</dd>
	<dt>$patient_3PA</dt>
		<dd>0</dd>
</dl>




```R
doublet.list[[1]]@assays
```


    $RNA
    Assay data with 22621 features for 10446 cells
    Top 10 variable features:
     ITLN1, IBSP, LUM, MT1G, PTGDS, COL1A1, CCL18, MT1H, MMP9, MMP7 
    
    $integrated
    Assay data with 3000 features for 10446 cells
    Top 10 variable features:
     ITLN1, JCHAIN, PTGDS, CCL18, S100A8, S100A9, SPP1, CXCL10, SFRP4,
    TPSAB1 




```R
#doubletrate-0.05
nExp_poi.list <- lapply(X = doublet.list, FUN = function(x) {
    round(0.05*ncol(x@assays$integrated@data))
    })
```


```R
nExp_poi.adj.list <- nExp_poi.list
for (i in 1:6) {
  nExp_poi.adj.list[[i]] <- round(nExp_poi.list[[i]]*(1-homotypic.list[[i]]))   
}
```


```R
Exp_poi.adj.num <- round(as.numeric(nExp_poi.list)*(1-as.numeric(homotypic.list)))
```


```R
for (i in 1:6) {
    doublet.list[[i]] <- doubletFinder_v3(doublet.list[[i]],
    PCs = 1:20,
    pN = 0.25,
    pK = as.numeric(opt_pk_list[[i]]),
    nExp = as.numeric(nExp_poi.adj.list[[i]]),
    reuse.pANN = FALSE,
    sct = FALSE                                         
    )}
```

    [1] "Creating 3482 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."
    [1] "Creating 1114 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."
    [1] "Creating 5025 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."
    [1] "Creating 1730 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."
    [1] "Creating 3821 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."
    [1] "Creating 1037 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."



```R
doublet.list[[1]]@meta.data
```


<table class="dataframe">
<caption>A data.frame: 10446 × 8</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_RNA</th><th scope=col>nFeature_RNA</th><th scope=col>sample</th><th scope=col>group</th><th scope=col>percent.mt</th><th scope=col>pANN_0.25_0.05_522</th><th scope=col>DF.classifications_0.25_0.05_522</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACCCAAGATTAGAC-1_AC1</th><td>AC1</td><td>11667</td><td>2905</td><td>patient_1AC</td><td>AC</td><td> 2.6056398</td><td>0.36063218</td><td>Doublet</td></tr>
	<tr><th scope=row>AAACCCAAGCATGTTC-1_AC1</th><td>AC1</td><td> 5385</td><td>1494</td><td>patient_1AC</td><td>AC</td><td> 3.8625812</td><td>0.09482759</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCAAGCCTGTCG-1_AC1</th><td>AC1</td><td> 1471</td><td> 575</td><td>patient_1AC</td><td>AC</td><td> 7.3419443</td><td>0.01149425</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCAAGGGTTTCT-1_AC1</th><td>AC1</td><td> 8957</td><td>2515</td><td>patient_1AC</td><td>AC</td><td> 1.8198057</td><td>0.23850575</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCAAGTGCACCC-1_AC1</th><td>AC1</td><td>11681</td><td>2586</td><td>patient_1AC</td><td>AC</td><td> 1.8919613</td><td>0.23994253</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCACAAGTTTGC-1_AC1</th><td>AC1</td><td> 2495</td><td> 794</td><td>patient_1AC</td><td>AC</td><td> 1.9639279</td><td>0.02873563</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCACATTCTCCG-1_AC1</th><td>AC1</td><td> 8930</td><td>2326</td><td>patient_1AC</td><td>AC</td><td> 1.9148936</td><td>0.18678161</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCAGTGTGTGTT-1_AC1</th><td>AC1</td><td>19264</td><td>3317</td><td>patient_1AC</td><td>AC</td><td> 1.9362542</td><td>0.27729885</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCAGTGTTTCTT-1_AC1</th><td>AC1</td><td>  932</td><td> 476</td><td>patient_1AC</td><td>AC</td><td> 0.1072961</td><td>0.10775862</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCAGTTACCCTC-1_AC1</th><td>AC1</td><td>10108</td><td>2366</td><td>patient_1AC</td><td>AC</td><td> 1.9885239</td><td>0.18534483</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCATCTCACTCG-1_AC1</th><td>AC1</td><td>  809</td><td> 519</td><td>patient_1AC</td><td>AC</td><td> 6.9221261</td><td>0.01293103</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCATCTCATTAC-1_AC1</th><td>AC1</td><td> 1798</td><td> 856</td><td>patient_1AC</td><td>AC</td><td> 3.0033370</td><td>0.09482759</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCCATCTCTCTAA-1_AC1</th><td>AC1</td><td>12381</td><td>2838</td><td>patient_1AC</td><td>AC</td><td> 1.8496083</td><td>0.11637931</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAAAGAAGTGTT-1_AC1</th><td>AC1</td><td> 9902</td><td>2348</td><td>patient_1AC</td><td>AC</td><td> 3.7871137</td><td>0.20114943</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAAAGAGTGAAG-1_AC1</th><td>AC1</td><td> 5238</td><td>1515</td><td>patient_1AC</td><td>AC</td><td> 5.7464681</td><td>0.36781609</td><td>Doublet</td></tr>
	<tr><th scope=row>AAACGAACACAGCCTG-1_AC1</th><td>AC1</td><td> 4594</td><td>1638</td><td>patient_1AC</td><td>AC</td><td> 5.6377884</td><td>0.20977011</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAACACTGGACC-1_AC1</th><td>AC1</td><td> 5729</td><td>1743</td><td>patient_1AC</td><td>AC</td><td> 2.0596963</td><td>0.31321839</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAAGTAGGCTGA-1_AC1</th><td>AC1</td><td> 3168</td><td>1172</td><td>patient_1AC</td><td>AC</td><td> 3.2196970</td><td>0.10919540</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAAGTCGTTGCG-1_AC1</th><td>AC1</td><td>14693</td><td>2959</td><td>patient_1AC</td><td>AC</td><td> 0.3130743</td><td>0.20258621</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAAGTTATGTGC-1_AC1</th><td>AC1</td><td> 8378</td><td>1997</td><td>patient_1AC</td><td>AC</td><td> 1.9813798</td><td>0.24856322</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAAGTTCTCGTC-1_AC1</th><td>AC1</td><td> 3601</td><td>1445</td><td>patient_1AC</td><td>AC</td><td> 2.5826159</td><td>0.31178161</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAAGTTTACCAG-1_AC1</th><td>AC1</td><td>  607</td><td> 386</td><td>patient_1AC</td><td>AC</td><td>10.3789127</td><td>0.14224138</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAATCGATTCCC-1_AC1</th><td>AC1</td><td>15027</td><td>2654</td><td>patient_1AC</td><td>AC</td><td> 2.3091768</td><td>0.31178161</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGAATCTCTTGCG-1_AC1</th><td>AC1</td><td> 1437</td><td> 699</td><td>patient_1AC</td><td>AC</td><td> 3.4794711</td><td>0.10344828</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGCTAGGAGCAAA-1_AC1</th><td>AC1</td><td> 2794</td><td> 980</td><td>patient_1AC</td><td>AC</td><td> 8.6614173</td><td>0.08764368</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGCTAGGTTCCAT-1_AC1</th><td>AC1</td><td> 4126</td><td>1293</td><td>patient_1AC</td><td>AC</td><td> 2.1085797</td><td>0.07040230</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGCTAGTCCCGGT-1_AC1</th><td>AC1</td><td> 2692</td><td>1054</td><td>patient_1AC</td><td>AC</td><td> 4.3462110</td><td>0.10344828</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGCTAGTTACGTC-1_AC1</th><td>AC1</td><td>11864</td><td>2568</td><td>patient_1AC</td><td>AC</td><td> 1.6604855</td><td>0.18390805</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGCTCACAAGCTT-1_AC1</th><td>AC1</td><td> 8898</td><td>1994</td><td>patient_1AC</td><td>AC</td><td> 2.3038885</td><td>0.18534483</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGCTGTAGCTTGT-1_AC1</th><td>AC1</td><td> 2732</td><td>1303</td><td>patient_1AC</td><td>AC</td><td> 3.3308931</td><td>0.22844828</td><td>Singlet</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>TTTGGAGCACTATCGA-1_AC1</th><td>AC1</td><td>  673</td><td> 238</td><td>patient_1AC</td><td>AC</td><td>13.075780</td><td>0.01436782</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGCAGAACGCA-1_AC1</th><td>AC1</td><td> 2347</td><td> 986</td><td>patient_1AC</td><td>AC</td><td> 1.832126</td><td>0.08764368</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGCAGCTGTCG-1_AC1</th><td>AC1</td><td>12090</td><td>2841</td><td>patient_1AC</td><td>AC</td><td> 2.109181</td><td>0.25287356</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGCAGCTTCCT-1_AC1</th><td>AC1</td><td> 2419</td><td>1008</td><td>patient_1AC</td><td>AC</td><td> 3.472509</td><td>0.16810345</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGGTACAACGG-1_AC1</th><td>AC1</td><td> 3657</td><td>1423</td><td>patient_1AC</td><td>AC</td><td> 3.254033</td><td>0.26005747</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGGTAGATCCT-1_AC1</th><td>AC1</td><td> 5653</td><td>1766</td><td>patient_1AC</td><td>AC</td><td> 2.600389</td><td>0.07327586</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGGTCAATCTG-1_AC1</th><td>AC1</td><td>14866</td><td>3315</td><td>patient_1AC</td><td>AC</td><td> 3.235571</td><td>0.23419540</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGGTCTTCCGT-1_AC1</th><td>AC1</td><td> 1847</td><td> 939</td><td>patient_1AC</td><td>AC</td><td> 3.898213</td><td>0.08333333</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGTCCACATAG-1_AC1</th><td>AC1</td><td> 2432</td><td> 900</td><td>patient_1AC</td><td>AC</td><td> 8.840461</td><td>0.15948276</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGTCCCGAGTG-1_AC1</th><td>AC1</td><td> 6993</td><td>1659</td><td>patient_1AC</td><td>AC</td><td> 3.160303</td><td>0.26436782</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGTCGCAGTTA-1_AC1</th><td>AC1</td><td> 6251</td><td>1579</td><td>patient_1AC</td><td>AC</td><td> 1.983683</td><td>0.12787356</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGAGTCTCCCATG-1_AC1</th><td>AC1</td><td>  761</td><td> 459</td><td>patient_1AC</td><td>AC</td><td>12.483574</td><td>0.18821839</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTAGTCATGCT-1_AC1</th><td>AC1</td><td> 6134</td><td>1819</td><td>patient_1AC</td><td>AC</td><td> 4.988588</td><td>0.28735632</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTCAATGTCTG-1_AC1</th><td>AC1</td><td> 4385</td><td>1564</td><td>patient_1AC</td><td>AC</td><td> 2.098062</td><td>0.10057471</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTCATTAAGCC-1_AC1</th><td>AC1</td><td> 2300</td><td>1003</td><td>patient_1AC</td><td>AC</td><td> 1.043478</td><td>0.08764368</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTGTGGCCTCA-1_AC1</th><td>AC1</td><td> 9926</td><td>2615</td><td>patient_1AC</td><td>AC</td><td> 3.324602</td><td>0.19971264</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTTCGATACTG-1_AC1</th><td>AC1</td><td>10969</td><td>2336</td><td>patient_1AC</td><td>AC</td><td> 2.288267</td><td>0.22844828</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTTCTCCATAT-1_AC1</th><td>AC1</td><td>14297</td><td>2814</td><td>patient_1AC</td><td>AC</td><td> 1.895503</td><td>0.30603448</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGAGCATGCGA-1_AC1</th><td>AC1</td><td> 6609</td><td>2006</td><td>patient_1AC</td><td>AC</td><td> 1.119685</td><td>0.28160920</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGAGGGTAATT-1_AC1</th><td>AC1</td><td> 8229</td><td>2067</td><td>patient_1AC</td><td>AC</td><td> 2.150930</td><td>0.26580460</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGCATGTGTCA-1_AC1</th><td>AC1</td><td> 5266</td><td>1767</td><td>patient_1AC</td><td>AC</td><td> 2.449677</td><td>0.10632184</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGCATTGGCAT-1_AC1</th><td>AC1</td><td> 5093</td><td>1405</td><td>patient_1AC</td><td>AC</td><td> 2.611427</td><td>0.08189655</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGGTGAGTCAG-1_AC1</th><td>AC1</td><td> 5330</td><td>1551</td><td>patient_1AC</td><td>AC</td><td> 1.613508</td><td>0.09913793</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGGTGATACCT-1_AC1</th><td>AC1</td><td>15519</td><td>3210</td><td>patient_1AC</td><td>AC</td><td> 3.034989</td><td>0.30459770</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGGTTCGGCCA-1_AC1</th><td>AC1</td><td> 6841</td><td>1965</td><td>patient_1AC</td><td>AC</td><td> 2.061102</td><td>0.28160920</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGGTTTACTGG-1_AC1</th><td>AC1</td><td>  537</td><td> 216</td><td>patient_1AC</td><td>AC</td><td>10.986965</td><td>0.01149425</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGTCAGTAGGG-1_AC1</th><td>AC1</td><td> 4405</td><td>1654</td><td>patient_1AC</td><td>AC</td><td> 2.247446</td><td>0.16235632</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGTCGAATGCT-1_AC1</th><td>AC1</td><td> 1236</td><td> 560</td><td>patient_1AC</td><td>AC</td><td> 9.142395</td><td>0.01149425</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGTCTCGTGAA-1_AC1</th><td>AC1</td><td> 1544</td><td> 765</td><td>patient_1AC</td><td>AC</td><td> 3.432642</td><td>0.12787356</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTTGTCTTCTGTA-1_AC1</th><td>AC1</td><td> 2825</td><td>1143</td><td>patient_1AC</td><td>AC</td><td> 2.902655</td><td>0.10775862</td><td>Singlet</td></tr>
</tbody>
</table>




```R
for (i in 1:6) {
   colnames(doublet.list[[i]]@meta.data)[8] <-"DoubletFinder"
}
for (i in 1:6) {
   colnames(doublet.list[[i]]@meta.data)[7] <-"DoubletScore"
}
```


```R
#merge
 doublefinders_seurat <- merge(doublet.list[[1]], y=c(doublet.list[[2]],doublet.list[[3]],doublet.list[[4]],
                               doublet.list[[5]],doublet.list[[6]]))
```


```R
immune.combined[['DoubletFinder']] <- doublefinders_seurat@meta.data$DoubletFinder
immune.combined[['DoubletScore']] <- doublefinders_seurat@meta.data$DoubletScore
```


```R
immune.combined@meta.data
```


<table class="dataframe">
<caption>A data.frame: 48624 × 8</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_RNA</th><th scope=col>nFeature_RNA</th><th scope=col>sample</th><th scope=col>group</th><th scope=col>percent.mt</th><th scope=col>DoubletFinder</th><th scope=col>DoubletScore</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACCCAAGATTAGAC-1_AC1</th><td>AC1</td><td>11667</td><td>2905</td><td>patient_1AC</td><td>AC</td><td> 2.6056398</td><td>Doublet</td><td>0.36063218</td></tr>
	<tr><th scope=row>AAACCCAAGCATGTTC-1_AC1</th><td>AC1</td><td> 5385</td><td>1494</td><td>patient_1AC</td><td>AC</td><td> 3.8625812</td><td>Singlet</td><td>0.09482759</td></tr>
	<tr><th scope=row>AAACCCAAGCCTGTCG-1_AC1</th><td>AC1</td><td> 1471</td><td> 575</td><td>patient_1AC</td><td>AC</td><td> 7.3419443</td><td>Singlet</td><td>0.01149425</td></tr>
	<tr><th scope=row>AAACCCAAGGGTTTCT-1_AC1</th><td>AC1</td><td> 8957</td><td>2515</td><td>patient_1AC</td><td>AC</td><td> 1.8198057</td><td>Singlet</td><td>0.23850575</td></tr>
	<tr><th scope=row>AAACCCAAGTGCACCC-1_AC1</th><td>AC1</td><td>11681</td><td>2586</td><td>patient_1AC</td><td>AC</td><td> 1.8919613</td><td>Singlet</td><td>0.23994253</td></tr>
	<tr><th scope=row>AAACCCACAAGTTTGC-1_AC1</th><td>AC1</td><td> 2495</td><td> 794</td><td>patient_1AC</td><td>AC</td><td> 1.9639279</td><td>Singlet</td><td>0.02873563</td></tr>
	<tr><th scope=row>AAACCCACATTCTCCG-1_AC1</th><td>AC1</td><td> 8930</td><td>2326</td><td>patient_1AC</td><td>AC</td><td> 1.9148936</td><td>Singlet</td><td>0.18678161</td></tr>
	<tr><th scope=row>AAACCCAGTGTGTGTT-1_AC1</th><td>AC1</td><td>19264</td><td>3317</td><td>patient_1AC</td><td>AC</td><td> 1.9362542</td><td>Singlet</td><td>0.27729885</td></tr>
	<tr><th scope=row>AAACCCAGTGTTTCTT-1_AC1</th><td>AC1</td><td>  932</td><td> 476</td><td>patient_1AC</td><td>AC</td><td> 0.1072961</td><td>Singlet</td><td>0.10775862</td></tr>
	<tr><th scope=row>AAACCCAGTTACCCTC-1_AC1</th><td>AC1</td><td>10108</td><td>2366</td><td>patient_1AC</td><td>AC</td><td> 1.9885239</td><td>Singlet</td><td>0.18534483</td></tr>
	<tr><th scope=row>AAACCCATCTCACTCG-1_AC1</th><td>AC1</td><td>  809</td><td> 519</td><td>patient_1AC</td><td>AC</td><td> 6.9221261</td><td>Singlet</td><td>0.01293103</td></tr>
	<tr><th scope=row>AAACCCATCTCATTAC-1_AC1</th><td>AC1</td><td> 1798</td><td> 856</td><td>patient_1AC</td><td>AC</td><td> 3.0033370</td><td>Singlet</td><td>0.09482759</td></tr>
	<tr><th scope=row>AAACCCATCTCTCTAA-1_AC1</th><td>AC1</td><td>12381</td><td>2838</td><td>patient_1AC</td><td>AC</td><td> 1.8496083</td><td>Singlet</td><td>0.11637931</td></tr>
	<tr><th scope=row>AAACGAAAGAAGTGTT-1_AC1</th><td>AC1</td><td> 9902</td><td>2348</td><td>patient_1AC</td><td>AC</td><td> 3.7871137</td><td>Singlet</td><td>0.20114943</td></tr>
	<tr><th scope=row>AAACGAAAGAGTGAAG-1_AC1</th><td>AC1</td><td> 5238</td><td>1515</td><td>patient_1AC</td><td>AC</td><td> 5.7464681</td><td>Doublet</td><td>0.36781609</td></tr>
	<tr><th scope=row>AAACGAACACAGCCTG-1_AC1</th><td>AC1</td><td> 4594</td><td>1638</td><td>patient_1AC</td><td>AC</td><td> 5.6377884</td><td>Singlet</td><td>0.20977011</td></tr>
	<tr><th scope=row>AAACGAACACTGGACC-1_AC1</th><td>AC1</td><td> 5729</td><td>1743</td><td>patient_1AC</td><td>AC</td><td> 2.0596963</td><td>Singlet</td><td>0.31321839</td></tr>
	<tr><th scope=row>AAACGAAGTAGGCTGA-1_AC1</th><td>AC1</td><td> 3168</td><td>1172</td><td>patient_1AC</td><td>AC</td><td> 3.2196970</td><td>Singlet</td><td>0.10919540</td></tr>
	<tr><th scope=row>AAACGAAGTCGTTGCG-1_AC1</th><td>AC1</td><td>14693</td><td>2959</td><td>patient_1AC</td><td>AC</td><td> 0.3130743</td><td>Singlet</td><td>0.20258621</td></tr>
	<tr><th scope=row>AAACGAAGTTATGTGC-1_AC1</th><td>AC1</td><td> 8378</td><td>1997</td><td>patient_1AC</td><td>AC</td><td> 1.9813798</td><td>Singlet</td><td>0.24856322</td></tr>
	<tr><th scope=row>AAACGAAGTTCTCGTC-1_AC1</th><td>AC1</td><td> 3601</td><td>1445</td><td>patient_1AC</td><td>AC</td><td> 2.5826159</td><td>Singlet</td><td>0.31178161</td></tr>
	<tr><th scope=row>AAACGAAGTTTACCAG-1_AC1</th><td>AC1</td><td>  607</td><td> 386</td><td>patient_1AC</td><td>AC</td><td>10.3789127</td><td>Singlet</td><td>0.14224138</td></tr>
	<tr><th scope=row>AAACGAATCGATTCCC-1_AC1</th><td>AC1</td><td>15027</td><td>2654</td><td>patient_1AC</td><td>AC</td><td> 2.3091768</td><td>Singlet</td><td>0.31178161</td></tr>
	<tr><th scope=row>AAACGAATCTCTTGCG-1_AC1</th><td>AC1</td><td> 1437</td><td> 699</td><td>patient_1AC</td><td>AC</td><td> 3.4794711</td><td>Singlet</td><td>0.10344828</td></tr>
	<tr><th scope=row>AAACGCTAGGAGCAAA-1_AC1</th><td>AC1</td><td> 2794</td><td> 980</td><td>patient_1AC</td><td>AC</td><td> 8.6614173</td><td>Singlet</td><td>0.08764368</td></tr>
	<tr><th scope=row>AAACGCTAGGTTCCAT-1_AC1</th><td>AC1</td><td> 4126</td><td>1293</td><td>patient_1AC</td><td>AC</td><td> 2.1085797</td><td>Singlet</td><td>0.07040230</td></tr>
	<tr><th scope=row>AAACGCTAGTCCCGGT-1_AC1</th><td>AC1</td><td> 2692</td><td>1054</td><td>patient_1AC</td><td>AC</td><td> 4.3462110</td><td>Singlet</td><td>0.10344828</td></tr>
	<tr><th scope=row>AAACGCTAGTTACGTC-1_AC1</th><td>AC1</td><td>11864</td><td>2568</td><td>patient_1AC</td><td>AC</td><td> 1.6604855</td><td>Singlet</td><td>0.18390805</td></tr>
	<tr><th scope=row>AAACGCTCACAAGCTT-1_AC1</th><td>AC1</td><td> 8898</td><td>1994</td><td>patient_1AC</td><td>AC</td><td> 2.3038885</td><td>Singlet</td><td>0.18534483</td></tr>
	<tr><th scope=row>AAACGCTGTAGCTTGT-1_AC1</th><td>AC1</td><td> 2732</td><td>1303</td><td>patient_1AC</td><td>AC</td><td> 3.3308931</td><td>Singlet</td><td>0.22844828</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>TTTACTGGTCCTGTCT-1_PA3</th><td>PA3</td><td> 2368</td><td> 961</td><td>patient_3PA</td><td>PA</td><td> 3.7162162</td><td>Singlet</td><td>0.11290323</td></tr>
	<tr><th scope=row>TTTACTGGTCTTCTAT-1_PA3</th><td>PA3</td><td>10680</td><td>2985</td><td>patient_3PA</td><td>PA</td><td> 3.2677903</td><td>Singlet</td><td>0.08870968</td></tr>
	<tr><th scope=row>TTTACTGGTGATGTAA-1_PA3</th><td>PA3</td><td> 3002</td><td>1388</td><td>patient_3PA</td><td>PA</td><td> 5.6295803</td><td>Singlet</td><td>0.04032258</td></tr>
	<tr><th scope=row>TTTACTGGTTTGCCGG-1_PA3</th><td>PA3</td><td> 8028</td><td>2549</td><td>patient_3PA</td><td>PA</td><td> 3.5251619</td><td>Singlet</td><td>0.16935484</td></tr>
	<tr><th scope=row>TTTACTGTCCTGGGAC-1_PA3</th><td>PA3</td><td> 8915</td><td>2717</td><td>patient_3PA</td><td>PA</td><td> 2.9613012</td><td>Singlet</td><td>0.09677419</td></tr>
	<tr><th scope=row>TTTATGCAGTAGCATA-1_PA3</th><td>PA3</td><td> 1648</td><td> 508</td><td>patient_3PA</td><td>PA</td><td>20.3883495</td><td>Singlet</td><td>0.04838710</td></tr>
	<tr><th scope=row>TTTATGCGTGCATTTG-1_PA3</th><td>PA3</td><td>13077</td><td>3278</td><td>patient_3PA</td><td>PA</td><td> 3.5329204</td><td>Singlet</td><td>0.01612903</td></tr>
	<tr><th scope=row>TTTCACAAGTCGGCAA-1_PA3</th><td>PA3</td><td> 8185</td><td>2778</td><td>patient_3PA</td><td>PA</td><td> 2.8955406</td><td>Singlet</td><td>0.08064516</td></tr>
	<tr><th scope=row>TTTCACAAGTTTGAGA-1_PA3</th><td>PA3</td><td> 9928</td><td>2459</td><td>patient_3PA</td><td>PA</td><td> 4.1800967</td><td>Singlet</td><td>0.11290323</td></tr>
	<tr><th scope=row>TTTCAGTAGGATTTCC-1_PA3</th><td>PA3</td><td>10956</td><td>2808</td><td>patient_3PA</td><td>PA</td><td> 4.8010223</td><td>Singlet</td><td>0.21774194</td></tr>
	<tr><th scope=row>TTTCAGTCAAAGGCAC-1_PA3</th><td>PA3</td><td> 7957</td><td>1696</td><td>patient_3PA</td><td>PA</td><td> 1.9354028</td><td>Singlet</td><td>0.03225806</td></tr>
	<tr><th scope=row>TTTCATGCAACCAATC-1_PA3</th><td>PA3</td><td> 5742</td><td>1822</td><td>patient_3PA</td><td>PA</td><td> 4.4235458</td><td>Singlet</td><td>0.06451613</td></tr>
	<tr><th scope=row>TTTCATGCACAATGAA-1_PA3</th><td>PA3</td><td> 5551</td><td>1302</td><td>patient_3PA</td><td>PA</td><td> 3.4588362</td><td>Singlet</td><td>0.04032258</td></tr>
	<tr><th scope=row>TTTCATGGTATCGAGG-1_PA3</th><td>PA3</td><td> 7830</td><td>2388</td><td>patient_3PA</td><td>PA</td><td> 2.7969349</td><td>Singlet</td><td>0.04838710</td></tr>
	<tr><th scope=row>TTTCCTCAGCACTCTA-1_PA3</th><td>PA3</td><td> 6351</td><td>1867</td><td>patient_3PA</td><td>PA</td><td> 3.3695481</td><td>Singlet</td><td>0.06451613</td></tr>
	<tr><th scope=row>TTTCCTCCATTGAGCT-1_PA3</th><td>PA3</td><td> 9493</td><td>2809</td><td>patient_3PA</td><td>PA</td><td> 5.8885495</td><td>Singlet</td><td>0.08064516</td></tr>
	<tr><th scope=row>TTTCCTCTCAACTTTC-1_PA3</th><td>PA3</td><td> 7943</td><td>2072</td><td>patient_3PA</td><td>PA</td><td> 1.6492509</td><td>Singlet</td><td>0.09677419</td></tr>
	<tr><th scope=row>TTTCCTCTCAATCCGA-1_PA3</th><td>PA3</td><td> 6591</td><td>1815</td><td>patient_3PA</td><td>PA</td><td> 3.6565013</td><td>Singlet</td><td>0.13709677</td></tr>
	<tr><th scope=row>TTTCGATTCGATACTG-1_PA3</th><td>PA3</td><td> 8176</td><td>2213</td><td>patient_3PA</td><td>PA</td><td> 4.0851272</td><td>Doublet</td><td>0.26612903</td></tr>
	<tr><th scope=row>TTTCGATTCTGACAGT-1_PA3</th><td>PA3</td><td> 4963</td><td>1766</td><td>patient_3PA</td><td>PA</td><td> 3.8283296</td><td>Singlet</td><td>0.03225806</td></tr>
	<tr><th scope=row>TTTGACTAGTGATAGT-1_PA3</th><td>PA3</td><td>12615</td><td>3387</td><td>patient_3PA</td><td>PA</td><td> 3.0757035</td><td>Singlet</td><td>0.07258065</td></tr>
	<tr><th scope=row>TTTGACTCAATACGAA-1_PA3</th><td>PA3</td><td> 1270</td><td> 561</td><td>patient_3PA</td><td>PA</td><td>15.9842520</td><td>Singlet</td><td>0.05645161</td></tr>
	<tr><th scope=row>TTTGACTCATGTGACT-1_PA3</th><td>PA3</td><td> 7301</td><td>2249</td><td>patient_3PA</td><td>PA</td><td> 3.5337625</td><td>Singlet</td><td>0.07258065</td></tr>
	<tr><th scope=row>TTTGACTGTTCGAAGG-1_PA3</th><td>PA3</td><td>  614</td><td> 420</td><td>patient_3PA</td><td>PA</td><td>16.1237785</td><td>Singlet</td><td>0.01612903</td></tr>
	<tr><th scope=row>TTTGACTGTTCGGCTG-1_PA3</th><td>PA3</td><td>11523</td><td>3054</td><td>patient_3PA</td><td>PA</td><td> 3.0634383</td><td>Singlet</td><td>0.16129032</td></tr>
	<tr><th scope=row>TTTGATCTCGTCAACA-1_PA3</th><td>PA3</td><td> 8924</td><td>2582</td><td>patient_3PA</td><td>PA</td><td> 3.6418646</td><td>Singlet</td><td>0.08870968</td></tr>
	<tr><th scope=row>TTTGGAGGTCGCATGC-1_PA3</th><td>PA3</td><td> 5562</td><td>2054</td><td>patient_3PA</td><td>PA</td><td> 3.9374326</td><td>Singlet</td><td>0.12903226</td></tr>
	<tr><th scope=row>TTTGGAGTCATGGATC-1_PA3</th><td>PA3</td><td>11168</td><td>2589</td><td>patient_3PA</td><td>PA</td><td> 4.9874642</td><td>Singlet</td><td>0.24193548</td></tr>
	<tr><th scope=row>TTTGGTTCAAGAAATC-1_PA3</th><td>PA3</td><td>11099</td><td>2595</td><td>patient_3PA</td><td>PA</td><td> 3.5768988</td><td>Singlet</td><td>0.17741935</td></tr>
	<tr><th scope=row>TTTGTTGTCTTAGCAG-1_PA3</th><td>PA3</td><td> 4365</td><td>1379</td><td>patient_3PA</td><td>PA</td><td> 0.7560137</td><td>Singlet</td><td>0.18548387</td></tr>
</tbody>
</table>




```R
#
saveRDS(immune.combined, file ='Step1_AddDoublet_Seurat.Rds')
```
