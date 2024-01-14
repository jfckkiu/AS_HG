```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(DoubletFinder)
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
    



```R
main_path <- '~/AS/AS_Mouse/AS_Mouse2/'
```


```R
immune.combined <- readRDS(paste0(main_path,'Step2_QC_CCA_Seurat.Rds'))
```


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
    Positive:  Bgn, Timp3, Sparc, Igfbp7, Sod3, Nedd4, Fstl1, Fxyd1, Aebp1, Cpe 
    	   Ogn, Fhl1, Cald1, Rbp1, Cd200, Tm4sf1, Col1a2, Snhg18, Fbln5, Nfib 
    	   Map1b, Cped1, Eln, Pcdh7, Ptrf, Ptgis, Myl9, Rarres2, Serpinh1, Mgp 
    Negative:  Ftl1, Ctss, Lyz2, Cd74, Apoe, Lgmn, H2-Aa, H2-Ab1, H2-Eb1, Ier5 
    	   Cd14, Ly6e, H2-DMa, Cd68, Ifi27l2a, Srgn, Tgfbi, Slfn2, Fcgr1, H3f3b 
    	   Ms4a7, Ctsb, H2-DMb1, Ctsz, Pltp, Sepp1, AI607873, Wfdc17, Slc15a3, Ccl6 
    PC_ 2 
    Positive:  Cd52, Tmsb10, Lsp1, Actg1, Cytip, Btg1, Fxyd5, Bcl2a1d, Gngt2, Napsa 
    	   Plbd1, S100a11, Traf1, Rgs1, Il1b, Tspan13, Cd9, Lmnb1, Id2, Lgals3 
    	   Cd72, Gapdh, Runx3, Bhlhe40, Itgax, AA467197, Srgn, Atox1, Glipr1, Itgb7 
    Negative:  Cbr2, F13a1, Folr2, Pf4, Sepp1, Gas6, Mrc1, Maf, Ltc4s, Fcgrt 
    	   Igfbp4, Cfh, C4b, Lyve1, Stab1, Cd163, Rnase4, Lyz2, Timp2, Ninj1 
    	   Ccl24, Ccl6, Txnip, Apoe, Fxyd2, Hpgd, Ifi27l2a, Ly6e, Ccl8, Igf1 
    PC_ 3 
    Positive:  Cd14, Cd83, Fth1, Ifrd1, Nfkbia, Kdm6b, Gadd45b, Cxcl2, Ccrl2, Nfe2l2 
    	   Junb, Bcl2a1b, Csrnp1, Nfkbiz, Tlr2, Zfp36, Rab20, Dusp1, Slc15a3, Pim1 
    	   Nlrp3, Tgif1, Maff, Cdkn1a, Nr4a1, Tnf, Cxcl16, Skil, Cebpb, Atf3 
    Negative:  Cd3g, Cd3d, Nkg7, Ptprcap, Thy1, Lck, Cd3e, Ms4a4b, Gimap4, Ctsw 
    	   Skap1, Lat, Sept1, Gimap1, Gimap3, Cxcr6, Gimap7, Ltb, Cd8a, Cd2 
    	   Ccl5, Cd8b1, Gimap6, Ikzf3, Ctla2a, Gimap5, Tmsb10, Gzmk, Klrd1, Bcl11b 
    PC_ 4 
    Positive:  Napsa, Ltb, Klrd1, Itgb7, Nkg7, H2-Oa, H2-DMb2, Ccl5, Ctsw, Thy1 
    	   Cd3d, Cd3g, Cd7, Cd8a, Cd8b1, Dpp4, Gpr171, Cd244, Cxcr6, Cd3e 
    	   Ptprcap, Flt3, Gimap4, Mcemp1, Sept1, Ffar2, Tmsb10, Ms4a4b, Gimap3, Pkp3 
    Negative:  Birc5, Top2a, Nusap1, 2810417H13Rik, Stmn1, Ccna2, Mki67, Ube2c, Cdca3, Cdk1 
    	   Prc1, Tpx2, Cdca8, Hmmr, Kif22, Ckap2l, Pbk, Tk1, Fam64a, Kif11 
    	   Ccnb2, Rrm2, Cenpf, Aurkb, Neil3, Casc5, Cdca2, Racgap1, Spc24, Cenpe 
    PC_ 5 
    Positive:  Lgals3, Psap, Ctss, Ctsz, Atox1, Gm2a, Anxa1, Napsa, Ifi30, Anxa2 
    	   Itgax, Cd9, Cd68, Cst3, Prdx1, Capg, Anxa5, Pkib, Gusb, Syngr1 
    	   Rnh1, Gngt2, Msrb1, Anpep, Gpr137b, Fgr, Gpnmb, Cd300a, Abcg1, Mpeg1 
    Negative:  Cd3g, Nkg7, Cd3d, Thy1, Lck, Cd3e, Ctsw, Gimap4, Skap1, Ptprcap 
    	   Lat, Ms4a4b, Gimap3, Cxcr6, Cd8a, Gimap7, Gimap1, Cd8b1, Ccl5, Sept1 
    	   Cd2, Gzmk, Gimap5, Ctla2a, Ikzf3, Tnfaip3, Cd247, Il2rb, Bcl11b, Nfkbia 
    
    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    14:34:16 UMAP embedding parameters a = 0.9922 b = 1.112
    
    14:34:16 Read 3572 rows and found 20 numeric columns
    
    14:34:16 Using Annoy for neighbor search, n_neighbors = 30
    
    14:34:16 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    14:34:17 Writing NN index file to temp file /tmp/RtmpRb8Ew0/file5575123c395
    
    14:34:17 Searching Annoy index using 1 thread, search_k = 3000
    
    14:34:17 Annoy recall = 100%
    
    14:34:17 Commencing smooth kNN distance calibration using 1 thread
    
    14:34:18 Initializing from normalized Laplacian + noise
    
    14:34:18 Commencing optimization for 500 epochs, with 143966 positive edges
    
    14:34:20 Optimization finished
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  Pf4, F13a1, Mgl2, Folr2, Ly6e, Pltp, Ifi27l2a, Cd14, Ltc4s, Marcksl1 
    	   Zfp36, Ier5, Lyve1, Nfkbia, Lcp1, AI607873, Wfdc17, Ccl24, Ccl2, Ier3 
    	   Ninj1, Hmox1, H3f3b, Junb, Ccl9, Ms4a7, Mafb, Cd163, Cd74, Ubc 
    Negative:  Igfbp7, Bgn, Sod3, Fstl1, Tm4sf1, Cpe, Aebp1, Lhfp, Fxyd1, Eln 
    	   Fhl1, Cald1, Sparc, Serpinh1, Ogn, Timp3, Col3a1, Htra1, Nedd4, Cd200 
    	   Col1a2, Snhg18, Rarres2, Nfib, Col4a1, Ptgis, Rbp1, Fmo2, Ltbp4, Prkcdbp 
    PC_ 2 
    Positive:  Folr2, Ninj1, F13a1, Pf4, Fcgrt, Clec10a, Lyve1, Ifitm3, Ccl24, Ltc4s 
    	   Pltp, Wfdc17, Tslp, Hmox1, Ccl6, Marcksl1, Cd209f, Mgl2, Ctsl, Mt1 
    	   Mmp9, Ednrb, Cd209g, Fxyd2, Pla2g2d, Jun, Cd163, Timd4, Siglech, C4b 
    Negative:  Nusap1, Top2a, Ccna2, Birc5, Cdca3, Prc1, Ube2c, 2810417H13Rik, Hmmr, Lockd 
    	   Stmn1, Tpx2, Cdk1, Mki67, Ccnb2, Cks1b, Cdca8, Racgap1, Pbk, Fam64a 
    	   Casc5, Kif22, Tk1, Lmnb1, Ncapd2, Ndc80, Ckap2l, Kif11, Cenpf, Nuf2 
    PC_ 3 
    Positive:  Nfkbia, Junb, Zfp36, Ier5, Nfkbiz, Cxcl2, Atf3, H3f3b, Socs3, Ifrd1 
    	   Ccl4, Btg2, Cd83, Kdm6b, Gadd45b, Pim1, Jun, Skil, Ubc, Ier3 
    	   Egr1, Tnf, Ccrl2, Ccl3, Cdkn1a, Mt1, Gadd45g, Gm17056, Slc15a3, Csrnp1 
    Negative:  Ccl6, Folr2, Fcna, Crip1, Ccl9, Fcgrt, Pltp, S100a10, Lyve1, F13a1 
    	   Ednrb, Fxyd2, Ninj1, Lgals1, Ltc4s, Wfdc17, Ccl24, Fgfr1, Emp1, Timd4 
    	   Rcn3, Cfp, S100a4, Pf4, C4b, Ecm1, Clec10a, S100a6, Ifitm3, Tppp3 
    PC_ 4 
    Positive:  Lsp1, Napsa, Lgals3, Cd52, Cytip, Itgax, Cadm1, Cd24a, Cd9, Mmp12 
    	   Bcl2a1d, Ccr2, Mpeg1, Plbd1, Tmsb10, Gngt2, Bcl2a1a, Itgb7, Ltb4r1, Sh2d1b1 
    	   Coro1a, Gm2a, Ear2, Tspan13, Abi3, Alcam, Fxyd5, H2-Oa, St3gal5, H2-DMb2 
    Negative:  Pf4, F13a1, Ier3, Mt1, Ly6e, Folr2, Jun, Socs3, Ltc4s, Junb 
    	   Ccl2, Zfp36, Ccl24, Hmox1, Ninj1, Marcksl1, Egr1, Fcgrt, AI607873, Cd163 
    	   Rnasel, Lyve1, Mt2, Ubc, Zfand5, Jund, Pla2g2d, Casp4, Prc1, Nusap1 
    PC_ 5 
    Positive:  Dcn, Abca8a, Lum, Cygb, Pi16, F3, Serping1, Dpt, Pdgfra, Fbln1 
    	   Bicc1, Cd34, Pdgfrl, Dpep1, Mmp2, Slc43a3, Ctsk, Smoc2, Pcolce, Clec11a 
    	   Hsd11b1, Aspn, Entpd2, C1s1, Tspan11, Ndufa4l2, Mdk, Gas1, Nbl1, Abi3bp 
    Negative:  Cnn1, Myh11, Postn, Npnt, Col18a1, Susd5, Tagln, Sost, Pdlim3, Mcam 
    	   Synpo2, Ncam1, Sh3bgr, Enah, Kcnmb1, Acta2, Itga8, Pcp4l1, Scube3, Sgcg 
    	   Efhd1, Cspg4, Bcam, Tpm2, Lmcd1, Tspan2, Rasl12, Tinagl1, Ppp1r14a, Nexn 
    
    14:34:25 UMAP embedding parameters a = 0.9922 b = 1.112
    
    14:34:25 Read 6208 rows and found 20 numeric columns
    
    14:34:25 Using Annoy for neighbor search, n_neighbors = 30
    
    14:34:25 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    14:34:25 Writing NN index file to temp file /tmp/RtmpRb8Ew0/file5571a2b8da5
    
    14:34:25 Searching Annoy index using 1 thread, search_k = 3000
    
    14:34:26 Annoy recall = 100%
    
    14:34:26 Commencing smooth kNN distance calibration using 1 thread
    
    14:34:27 Initializing from normalized Laplacian + noise
    
    14:34:27 Commencing optimization for 500 epochs, with 255566 positive edges
    
    14:34:31 Optimization finished
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  Lyz2, Ctss, Pf4, F13a1, Ctsb, Ftl1, Apoe, Folr2, Cbr2, Mgl2 
    	   Mafb, Actb, Ccl9, Cyba, Aif1, Ifi27l2a, Ms4a7, Lyve1, Ccl6, Fcrls 
    	   Wfdc17, Ltc4s, Pltp, Itgam, H3f3b, Cd163, Cd68, Ifi203, Plek, Ccl24 
    Negative:  Bgn, Nfib, Fstl1, Timp3, Cpe, Eln, Aebp1, Sparc, Ccdc80, Nedd4 
    	   Cd200, Pam, Tm4sf1, Serpinh1, Plpp3, Rbp1, Ltbp4, Gxylt2, Sod3, Id3 
    	   Meg3, Mfap5, Fhl1, Mgp, Col5a2, Col4a1, Cald1, Gata6, Lox, Emp2 
    PC_ 2 
    Positive:  Mmrn2, Pecam1, Cytl1, Gja5, Bmx, Cdh5, Sdpr, Myct1, Esam, Vwf 
    	   Nos3, Tek, Epas1, Podxl, Clu, Cldn5, Ecscr, Egfl7, Palmd, Slco2a1 
    	   Tinagl1, Tmem158, Clec14a, Rasip1, Tie1, Edn1, Heg1, Car8, Procr, Efnb2 
    Negative:  Ogn, Col1a1, Smoc2, Adamts2, Abca8a, F3, Dpt, Lum, Cygb, Svep1 
    	   Col6a1, Loxl1, Bicc1, Lama2, Fbln1, Pdgfrl, Mfap4, Col3a1, Prrx1, Scn7a 
    	   Gdf10, Fmod, Podn, Pdgfra, C1s1, Abcc9, Col1a2, Wnt5a, Dpep1, Nr2f2 
    PC_ 3 
    Positive:  Cst3, Ftl1, Pf4, Apoe, Ctsb, Cfh, F13a1, Cbr2, Lyz2, Ifitm3 
    	   Ninj1, Folr2, Lyve1, Ctss, Mafb, Pltp, Ltc4s, Cd163, Cd63, Ccl24 
    	   Ccl9, Mgl2, Fcrls, Dstn, Ifi27l2a, Fth1, Ccl6, Klf2, Fxyd2, Ms4a7 
    Negative:  Ptprcap, Cd3g, Sept1, Cxcr6, Il7r, Bcl11b, Ltb, Cd3e, Cd3d, Il18r1 
    	   Lat, Itgb7, Thy1, Podnl1, Actn2, Ets1, Cd163l1, Blk, Skap1, Tmsb10 
    	   Ikzf3, Gimap3, Il17re, Ltb4r1, Cd247, Ramp1, Tbc1d10c, Ly6g5b, Rgcc, Cd52 
    PC_ 4 
    Positive:  Cnn1, Myh11, Tagln, Itga8, Npnt, Myl9, Acta2, Serpine2, Sost, Col18a1 
    	   Susd5, Tpm2, Synpo2, Postn, Mcam, Sh3bgr, Nov, Thsd4, Pcdh7, Lmod1 
    	   Gucy1a3, Dsp, Kcnmb1, Efhd1, Ppp1r14a, Ptprz1, Mylk4, 1190005I06Rik, Adcy5, Scube3 
    Negative:  Smoc2, Abca8a, Cygb, Svep1, Lum, Dcn, Dpt, F3, Scn7a, Fbln1 
    	   Dpep1, Pdgfra, Bicc1, Gdf10, Cpxm1, Pdgfrl, Fmod, C1s1, Serping1, Mmp2 
    	   Pcolce, Abca8b, Podn, Tshz2, Inmt, Cxcl12, Wnt5a, Mdk, Srpx, Hsd11b1 
    PC_ 5 
    Positive:  Ltb, Il18r1, Il7r, Cxcr6, Ptprcap, Itgb7, Actn2, Ly6g5b, Ltb4r1, Cd3g 
    	   Podnl1, Sept1, Cd3d, Il17re, Cd163l1, Blk, Gimap3, Cd3e, Ikzf3, Lat 
    	   Ramp1, Pkp3, F2r, Il23r, Zap70, Il27ra, Aqp3, Il1r1, Gimap1, Amica1 
    Negative:  Nusap1, Top2a, Mki67, Prc1, Cdca3, Casc5, Birc5, Pbk, Tpx2, Hmmr 
    	   Cenpf, Stmn1, Hist1h1b, Cenpe, Kif11, Esco2, 2810417H13Rik, Ckap2, Rad51ap1, Sgol2a 
    	   Ckap2l, Kif15, Ccnb2, Ube2c, Sgol1, Neil3, Aspm, Cdca8, Lockd, Cdca2 
    
    14:34:33 UMAP embedding parameters a = 0.9922 b = 1.112
    
    14:34:33 Read 2203 rows and found 20 numeric columns
    
    14:34:33 Using Annoy for neighbor search, n_neighbors = 30
    
    14:34:33 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    14:34:33 Writing NN index file to temp file /tmp/RtmpRb8Ew0/file557473b67d4
    
    14:34:33 Searching Annoy index using 1 thread, search_k = 3000
    
    14:34:33 Annoy recall = 100%
    
    14:34:33 Commencing smooth kNN distance calibration using 1 thread
    
    14:34:34 Initializing from normalized Laplacian + noise
    
    14:34:34 Commencing optimization for 500 epochs, with 85262 positive edges
    
    14:34:35 Optimization finished
    



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
    
    Loading required package: dotCall64
    
    Loading required package: grid
    
    Spam version 2.7-0 (2021-06-25) is loaded.
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
    



    
![png](Step7.2._Doublet_Finder_files/Step7.2._Doublet_Finder_6_37.png)
    


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



    
![png](Step7.2._Doublet_Finder_files/Step7.2._Doublet_Finder_6_39.png)
    


    NULL



<dl>
	<dt>$GSE116240_12W</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 31 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.005</td><td>0.7936991</td><td>0.0024336344</td><td> 326.1373</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.7960958</td><td>0.0032242930</td><td> 246.9055</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.7982389</td><td>0.0031092982</td><td> 256.7264</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.8073271</td><td>0.0003535185</td><td>2283.6915</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.8028922</td><td>0.0002042705</td><td>3930.5342</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.7891753</td><td>0.0004115064</td><td>1917.7714</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.7636648</td><td>0.0006516029</td><td>1171.9788</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.7517697</td><td>0.0006775960</td><td>1109.4660</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.7375980</td><td>0.0009599494</td><td> 768.3717</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.7358395</td><td>0.0017188162</td><td> 428.1083</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.7288297</td><td>0.0021249511</td><td> 342.9866</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.7239551</td><td>0.0029746836</td><td> 243.3721</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.7170198</td><td>0.0045277206</td><td> 158.3622</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.7092186</td><td>0.0056434447</td><td> 125.6712</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.7032444</td><td>0.0049109931</td><td> 143.1980</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.6956643</td><td>0.0046106916</td><td> 150.8807</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.6921461</td><td>0.0037456826</td><td> 184.7850</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.6948212</td><td>0.0035691661</td><td> 194.6733</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.6967897</td><td>0.0035329418</td><td> 197.2265</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.6981054</td><td>0.0033447429</td><td> 208.7172</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.6996407</td><td>0.0032327414</td><td> 216.4233</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.7015083</td><td>0.0032796899</td><td> 213.8947</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.7057069</td><td>0.0025654275</td><td> 275.0836</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.7083597</td><td>0.0020695453</td><td> 342.2779</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.7148281</td><td>0.0017371902</td><td> 411.4852</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.7175904</td><td>0.0011431947</td><td> 627.7062</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.7174324</td><td>0.0008145021</td><td> 880.8232</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.7187694</td><td>0.0009691388</td><td> 741.6579</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.7200258</td><td>0.0010421603</td><td> 690.8973</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.7159845</td><td>0.0015150875</td><td> 472.5698</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.7140345</td><td>0.0020738678</td><td> 344.3009</td></tr>
</tbody>
</table>
</dd>
	<dt>$GSE154817_C</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 31 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.005</td><td>0.7529766</td><td>0.0032198039</td><td> 233.8579</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.7549381</td><td>0.0019603637</td><td> 385.1010</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.7518757</td><td>0.0010230475</td><td> 734.9372</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.7679656</td><td>0.0020840557</td><td> 368.4957</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.7606316</td><td>0.0012175672</td><td> 624.7143</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.7452089</td><td>0.0008816118</td><td> 845.2800</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.7313155</td><td>0.0008750443</td><td> 835.7469</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.7204981</td><td>0.0007617549</td><td> 945.8398</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.7118188</td><td>0.0007733575</td><td> 920.4267</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.7050692</td><td>0.0006562146</td><td>1074.4492</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.6975767</td><td>0.0005055337</td><td>1379.8817</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.6920384</td><td>0.0003983642</td><td>1737.2003</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.6880263</td><td>0.0004174763</td><td>1648.0609</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.6865138</td><td>0.0005928327</td><td>1158.0229</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.6850830</td><td>0.0007372067</td><td> 929.2957</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.6804553</td><td>0.0004678355</td><td>1454.4757</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.6752684</td><td>0.0003516128</td><td>1920.4886</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.6723439</td><td>0.0002243355</td><td>2997.0458</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.6699126</td><td>0.0001163275</td><td>5758.8494</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.6698819</td><td>0.0001388066</td><td>4826.0088</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.6676083</td><td>0.0001649459</td><td>4047.4390</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.6684882</td><td>0.0002359379</td><td>2833.3228</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.6687935</td><td>0.0002962084</td><td>2257.8480</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.6693256</td><td>0.0002941049</td><td>2275.8058</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.6684238</td><td>0.0002727823</td><td>2450.3928</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.6687287</td><td>0.0002631249</td><td>2541.4875</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.6719776</td><td>0.0003058507</td><td>2197.0768</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.6727506</td><td>0.0002838598</td><td>2370.0105</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.6769447</td><td>0.0002891960</td><td>2340.7813</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.6794631</td><td>0.0002383631</td><td>2850.5381</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.6821598</td><td>0.0001754819</td><td>3887.3507</td></tr>
</tbody>
</table>
</dd>
	<dt>$GSE154817_3W</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 31 × 5</caption>
<thead>
	<tr><th scope=col>ParamID</th><th scope=col>pK</th><th scope=col>MeanBC</th><th scope=col>VarBC</th><th scope=col>BCmetric</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 1</td><td>0.005</td><td>0.7414644</td><td>1.155075e-03</td><td>  641.91897</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.7756462</td><td>3.009956e-04</td><td> 2576.93487</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.7985325</td><td>1.758297e-03</td><td>  454.15111</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.7547415</td><td>2.508156e-03</td><td>  300.91494</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.7024970</td><td>1.441040e-02</td><td>   48.74930</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.7099211</td><td>1.121019e-02</td><td>   63.32818</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.6939579</td><td>5.859282e-03</td><td>  118.43736</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.6891059</td><td>2.713212e-03</td><td>  253.98159</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.6733169</td><td>2.860550e-03</td><td>  235.38024</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.6652440</td><td>3.223636e-03</td><td>  206.36448</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.6516442</td><td>2.951993e-03</td><td>  220.74719</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.6440159</td><td>2.824000e-03</td><td>  228.05091</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.6346132</td><td>2.639214e-03</td><td>  240.45539</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.6318293</td><td>1.683998e-03</td><td>  375.19607</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.6201970</td><td>1.097624e-03</td><td>  565.03570</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.6163714</td><td>8.105115e-04</td><td>  760.47208</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.6166468</td><td>7.279971e-04</td><td>  847.04568</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.6122884</td><td>6.488248e-04</td><td>  943.68841</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.6097996</td><td>3.031798e-04</td><td> 2011.34630</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.6051726</td><td>5.948409e-05</td><td>10173.68819</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.6066827</td><td>5.952745e-05</td><td>10191.64733</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.6123026</td><td>9.513609e-05</td><td> 6436.07121</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.6133412</td><td>5.510014e-05</td><td>11131.39055</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.6174696</td><td>1.054562e-04</td><td> 5855.22089</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.6180139</td><td>1.300726e-04</td><td> 4751.30099</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.6198982</td><td>2.609443e-04</td><td> 2375.59630</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.6181564</td><td>3.165545e-04</td><td> 1952.76423</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.6192165</td><td>3.626471e-04</td><td> 1707.49048</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.6188176</td><td>2.858874e-04</td><td> 2164.55009</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.6188099</td><td>4.237761e-04</td><td> 1460.22845</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.6216962</td><td>4.085279e-04</td><td> 1521.79642</td></tr>
</tbody>
</table>
</dd>
</dl>




    
![png](Step7.2._Doublet_Finder_files/Step7.2._Doublet_Finder_6_42.png)
    



```R
#find optimalPK
opt_pk_list <- lapply (X = PK.list, FUN = function(x) { 
    as.numeric(as.vector(x$pK[which.max(x$BCmetric)]))
})
opt_pk_list
```


<dl>
	<dt>$GSE116240_12W</dt>
		<dd>0.04</dd>
	<dt>$GSE154817_C</dt>
		<dd>0.18</dd>
	<dt>$GSE154817_3W</dt>
		<dd>0.22</dd>
</dl>




```R
#find optimal_pk for 1:9 by artifical
i = 1
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
	<tr><td> 1</td><td>0.005</td><td>0.7936991</td><td>0.0024336344</td><td> 326.1373</td></tr>
	<tr><td> 2</td><td>0.01 </td><td>0.7960958</td><td>0.0032242930</td><td> 246.9055</td></tr>
	<tr><td> 3</td><td>0.02 </td><td>0.7982389</td><td>0.0031092982</td><td> 256.7264</td></tr>
	<tr><td> 4</td><td>0.03 </td><td>0.8073271</td><td>0.0003535185</td><td>2283.6915</td></tr>
	<tr><td> 5</td><td>0.04 </td><td>0.8028922</td><td>0.0002042705</td><td>3930.5342</td></tr>
	<tr><td> 6</td><td>0.05 </td><td>0.7891753</td><td>0.0004115064</td><td>1917.7714</td></tr>
	<tr><td> 7</td><td>0.06 </td><td>0.7636648</td><td>0.0006516029</td><td>1171.9788</td></tr>
	<tr><td> 8</td><td>0.07 </td><td>0.7517697</td><td>0.0006775960</td><td>1109.4660</td></tr>
	<tr><td> 9</td><td>0.08 </td><td>0.7375980</td><td>0.0009599494</td><td> 768.3717</td></tr>
	<tr><td>10</td><td>0.09 </td><td>0.7358395</td><td>0.0017188162</td><td> 428.1083</td></tr>
	<tr><td>11</td><td>0.1  </td><td>0.7288297</td><td>0.0021249511</td><td> 342.9866</td></tr>
	<tr><td>12</td><td>0.11 </td><td>0.7239551</td><td>0.0029746836</td><td> 243.3721</td></tr>
	<tr><td>13</td><td>0.12 </td><td>0.7170198</td><td>0.0045277206</td><td> 158.3622</td></tr>
	<tr><td>14</td><td>0.13 </td><td>0.7092186</td><td>0.0056434447</td><td> 125.6712</td></tr>
	<tr><td>15</td><td>0.14 </td><td>0.7032444</td><td>0.0049109931</td><td> 143.1980</td></tr>
	<tr><td>16</td><td>0.15 </td><td>0.6956643</td><td>0.0046106916</td><td> 150.8807</td></tr>
	<tr><td>17</td><td>0.16 </td><td>0.6921461</td><td>0.0037456826</td><td> 184.7850</td></tr>
	<tr><td>18</td><td>0.17 </td><td>0.6948212</td><td>0.0035691661</td><td> 194.6733</td></tr>
	<tr><td>19</td><td>0.18 </td><td>0.6967897</td><td>0.0035329418</td><td> 197.2265</td></tr>
	<tr><td>20</td><td>0.19 </td><td>0.6981054</td><td>0.0033447429</td><td> 208.7172</td></tr>
	<tr><td>21</td><td>0.2  </td><td>0.6996407</td><td>0.0032327414</td><td> 216.4233</td></tr>
	<tr><td>22</td><td>0.21 </td><td>0.7015083</td><td>0.0032796899</td><td> 213.8947</td></tr>
	<tr><td>23</td><td>0.22 </td><td>0.7057069</td><td>0.0025654275</td><td> 275.0836</td></tr>
	<tr><td>24</td><td>0.23 </td><td>0.7083597</td><td>0.0020695453</td><td> 342.2779</td></tr>
	<tr><td>25</td><td>0.24 </td><td>0.7148281</td><td>0.0017371902</td><td> 411.4852</td></tr>
	<tr><td>26</td><td>0.25 </td><td>0.7175904</td><td>0.0011431947</td><td> 627.7062</td></tr>
	<tr><td>27</td><td>0.26 </td><td>0.7174324</td><td>0.0008145021</td><td> 880.8232</td></tr>
	<tr><td>28</td><td>0.27 </td><td>0.7187694</td><td>0.0009691388</td><td> 741.6579</td></tr>
	<tr><td>29</td><td>0.28 </td><td>0.7200258</td><td>0.0010421603</td><td> 690.8973</td></tr>
	<tr><td>30</td><td>0.29 </td><td>0.7159845</td><td>0.0015150875</td><td> 472.5698</td></tr>
	<tr><td>31</td><td>0.3  </td><td>0.7140345</td><td>0.0020738678</td><td> 344.3009</td></tr>
</tbody>
</table>




    
![png](Step7.2._Doublet_Finder_files/Step7.2._Doublet_Finder_8_1.png)
    



```R
opt_pk_list[[1]] <-0.04
opt_pk_list[[2]] <-0.02
opt_pk_list[[3]] <-0.01
```


```R
homotypic.list <- lapply(X = doublet.list, FUN = function(x) {
        x@meta.data$celltype %>%
        modelHomotypic()})
homotypic.list 
```


<dl>
	<dt>$GSE116240_12W</dt>
		<dd>0</dd>
	<dt>$GSE154817_C</dt>
		<dd>0</dd>
	<dt>$GSE154817_3W</dt>
		<dd>0</dd>
</dl>




```R
doublet.list[[1]]@assays
```


    $RNA
    Assay data with 16581 features for 3572 cells
    Top 10 variable features:
     Retnla, Fabp5, Ccl5, Spp1, Gpnmb, Ifitm1, Cxcl1, Plac8, Cxcl13, Acta2 
    
    $integrated
    Assay data with 3000 features for 3572 cells
    Top 10 variable features:
     Fabp5, Ccl5, Igfbp7, Retnla, Sparc, Gpnmb, Tagln, Acta2, Col1a2, Bgn 




```R
#doubletrate-0.05
nExp_poi.list <- lapply(X = doublet.list, FUN = function(x) {
    round(0.05*ncol(x@assays$integrated@data))
    })
```


```R
nExp_poi.adj.list <- nExp_poi.list
for (i in 1:3) {
  nExp_poi.adj.list[[i]] <- round(nExp_poi.list[[i]]*(1-homotypic.list[[i]]))   
}
```


```R
Exp_poi.adj.num <- round(as.numeric(nExp_poi.list)*(1-as.numeric(homotypic.list)))
```


```R
for (i in 1:3) {
    doublet.list[[i]] <- doubletFinder_v3(doublet.list[[i]],
    PCs = 1:20,
    pN = 0.25,
    pK = as.numeric(opt_pk_list[[i]]),
    nExp = as.numeric(nExp_poi.adj.list[[i]]),
    reuse.pANN = FALSE,
    sct = FALSE                                         
    )}
```

    [1] "Creating 1191 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."
    [1] "Creating 2069 artificial doublets..."
    [1] "Creating Seurat object..."
    [1] "Normalizing Seurat object..."
    [1] "Finding variable genes..."
    [1] "Scaling data..."


    Centering and scaling data matrix
    


    [1] "Running PCA..."
    [1] "Calculating PC distance matrix..."
    [1] "Computing pANN..."
    [1] "Classifying doublets.."
    [1] "Creating 734 artificial doublets..."
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
<caption>A data.frame: 3572 × 9</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_RNA</th><th scope=col>nFeature_RNA</th><th scope=col>sample</th><th scope=col>group</th><th scope=col>Dataset</th><th scope=col>percent.mt</th><th scope=col>pANN_0.25_0.04_179</th><th scope=col>DF.classifications_0.25_0.04_179</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACCTGAGATGCCTT-1_GSE116240_12W</th><td>GSE116240</td><td> 8162</td><td>2217</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.3278608</td><td>0.26701571</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGAGCTAGTGG-1_GSE116240_12W</th><td>GSE116240</td><td> 3117</td><td>1209</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>3.4648701</td><td>0.35078534</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGCAAGCTGTT-1_GSE116240_12W</th><td>GSE116240</td><td> 8290</td><td>2084</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.4716526</td><td>0.14136126</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGCACCGAATT-1_GSE116240_12W</th><td>GSE116240</td><td> 5460</td><td>1961</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.7032967</td><td>0.21989529</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGCACGCTTTC-1_GSE116240_12W</th><td>GSE116240</td><td> 2434</td><td> 908</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>0.3286771</td><td>0.13089005</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGCAGGACCCT-1_GSE116240_12W</th><td>GSE116240</td><td>11581</td><td>2985</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.0205509</td><td>0.43455497</td><td>Doublet</td></tr>
	<tr><th scope=row>AAACCTGGTCACAAGG-1_GSE116240_12W</th><td>GSE116240</td><td> 6385</td><td>2028</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>3.1480031</td><td>0.18324607</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGGTCACCTAA-1_GSE116240_12W</th><td>GSE116240</td><td> 4679</td><td>1727</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.5387903</td><td>0.04188482</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGGTCCGTGAC-1_GSE116240_12W</th><td>GSE116240</td><td>12945</td><td>2997</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.6763229</td><td>0.22513089</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGGTCTCACCT-1_GSE116240_12W</th><td>GSE116240</td><td> 2832</td><td>1154</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.3658192</td><td>0.32984293</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGTCGCCTGAG-1_GSE116240_12W</th><td>GSE116240</td><td> 4035</td><td>1434</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.2057001</td><td>0.07329843</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACCTGTCGTAGGTT-1_GSE116240_12W</th><td>GSE116240</td><td> 5497</td><td>1833</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>3.4018556</td><td>0.21989529</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGGGAGAAGATTC-1_GSE116240_12W</th><td>GSE116240</td><td> 5753</td><td>2014</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>0.9386407</td><td>0.31413613</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGGGAGTGCCATT-1_GSE116240_12W</th><td>GSE116240</td><td> 7279</td><td>2205</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.4974584</td><td>0.08376963</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGGGCACTATCTT-1_GSE116240_12W</th><td>GSE116240</td><td> 5947</td><td>1463</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.7992265</td><td>0.32460733</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGGGGTGTGCCTG-1_GSE116240_12W</th><td>GSE116240</td><td>11091</td><td>2345</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.9835903</td><td>0.30890052</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGGGTCATGCATG-1_GSE116240_12W</th><td>GSE116240</td><td> 8783</td><td>2297</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.5598315</td><td>0.33507853</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGGGTCGTGGACC-1_GSE116240_12W</th><td>GSE116240</td><td>10611</td><td>2514</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.2806522</td><td>0.31937173</td><td>Singlet</td></tr>
	<tr><th scope=row>AAACGGGTCTTCATGT-1_GSE116240_12W</th><td>GSE116240</td><td>13222</td><td>2604</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.4807140</td><td>0.26701571</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGAGAGGTACC-1_GSE116240_12W</th><td>GSE116240</td><td> 5526</td><td>1956</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.4791893</td><td>0.26701571</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGAGCGTGAAC-1_GSE116240_12W</th><td>GSE116240</td><td> 6537</td><td>2125</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.2946306</td><td>0.22513089</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGAGTGATCGG-1_GSE116240_12W</th><td>GSE116240</td><td>10948</td><td>2636</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>4.5213738</td><td>0.14136126</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGCATGAAGTA-1_GSE116240_12W</th><td>GSE116240</td><td>14402</td><td>3041</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.6317178</td><td>0.30890052</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGGTCTCTTAT-1_GSE116240_12W</th><td>GSE116240</td><td>12378</td><td>3033</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.3814833</td><td>0.40314136</td><td>Doublet</td></tr>
	<tr><th scope=row>AAAGATGGTGAAAGAG-1_GSE116240_12W</th><td>GSE116240</td><td>17873</td><td>3544</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>0.9455603</td><td>0.27225131</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGGTGTGACGA-1_GSE116240_12W</th><td>GSE116240</td><td>14029</td><td>3088</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.2738613</td><td>0.21989529</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGTCACATACG-1_GSE116240_12W</th><td>GSE116240</td><td> 8335</td><td>2280</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.2315537</td><td>0.17277487</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGATGTCGGAGGTA-1_GSE116240_12W</th><td>GSE116240</td><td> 4031</td><td>1697</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>2.7784669</td><td>0.06806283</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGCAAAGCTCCTTC-1_GSE116240_12W</th><td>GSE116240</td><td> 9605</td><td>2152</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>3.3732431</td><td>0.13612565</td><td>Singlet</td></tr>
	<tr><th scope=row>AAAGCAAAGGAACTGC-1_GSE116240_12W</th><td>GSE116240</td><td> 9817</td><td>2354</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>1.7724356</td><td>0.25130890</td><td>Singlet</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>TTTATGCTCGATGAGG-1_GSE116240_12W</th><td>GSE116240</td><td> 3305</td><td>1636</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>11.921331</td><td>0.29319372</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTATGCTCTAACCGA-1_GSE116240_12W</th><td>GSE116240</td><td> 7825</td><td>2459</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.750799</td><td>0.14136126</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTCCTCAGATGGGTC-1_GSE116240_12W</th><td>GSE116240</td><td> 8675</td><td>2405</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.201729</td><td>0.47643979</td><td>Doublet</td></tr>
	<tr><th scope=row>TTTCCTCCAAGTAGTA-1_GSE116240_12W</th><td>GSE116240</td><td>11876</td><td>3412</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 4.942742</td><td>0.47120419</td><td>Doublet</td></tr>
	<tr><th scope=row>TTTCCTCCACGGCCAT-1_GSE116240_12W</th><td>GSE116240</td><td>11557</td><td>3017</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.600761</td><td>0.21465969</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTCCTCCAGGAATGC-1_GSE116240_12W</th><td>GSE116240</td><td> 7896</td><td>2090</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.482270</td><td>0.14659686</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTCCTCGTTACGACT-1_GSE116240_12W</th><td>GSE116240</td><td>13259</td><td>3119</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 3.197828</td><td>0.20942408</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTCCTCGTTCCTCCA-1_GSE116240_12W</th><td>GSE116240</td><td> 8768</td><td>2699</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.474909</td><td>0.08900524</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTCCTCTCCTCATTA-1_GSE116240_12W</th><td>GSE116240</td><td>13007</td><td>3063</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.045053</td><td>0.49214660</td><td>Doublet</td></tr>
	<tr><th scope=row>TTTCCTCTCTGAAAGA-1_GSE116240_12W</th><td>GSE116240</td><td> 5975</td><td>2063</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.958159</td><td>0.22513089</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGCGCCAAGCGAGT-1_GSE116240_12W</th><td>GSE116240</td><td> 2860</td><td>1132</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.622378</td><td>0.35078534</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGCGCCACCAGATT-1_GSE116240_12W</th><td>GSE116240</td><td> 8023</td><td>2454</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 4.113175</td><td>0.24607330</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGCGCCATGTTCCC-1_GSE116240_12W</th><td>GSE116240</td><td>12195</td><td>2970</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.156622</td><td>0.28272251</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGCGCGTACACCGC-1_GSE116240_12W</th><td>GSE116240</td><td>10749</td><td>3016</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.511862</td><td>0.35078534</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGCGCTCCACGTTC-1_GSE116240_12W</th><td>GSE116240</td><td> 8859</td><td>2671</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 3.070324</td><td>0.24607330</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGCGCTCCGTCATC-1_GSE116240_12W</th><td>GSE116240</td><td>12353</td><td>2938</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.586659</td><td>0.21465969</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTAGACCTAGG-1_GSE116240_12W</th><td>GSE116240</td><td> 6789</td><td>2184</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.091619</td><td>0.39790576</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTAGACTGGGT-1_GSE116240_12W</th><td>GSE116240</td><td> 4506</td><td>1590</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.797603</td><td>0.12565445</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTCAAAGCAAT-1_GSE116240_12W</th><td>GSE116240</td><td> 3281</td><td>1351</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td>16.092655</td><td>0.13612565</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTGTGTATGGG-1_GSE116240_12W</th><td>GSE116240</td><td> 8742</td><td>2358</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.956074</td><td>0.21989529</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTGTTAGATGA-1_GSE116240_12W</th><td>GSE116240</td><td> 3315</td><td>1361</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 4.675716</td><td>0.11518325</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTGTTTAAGCC-1_GSE116240_12W</th><td>GSE116240</td><td>10808</td><td>2352</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.508142</td><td>0.10471204</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTTCGATCCCT-1_GSE116240_12W</th><td>GSE116240</td><td> 7667</td><td>2548</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.604278</td><td>0.11518325</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTTCGTGTAGT-1_GSE116240_12W</th><td>GSE116240</td><td> 8594</td><td>2794</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 6.388178</td><td>0.31937173</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGGTTTCTTGAGGT-1_GSE116240_12W</th><td>GSE116240</td><td> 9478</td><td>2283</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.719772</td><td>0.31937173</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTCAAGCGGCTTC-1_GSE116240_12W</th><td>GSE116240</td><td> 5576</td><td>1667</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.918938</td><td>0.12565445</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTCAAGTCCTCCT-1_GSE116240_12W</th><td>GSE116240</td><td> 9268</td><td>2731</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 1.143720</td><td>0.20418848</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTCACAGGCTCAC-1_GSE116240_12W</th><td>GSE116240</td><td> 7087</td><td>2012</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.356427</td><td>0.14136126</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTCACATCACCCT-1_GSE116240_12W</th><td>GSE116240</td><td>11528</td><td>2922</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.055864</td><td>0.27748691</td><td>Singlet</td></tr>
	<tr><th scope=row>TTTGTCAGTCGGCATC-1_GSE116240_12W</th><td>GSE116240</td><td>13308</td><td>3456</td><td>GSE116240_12W</td><td>12W</td><td>GSE116240</td><td> 2.983168</td><td>0.29842932</td><td>Singlet</td></tr>
</tbody>
</table>




```R
for (i in 1:3) {
   colnames(doublet.list[[i]]@meta.data)[9] <-"DoubletFinder"
}
for (i in 1:3) {
   colnames(doublet.list[[i]]@meta.data)[8] <-"DoubletScore"
}
```


```R
#merge
 doublefinders_seurat <- merge(doublet.list[[1]], y=c(doublet.list[[2]],doublet.list[[3]]))
```


```R
immune.combined[['DoubletFinder']] <- doublefinders_seurat@meta.data$DoubletFinder
immune.combined[['DoubletScore']] <- doublefinders_seurat@meta.data$DoubletScore
```


```R
#
saveRDS(immune.combined, file ='Step3_AddDoublet_Seurat.Rds')
```
