```R
###
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
```

    Attaching SeuratObject
    
    Warning message in system("timedatectl", intern = TRUE):
    “running command 'timedatectl' had status 1”
    ── [1mAttaching packages[22m ─────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
    
    [32m✔[39m [34mggplot2[39m 3.3.5     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.7     [32m✔[39m [34mdplyr  [39m 1.0.7
    [32m✔[39m [34mtidyr  [39m 1.1.4     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 2.1.1     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    
    Loading required package: grid
    
    ========================================
    ComplexHeatmap version 2.13.1
    Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    Github page: https://github.com/jokergoo/ComplexHeatmap
    Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    
    If you use it in published research, please cite:
    Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
      genomic data. Bioinformatics 2016.
    
    The new InteractiveComplexHeatmap package can directly export static 
    complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    
    This message can be suppressed by:
      suppressPackageStartupMessages(library(ComplexHeatmap))
    ========================================
    
    



```R
markers <-  data.table::fread('./Findmarkers_res_celltype_exclude_pro.csv')
```

    Warning message in data.table::fread("./Findmarkers_res_celltype_exclude_pro.csv"):
    “Detected 7 column names but the data has 8 columns (i.e. invalid file). Added 1 extra default column name for the first column which is guessed to be row names or an index. Use setnames() afterwards if this guess is not correct, or fix the file write command that created the file to create a valid file.”



```R
immune.combined <- readRDS('~/AS/AS_Mouse/AS_Mouse2/Final_result/Figure8_Human_Plaques/Step2_Seurat_Exclude_Proliferating.Rds')
```


```R
features <- c("S100A8","LYZ",'CD68',"ACTA2","MYL9",'RGS5',
             "PECAM1", "ACKR1","VWF",'IL7R','CD3D','TRBC2',
              "XCL1","GNLY","KLRD1","CD79A","CD19","VPREB3",
              "JCHAIN","IGHG4","LAMP5","CPA3","TPSD1","TPSAB1")
```


```R
cluster_info <- immune.combined@meta.data %>%
 arrange(celltype_Brief) %>%
 pull(celltype_Brief) %>%
 magrittr::set_names(rownames(immune.combined@meta.data))
cluster_info
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>AAACCCAAGATTAGAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCAAGCATGTTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCAAGCCTGTCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCAAGGGTTTCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCAAGTGCACCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCACAAGTTTGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCACATTCTCCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCAGTGTGTGTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCAGTGTTTCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCAGTTACCCTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCATCTCACTCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCATCTCATTAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACCCATCTCTCTAA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAAAGAAGTGTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAAAGAGTGAAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAACACAGCCTG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAACACTGGACC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAAGTAGGCTGA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAAGTCGTTGCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAAGTTATGTGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAAGTTCTCGTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAAGTTTACCAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAATCGATTCCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGAATCTCTTGCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTAGGAGCAAA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTAGGTTCCAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTAGTCCCGGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTAGTTACGTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTCACAAGCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTGTAGCTTGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTGTAGTTACC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTGTATAGCTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTGTATAGGAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTGTGGCGCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTGTTGATGTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTGTTGGAGGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTTCACAAGGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTTCAGACATC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTTCGACCAAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAACGCTTCGGTCACG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGAACAGAGGATCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGAACAGAGGCGGA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGAACAGCAACAGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGAACAGGCACTCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGAACAGTAGCTCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGAACCACTACCCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGAACGTTTGCCGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATAGAGGGTCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATAGGCGAACT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATAGTCACACT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATAGTCTGTAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATAGTGAGTGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATCAGGCTATT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATCAGTCAGAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATCATTCTGTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATGTTGGACTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATTCACGACTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATTCTGAACGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGATTCTTGGCTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCAGAATTGCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCAGCACACAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCGTATGTCAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCGTTGCCATA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCTCAAGCCTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCTCAGTCACA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCTCCGTTGGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGGCTCCTTGACC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTAAGGCACTAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTAAGTAAATGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTAAGTGCAGCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTACATTGCCTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTAGTCGAGCAA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTAGTGACTATC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTATCCGTCACT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTATCGAATCCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTATCTCTAAGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTATCTTACGTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGGTATCTTGGGCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCAGAAATCCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCAGGTAACTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCCACAAGCCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCGTCACAGTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCGTGCTCTCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCGTGGATGAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCGTTCCGCAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCTCCCGAACG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCTCGCTTGCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCTCGGTAGAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCTCTAGTCAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTCCTCTCAGGCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGAAGACGAGCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGAAGCTGCGAA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGAAGTCCCGAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGACAACAGCCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGACAAGTGGGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGAGTACAAGCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGAGTAGATTGA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAAGTGAGTGATGAAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAAGAAACTGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAAGATAACGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAAGCAATTAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAAGCCATATC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAAGCCATGCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAAGGCCTTGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAAGGGACACT-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGACAGCTGCCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGACAGTAGAGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGACATGCGGTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGACATTGTCGA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGAGTATATGGA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGATCCCATTTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGATCGTTAGTG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGATCTCCTGAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGATCTCGACGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AAATGGATCTTTCCGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGAGCGCCTAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGAGGGCTTCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGAGGTCATTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGCAATCACGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGCACAAAGTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGCACGCGCTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGCATTAGGAA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGGTACACGCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGGTAGAAACT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGGTTGTACGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGTCAATCTTC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGTCAGCACCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGTCCGTCCTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAAGTCCGTGACG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCAGCGCCGTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCAGGCATCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCAGTGCTCAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCAGTTCTCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCCAACAACAA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCGTCCAGTTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCGTCCCGGTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCGTTACGATC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCGTTCCGCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCTCAAACCTG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCTCATTACCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAACCTCTCCCATG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGAAGGCGCTCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGAAGGTAATCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGACAGAGAAAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGACAGAGGTTG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGACAGATCCAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGACATGTGGTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGAGTAGCGTTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGAGTCCTTTGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAAGAGTTAAGGAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACAGCTCGAAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACAGGTGATCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACAGTCATGGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACAGTCTGCGC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACCAACAGCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACCATCCTAAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACCATTCTCTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACTCACTGAAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACTCAGAATAG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACTCGCGATCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACTCTAACGCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACACACTCTAAGCGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGAGAAGTATC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGAGAATCGAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGCACGCTGTG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGCACGTACAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGCATAAGCGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGCATTCTTCA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGGTGCGTGCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGTCCAAGAGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGTCCATTGCC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGTCGATACGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACAGGGTCTACGGTA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACAGAAACTAC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACAGAGAGCCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACAGCCAGACA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACCAATGAGCG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACCACGCTTAA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACCAGCGAACA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACCATCCAACA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACGTCACCCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACTCATTTGCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACTCCACGTGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCAACTCCCTGTTG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACAAGACATAGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACAAGATCACCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACAAGCGGGTAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACACACGAGGAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACACAGCTCTGG-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACAGTCACTACA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACAGTCGTAATC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACAGTTTAGACC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACATCGAATGCT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCACATCGCAATGT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCATGAGATCCCAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCATGAGATGCCGA-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCATGAGCACTCAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCATGCAACAAGAT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCATGCAGCGTACC-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCATGCAGCTCCTT-1_AC1</dt><dd>Mono/Mø</dd><dt>AACCATGGTAAGTTGA-1_AC1</dt><dd>⋯</dd><dt>AACCATGGTCCGAAAG-1_AC1</dt><dd>MC</dd><dt>AACCATGGTGAATTAG-1_AC1</dt><dd>MC</dd><dt>AACCATGGTGATCGTT-1_AC1</dt><dd>MC</dd><dt>AACCATGGTGCATACT-1_AC1</dt><dd>MC</dd><dt>AACCATGGTGCTCTCT-1_AC1</dt><dd>MC</dd><dt>AACCATGTCTCATGCC-1_AC1</dt><dd>MC</dd><dt>AACCATGTCTTCTGTA-1_AC1</dt><dd>MC</dd><dt>AACCCAAAGACCATGG-1_AC1</dt><dd>MC</dd><dt>AACCCAAAGACCTTTG-1_AC1</dt><dd>MC</dd><dt>AACCCAAAGATCGCTT-1_AC1</dt><dd>MC</dd><dt>AACCCAACAAGTACCT-1_AC1</dt><dd>MC</dd><dt>AACCCAACACCTATCC-1_AC1</dt><dd>MC</dd><dt>AACCCAACACTGCTTC-1_AC1</dt><dd>MC</dd><dt>AACCCAACAGACAATA-1_AC1</dt><dd>MC</dd><dt>AACCCAACAGTCCCGA-1_AC1</dt><dd>MC</dd><dt>AACCCAACATATCTGG-1_AC1</dt><dd>MC</dd><dt>AACCCAACATGACAAA-1_AC1</dt><dd>MC</dd><dt>AACCCAACATTAAGCC-1_AC1</dt><dd>MC</dd><dt>AACCCAACATTACTCT-1_AC1</dt><dd>MC</dd><dt>AACCCAAGTCTCACGG-1_AC1</dt><dd>MC</dd><dt>AACCCAAGTGAGCAGT-1_AC1</dt><dd>MC</dd><dt>AACCCAATCCAGCAAT-1_AC1</dt><dd>MC</dd><dt>AACCCAATCTCCCATG-1_AC1</dt><dd>MC</dd><dt>AACCTGAAGACGCTCC-1_AC1</dt><dd>MC</dd><dt>AACCTGACAACTTCTT-1_AC1</dt><dd>MC</dd><dt>AACCTGACACAATGAA-1_AC1</dt><dd>MC</dd><dt>AACCTGACACAGTCAT-1_AC1</dt><dd>MC</dd><dt>AACCTGACAGCGAACA-1_AC1</dt><dd>MC</dd><dt>AACCTGACATAATGCC-1_AC1</dt><dd>MC</dd><dt>AACCTGAGTATCGTGT-1_AC1</dt><dd>MC</dd><dt>AACCTGATCCCAAGCG-1_AC1</dt><dd>MC</dd><dt>AACCTGATCCCATGGG-1_AC1</dt><dd>MC</dd><dt>AACCTGATCCGATGTA-1_AC1</dt><dd>MC</dd><dt>AACCTGATCCGTTTCG-1_AC1</dt><dd>MC</dd><dt>AACCTTTAGAGGGCGA-1_AC1</dt><dd>MC</dd><dt>AACCTTTAGGGTTAAT-1_AC1</dt><dd>MC</dd><dt>AACCTTTGTACTGGGA-1_AC1</dt><dd>MC</dd><dt>AACCTTTTCGGAAGGT-1_AC1</dt><dd>MC</dd><dt>AACCTTTTCTAGGCCG-1_AC1</dt><dd>MC</dd><dt>AACCTTTTCTCATGGA-1_AC1</dt><dd>MC</dd><dt>AACGAAACAACCGATT-1_AC1</dt><dd>MC</dd><dt>AACGAAACACCAATTG-1_AC1</dt><dd>MC</dd><dt>AACGAAACATGCCGAC-1_AC1</dt><dd>MC</dd><dt>AACGAAAGTCCTTTGC-1_AC1</dt><dd>MC</dd><dt>AACGAAAGTCGCCTAG-1_AC1</dt><dd>MC</dd><dt>AACGAAAGTGAGCAGT-1_AC1</dt><dd>MC</dd><dt>AACGAAAGTTAGAAAC-1_AC1</dt><dd>MC</dd><dt>AACGAAATCAGCAGAG-1_AC1</dt><dd>MC</dd><dt>AACGGGAAGAGAGTGA-1_AC1</dt><dd>MC</dd><dt>AACGGGAAGATGATTG-1_AC1</dt><dd>MC</dd><dt>AACGGGACAGACCAAG-1_AC1</dt><dd>MC</dd><dt>AACGGGACAGGACAGT-1_AC1</dt><dd>MC</dd><dt>AACGGGACATCCTCAC-1_AC1</dt><dd>MC</dd><dt>AACGGGAGTGGGTCAA-1_AC1</dt><dd>MC</dd><dt>AACGGGATCGAGTCTA-1_AC1</dt><dd>MC</dd><dt>AACGTCAAGACATGCG-1_AC1</dt><dd>MC</dd><dt>AACGTCAAGCGTGTTT-1_AC1</dt><dd>MC</dd><dt>AACGTCAAGTCTCCTC-1_AC1</dt><dd>MC</dd><dt>AACGTCAAGTGGTTAA-1_AC1</dt><dd>MC</dd><dt>AACGTCACAAACCGGA-1_AC1</dt><dd>MC</dd><dt>AACGTCAGTCACGACC-1_AC1</dt><dd>MC</dd><dt>AACGTCAGTGACAGGT-1_AC1</dt><dd>MC</dd><dt>AACGTCAGTGGACCAA-1_AC1</dt><dd>MC</dd><dt>AACGTCAGTTACGTAC-1_AC1</dt><dd>MC</dd><dt>AACGTCATCAATCTTC-1_AC1</dt><dd>MC</dd><dt>AACGTCATCCATCTGC-1_AC1</dt><dd>MC</dd><dt>AACGTCATCGTTCGCT-1_AC1</dt><dd>MC</dd><dt>AACTTCTAGAAGGTAG-1_AC1</dt><dd>MC</dd><dt>AACTTCTGTACGACTT-1_AC1</dt><dd>MC</dd><dt>AACTTCTGTCCCACGA-1_AC1</dt><dd>MC</dd><dt>AACTTCTGTCGATTTG-1_AC1</dt><dd>MC</dd><dt>AACTTCTGTCTACGAT-1_AC1</dt><dd>MC</dd><dt>AACTTCTTCCAAGAGG-1_AC1</dt><dd>MC</dd><dt>AACTTCTTCCAGCACG-1_AC1</dt><dd>MC</dd><dt>AACTTCTTCGGATAAA-1_AC1</dt><dd>MC</dd><dt>AAGAACAAGGATGTTA-1_AC1</dt><dd>MC</dd><dt>AAGAACACACTGGAAG-1_AC1</dt><dd>MC</dd><dt>AAGAACACAGTGTATC-1_AC1</dt><dd>MC</dd><dt>AAGAACAGTCCAGGTC-1_AC1</dt><dd>MC</dd><dt>AAGAACAGTGTAGCAG-1_AC1</dt><dd>MC</dd><dt>AAGAACATCGGTCGAC-1_AC1</dt><dd>MC</dd><dt>AAGACAAAGCCACAAG-1_AC1</dt><dd>MC</dd><dt>AAGACAACAATCACGT-1_AC1</dt><dd>MC</dd><dt>AAGACAACATCAGCTA-1_AC1</dt><dd>MC</dd><dt>AAGACAAGTGTGAGCA-1_AC1</dt><dd>MC</dd><dt>AAGACAAGTTCGGTAT-1_AC1</dt><dd>MC</dd><dt>AAGACAATCTCAGAAC-1_AC1</dt><dd>MC</dd><dt>AAGACTCAGAAACACT-1_AC1</dt><dd>MC</dd><dt>AAGACTCAGAGGGTAA-1_AC1</dt><dd>MC</dd><dt>AAGACTCAGCGAAACC-1_AC1</dt><dd>MC</dd><dt>AAGACTCAGCGTCTCG-1_AC1</dt><dd>MC</dd><dt>AAGACTCCACAAAGCG-1_AC1</dt><dd>MC</dd><dt>AAGACTCCACCCTAGG-1_AC1</dt><dd>MC</dd><dt>AAGACTCGTATACAGA-1_AC1</dt><dd>MC</dd><dt>AAGACTCGTCACCGCA-1_AC1</dt><dd>MC</dd><dt>AAGACTCTCACTTGTT-1_AC1</dt><dd>MC</dd><dt>AAGACTCTCGTAGGGA-1_AC1</dt><dd>MC</dd><dt>AAGATAGAGCTGGCTC-1_AC1</dt><dd>MC</dd><dt>AAGATAGAGGGCTAAC-1_AC1</dt><dd>MC</dd><dt>AAGATAGCAGGCTCTG-1_AC1</dt><dd>MC</dd><dt>AAGATAGGTATTTCTC-1_AC1</dt><dd>MC</dd><dt>AAGATAGGTGTTAGCT-1_AC1</dt><dd>MC</dd><dt>AAGATAGTCGCCTTTG-1_AC1</dt><dd>MC</dd><dt>AAGATAGTCGGAATTC-1_AC1</dt><dd>MC</dd><dt>AAGATAGTCTTCGATT-1_AC1</dt><dd>MC</dd><dt>AAGCATCAGGACAGTC-1_AC1</dt><dd>MC</dd><dt>AAGCATCCAAGGCCTC-1_AC1</dt><dd>MC</dd><dt>AAGCATCCACGTCTCT-1_AC1</dt><dd>MC</dd><dt>AAGCATCGTAACTAAG-1_AC1</dt><dd>MC</dd><dt>AAGCCATAGCAGCCCT-1_AC1</dt><dd>MC</dd><dt>AAGCCATAGTATGAGT-1_AC1</dt><dd>MC</dd><dt>AAGCCATCACCAGACC-1_AC1</dt><dd>MC</dd><dt>AAGCCATCAGGGAATC-1_AC1</dt><dd>MC</dd><dt>AAGCCATCAGTTGAAA-1_AC1</dt><dd>MC</dd><dt>AAGCCATTCAAGGTGG-1_AC1</dt><dd>MC</dd><dt>AAGCCATTCCGTACGG-1_AC1</dt><dd>MC</dd><dt>AAGCCATTCGTAACAC-1_AC1</dt><dd>MC</dd><dt>AAGCGAGCACACCTGG-1_AC1</dt><dd>MC</dd><dt>AAGCGAGCACGCAAAG-1_AC1</dt><dd>MC</dd><dt>AAGCGAGCAGCGCGTT-1_AC1</dt><dd>MC</dd><dt>AAGCGAGGTAAGAACT-1_AC1</dt><dd>MC</dd><dt>AAGCGAGGTCGGAACA-1_AC1</dt><dd>MC</dd><dt>AAGCGTTAGCCGCACT-1_AC1</dt><dd>MC</dd><dt>AAGCGTTAGCGATCGA-1_AC1</dt><dd>MC</dd><dt>AAGCGTTAGTAGACAT-1_AC1</dt><dd>MC</dd><dt>AAGCGTTAGTCTAGCT-1_AC1</dt><dd>MC</dd><dt>AAGCGTTCAACACTAC-1_AC1</dt><dd>MC</dd><dt>AAGCGTTCACCATAAC-1_AC1</dt><dd>MC</dd><dt>AAGCGTTGTCGAGATG-1_AC1</dt><dd>MC</dd><dt>AAGCGTTGTTCTGACA-1_AC1</dt><dd>MC</dd><dt>AAGCGTTGTTTGGCTA-1_AC1</dt><dd>MC</dd><dt>AAGCGTTTCATGTCAG-1_AC1</dt><dd>MC</dd><dt>AAGCGTTTCCCTTTGG-1_AC1</dt><dd>MC</dd><dt>AAGCGTTTCCGTGGGT-1_AC1</dt><dd>MC</dd><dt>AAGCGTTTCTACTATC-1_AC1</dt><dd>MC</dd><dt>AAGGAATAGCCTAGGA-1_AC1</dt><dd>MC</dd><dt>AAGGAATAGCTCCGAC-1_AC1</dt><dd>MC</dd><dt>AAGGAATAGTCATGCT-1_AC1</dt><dd>MC</dd><dt>AAGGAATCAAGCGATG-1_AC1</dt><dd>MC</dd><dt>AAGGAATCACCCAATA-1_AC1</dt><dd>MC</dd><dt>AAGGAATCATCCGAGC-1_AC1</dt><dd>MC</dd><dt>AAGGAATCATCGGTTA-1_AC1</dt><dd>MC</dd><dt>AAGGAATCATTAAGCC-1_AC1</dt><dd>MC</dd><dt>AAGGAATGTCACCACG-1_AC1</dt><dd>MC</dd><dt>AAGGAATTCATAGCAC-1_AC1</dt><dd>MC</dd><dt>AAGGAATTCCTCTAGC-1_AC1</dt><dd>MC</dd><dt>AAGGAATTCGAAATCC-1_AC1</dt><dd>MC</dd><dt>AAGGAATTCTTCTTCC-1_AC1</dt><dd>MC</dd><dt>AAGGTAAAGCATCAGG-1_AC1</dt><dd>MC</dd><dt>AAGGTAAAGTCTGTAC-1_AC1</dt><dd>MC</dd><dt>AAGGTAAAGTGGAAAG-1_AC1</dt><dd>MC</dd><dt>AAGGTAACATAGCTGT-1_AC1</dt><dd>MC</dd><dt>AAGGTAAGTCTTCCGT-1_AC1</dt><dd>MC</dd><dt>AAGGTAAGTGAGCTCC-1_AC1</dt><dd>MC</dd><dt>AAGGTAATCAAAGGTA-1_AC1</dt><dd>MC</dd><dt>AAGTACCAGCTAGTTC-1_AC1</dt><dd>MC</dd><dt>AAGTACCCAACAGTGG-1_AC1</dt><dd>MC</dd><dt>AAGTACCCATATCTCT-1_AC1</dt><dd>MC</dd><dt>AAGTACCCATGGAAGC-1_AC1</dt><dd>MC</dd><dt>AAGTACCGTAGACACG-1_AC1</dt><dd>MC</dd><dt>AAGTACCTCAGTGTCA-1_AC1</dt><dd>MC</dd><dt>AAGTACCTCTCCTACG-1_AC1</dt><dd>MC</dd><dt>AAGTACCTCTTTGCTA-1_AC1</dt><dd>MC</dd><dt>AAGTCGTAGCTGAGTG-1_AC1</dt><dd>MC</dd><dt>AAGTCGTAGTCAGCCC-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCAAGCAGGT-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCAAGCTGTT-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCAATGTTGC-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCACACCGAC-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCAGACTCTA-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCAGAGTCAG-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCAGGTAGTG-1_AC1</dt><dd>MC</dd><dt>AAGTCGTCAGTGTGGA-1_AC1</dt><dd>MC</dd><dt>AAGTCGTGTATCTCGA-1_AC1</dt><dd>MC</dd><dt>AAGTCGTGTCAGGTAG-1_AC1</dt><dd>MC</dd><dt>AAGTCGTGTCATGACT-1_AC1</dt><dd>MC</dd><dt>AAGTCGTGTGGAACAC-1_AC1</dt><dd>MC</dd><dt>AAGTCGTTCAAATAGG-1_AC1</dt><dd>MC</dd><dt>AAGTCGTTCGTTCTCG-1_AC1</dt><dd>MC</dd><dt>AAGTGAAAGCTCGTGC-1_AC1</dt><dd>MC</dd><dt>AAGTGAACACATTGTG-1_AC1</dt><dd>MC</dd><dt>AAGTGAACATGCCGGT-1_AC1</dt><dd>MC</dd><dt>AAGTGAAGTCCACAGC-1_AC1</dt><dd>MC</dd><dt>AAGTGAAGTTTCCAAG-1_AC1</dt><dd>MC</dd><dt>AAGTGAATCAAGAGTA-1_AC1</dt><dd>MC</dd><dt>AAGTGAATCACTACTT-1_AC1</dt><dd>MC</dd><dt>AAGTGAATCCGTGGTG-1_AC1</dt><dd>MC</dd><dt>AAGTGAATCGTCTAAG-1_AC1</dt><dd>MC</dd><dt>AAGTGAATCTGGCCTT-1_AC1</dt><dd>MC</dd><dt>AAGTTCGAGCGAAACC-1_AC1</dt><dd>MC</dd><dt>AAGTTCGAGTAGACAT-1_AC1</dt><dd>MC</dd><dt>AAGTTCGAGTCAATCC-1_AC1</dt><dd>MC</dd><dt>AAGTTCGAGTCATGGG-1_AC1</dt><dd>MC</dd><dt>AAGTTCGAGTCTCCTC-1_AC1</dt><dd>MC</dd><dt>AAGTTCGTCCAACCGG-1_AC1</dt><dd>MC</dd><dt>AAGTTCGTCCGCTGTT-1_AC1</dt><dd>MC</dd><dt>AAGTTCGTCTGGTTGA-1_AC1</dt><dd>MC</dd><dt>AATAGAGAGAGTGACC-1_AC1</dt><dd>MC</dd><dt>AATAGAGAGCGAGTAC-1_AC1</dt><dd>MC</dd><dt>AATAGAGAGGCATGCA-1_AC1</dt><dd>MC</dd></dl>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Mono/Mø'</li><li>'SMC'</li><li>'EC'</li><li>'T'</li><li>'NK'</li><li>'B'</li><li>'PC'</li><li>'MC'</li></ol>
</details>



```R
mat <- immune.combined@assays$RNA@data
mat <- as.matrix(mat[features, names(cluster_info)])
immune.combined@meta.data %>% select(celltype_Brief)
```


<table class="dataframe">
<caption>A data.frame: 48235 × 1</caption>
<thead>
	<tr><th></th><th scope=col>celltype_Brief</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACCCAAGATTAGAC-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCAAGCATGTTC-1_AC1</th><td>T      </td></tr>
	<tr><th scope=row>AAACCCAAGCCTGTCG-1_AC1</th><td>T      </td></tr>
	<tr><th scope=row>AAACCCAAGGGTTTCT-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCAAGTGCACCC-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCACAAGTTTGC-1_AC1</th><td>B      </td></tr>
	<tr><th scope=row>AAACCCACATTCTCCG-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCAGTGTGTGTT-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCAGTGTTTCTT-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCAGTTACCCTC-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCATCTCACTCG-1_AC1</th><td>T      </td></tr>
	<tr><th scope=row>AAACCCATCTCATTAC-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACCCATCTCTCTAA-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGAAAGAAGTGTT-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGAAAGAGTGAAG-1_AC1</th><td>SMC    </td></tr>
	<tr><th scope=row>AAACGAACACAGCCTG-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGAACACTGGACC-1_AC1</th><td>EC     </td></tr>
	<tr><th scope=row>AAACGAAGTAGGCTGA-1_AC1</th><td>T      </td></tr>
	<tr><th scope=row>AAACGAAGTCGTTGCG-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGAAGTTATGTGC-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGAAGTTCTCGTC-1_AC1</th><td>SMC    </td></tr>
	<tr><th scope=row>AAACGAAGTTTACCAG-1_AC1</th><td>SMC    </td></tr>
	<tr><th scope=row>AAACGAATCGATTCCC-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGAATCTCTTGCG-1_AC1</th><td>NK     </td></tr>
	<tr><th scope=row>AAACGCTAGGAGCAAA-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGCTAGGTTCCAT-1_AC1</th><td>T      </td></tr>
	<tr><th scope=row>AAACGCTAGTCCCGGT-1_AC1</th><td>T      </td></tr>
	<tr><th scope=row>AAACGCTAGTTACGTC-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGCTCACAAGCTT-1_AC1</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>AAACGCTGTAGCTTGT-1_AC1</th><td>SMC    </td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td></tr>
	<tr><th scope=row>TTTACTGGTCCTGTCT-1_PA3</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>TTTACTGGTCTTCTAT-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTACTGGTGATGTAA-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTACTGGTTTGCCGG-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTACTGTCCTGGGAC-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTATGCAGTAGCATA-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTATGCGTGCATTTG-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTCACAAGTCGGCAA-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTCACAAGTTTGAGA-1_PA3</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>TTTCAGTAGGATTTCC-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTCAGTCAAAGGCAC-1_PA3</th><td>T      </td></tr>
	<tr><th scope=row>TTTCATGCAACCAATC-1_PA3</th><td>T      </td></tr>
	<tr><th scope=row>TTTCATGCACAATGAA-1_PA3</th><td>T      </td></tr>
	<tr><th scope=row>TTTCATGGTATCGAGG-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTCCTCAGCACTCTA-1_PA3</th><td>T      </td></tr>
	<tr><th scope=row>TTTCCTCCATTGAGCT-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTCCTCTCAACTTTC-1_PA3</th><td>T      </td></tr>
	<tr><th scope=row>TTTCCTCTCAATCCGA-1_PA3</th><td>B      </td></tr>
	<tr><th scope=row>TTTCGATTCGATACTG-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTCGATTCTGACAGT-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTGACTAGTGATAGT-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTGACTCAATACGAA-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTGACTCATGTGACT-1_PA3</th><td>T      </td></tr>
	<tr><th scope=row>TTTGACTGTTCGAAGG-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTGACTGTTCGGCTG-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTGATCTCGTCAACA-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTGGAGGTCGCATGC-1_PA3</th><td>EC     </td></tr>
	<tr><th scope=row>TTTGGAGTCATGGATC-1_PA3</th><td>SMC    </td></tr>
	<tr><th scope=row>TTTGGTTCAAGAAATC-1_PA3</th><td>Mono/Mø</td></tr>
	<tr><th scope=row>TTTGTTGTCTTAGCAG-1_PA3</th><td>T      </td></tr>
</tbody>
</table>




```R
immune.combined@meta.data %>% pull(sample)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>'patient_1AC'</li><li>⋯</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li><li>'patient_3PA'</li></ol>




```R
zscore_matrix1 <- t(mat) %>%
as.data.frame() %>%
rownames_to_column(var ='cell') %>%
bind_cols(celltype =immune.combined@meta.data %>% pull(celltype_Brief), sample =immune.combined@meta.data %>% pull(sample)) %>%
pivot_longer(-c("cell","celltype",'sample'),names_to = "gene") %>% 
group_by(celltype, gene, sample) %>% 
summarize(mean_celltype_expression = mean(value)) %>%
ungroup()
```

    `summarise()` has grouped output by 'celltype', 'gene'. You can override using the `.groups` argument.



```R
mean_gene_matrix <- zscore_matrix1 %>%
group_by(gene) %>% 
summarise(mean_expression = mean(mean_celltype_expression),sd_expression = sd(mean_celltype_expression))
mean_gene_matrix
```


<table class="dataframe">
<caption>A tibble: 24 × 3</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>mean_expression</th><th scope=col>sd_expression</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>ACKR1 </td><td>0.23898220</td><td>0.7314954</td></tr>
	<tr><td>ACTA2 </td><td>0.38868561</td><td>0.8159027</td></tr>
	<tr><td>CD19  </td><td>0.05389473</td><td>0.1561799</td></tr>
	<tr><td>CD3D  </td><td>0.26231918</td><td>0.5472119</td></tr>
	<tr><td>CD68  </td><td>0.37543524</td><td>0.5767371</td></tr>
	<tr><td>CD79A </td><td>0.39118157</td><td>0.8543917</td></tr>
	<tr><td>CPA3  </td><td>0.39317289</td><td>1.0556342</td></tr>
	<tr><td>GNLY  </td><td>0.41040524</td><td>0.8177187</td></tr>
	<tr><td>IGHG4 </td><td>0.09052218</td><td>0.2668570</td></tr>
	<tr><td>IL7R  </td><td>0.32146007</td><td>0.6122635</td></tr>
	<tr><td>JCHAIN</td><td>0.48166283</td><td>1.0425785</td></tr>
	<tr><td>KLRD1 </td><td>0.31905467</td><td>0.6920305</td></tr>
	<tr><td>LAMP5 </td><td>0.10539433</td><td>0.3302265</td></tr>
	<tr><td>LYZ   </td><td>0.40288202</td><td>0.8807704</td></tr>
	<tr><td>MYL9  </td><td>0.45564855</td><td>0.8973546</td></tr>
	<tr><td>PECAM1</td><td>0.34954397</td><td>0.6037583</td></tr>
	<tr><td>RGS5  </td><td>0.19971368</td><td>0.4737646</td></tr>
	<tr><td>S100A8</td><td>0.25882727</td><td>0.6128939</td></tr>
	<tr><td>TPSAB1</td><td>0.49943878</td><td>1.3650055</td></tr>
	<tr><td>TPSD1 </td><td>0.04832435</td><td>0.1937017</td></tr>
	<tr><td>TRBC2 </td><td>0.47615446</td><td>0.6366003</td></tr>
	<tr><td>VPREB3</td><td>0.09914959</td><td>0.2576911</td></tr>
	<tr><td>VWF   </td><td>0.24060745</td><td>0.6179149</td></tr>
	<tr><td>XCL1  </td><td>0.29367187</td><td>0.6917485</td></tr>
</tbody>
</table>




```R
zscore_matrix_final <- zscore_matrix1 %>%
group_by(gene) %>% 
mutate(mean_expression = mean(mean_celltype_expression),sd_expression = sd(mean_celltype_expression)) %>%
group_by(celltype, gene) %>%
summarise(mean_celltype_expression = mean(mean_celltype_expression)) %>%
ungroup() %>%
left_join(mean_gene_matrix,by ='gene') %>%
mutate(zscore = (mean_celltype_expression-mean_expression)/sd_expression)
```

    `summarise()` has grouped output by 'celltype'. You can override using the `.groups` argument.



```R
zscore_matrix_final$gene <- factor(zscore_matrix_final$gene,levels=  c("S100A8","LYZ",'CD68',"ACTA2","MYL9",'RGS5',
             "PECAM1", "ACKR1","VWF",'IL7R','CD3D','TRBC2',
              "XCL1","GNLY","KLRD1","CD79A","CD19","VPREB3",
              "JCHAIN","IGHG4","LAMP5","CPA3","TPSD1","TPSAB1"))
```


```R
zscore_matrix_final$celltype
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>T</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>NK</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>B</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>PC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li><li>MC</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Mono/Mø'</li><li>'SMC'</li><li>'EC'</li><li>'T'</li><li>'NK'</li><li>'B'</li><li>'PC'</li><li>'MC'</li></ol>
</details>



```R
zscore_matrix_final$celltype <-factor(zscore_matrix_final$celltype, levels = c('Mono/Mø','SMC','EC','T',
                                                                              'NK','B','PC','MC'))
```


```R
mat_heatmap <- zscore_matrix_final %>%
select(celltype, gene, zscore) %>%
pivot_wider(names_from = 'celltype', values_from = 'zscore') %>%
as.data.frame() %>%
column_to_rownames(var = 'gene')
```


```R
immune.combined@meta.data$celltype_Brief
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>Mono/Mø</li><li>T</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>B</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>NK</li><li>Mono/Mø</li><li>T</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>T</li><li>SMC</li><li>T</li><li>EC</li><li>T</li><li>T</li><li>T</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>EC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>NK</li><li>NK</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>T</li><li>NK</li><li>Mono/Mø</li><li>T</li><li>SMC</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>EC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>T</li><li>NK</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>NK</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>EC</li><li>T</li><li>T</li><li>T</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>T</li><li>EC</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>T</li><li>T</li><li>SMC</li><li>T</li><li>EC</li><li>Mono/Mø</li><li>T</li><li>SMC</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>T</li><li>T</li><li>Mono/Mø</li><li>SMC</li><li>T</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>MC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>SMC</li><li>T</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>MC</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>NK</li><li>Mono/Mø</li><li>Mono/Mø</li><li>EC</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>EC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>T</li><li>SMC</li><li>EC</li><li>T</li><li>Mono/Mø</li><li>NK</li><li>T</li><li>⋯</li><li>EC</li><li>EC</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>MC</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>SMC</li><li>NK</li><li>SMC</li><li>EC</li><li>Mono/Mø</li><li>EC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>EC</li><li>T</li><li>EC</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>T</li><li>SMC</li><li>SMC</li><li>SMC</li><li>EC</li><li>EC</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>EC</li><li>SMC</li><li>EC</li><li>Mono/Mø</li><li>SMC</li><li>T</li><li>EC</li><li>Mono/Mø</li><li>SMC</li><li>EC</li><li>T</li><li>Mono/Mø</li><li>SMC</li><li>EC</li><li>T</li><li>T</li><li>EC</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>SMC</li><li>SMC</li><li>EC</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>MC</li><li>SMC</li><li>T</li><li>Mono/Mø</li><li>NK</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>NK</li><li>SMC</li><li>SMC</li><li>SMC</li><li>EC</li><li>EC</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>SMC</li><li>EC</li><li>SMC</li><li>EC</li><li>MC</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>EC</li><li>T</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>EC</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>EC</li><li>SMC</li><li>T</li><li>Mono/Mø</li><li>SMC</li><li>T</li><li>MC</li><li>EC</li><li>T</li><li>Mono/Mø</li><li>T</li><li>EC</li><li>SMC</li><li>EC</li><li>T</li><li>B</li><li>SMC</li><li>SMC</li><li>EC</li><li>SMC</li><li>SMC</li><li>EC</li><li>T</li><li>EC</li><li>SMC</li><li>SMC</li><li>EC</li><li>SMC</li><li>EC</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>SMC</li><li>SMC</li><li>EC</li><li>MC</li><li>SMC</li><li>EC</li><li>SMC</li><li>SMC</li><li>T</li><li>Mono/Mø</li><li>SMC</li><li>Mono/Mø</li><li>Mono/Mø</li><li>SMC</li><li>EC</li><li>SMC</li><li>Mono/Mø</li><li>T</li><li>SMC</li><li>T</li><li>T</li><li>NK</li><li>T</li><li>Mono/Mø</li><li>NK</li><li>EC</li><li>SMC</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>SMC</li><li>SMC</li><li>EC</li><li>SMC</li><li>SMC</li><li>Mono/Mø</li><li>EC</li><li>EC</li><li>EC</li><li>EC</li><li>SMC</li><li>EC</li><li>SMC</li><li>Mono/Mø</li><li>SMC</li><li>T</li><li>T</li><li>T</li><li>EC</li><li>T</li><li>SMC</li><li>T</li><li>B</li><li>SMC</li><li>EC</li><li>SMC</li><li>EC</li><li>T</li><li>SMC</li><li>SMC</li><li>EC</li><li>EC</li><li>SMC</li><li>Mono/Mø</li><li>T</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Mono/Mø'</li><li>'SMC'</li><li>'EC'</li><li>'T'</li><li>'NK'</li><li>'B'</li><li>'PC'</li><li>'MC'</li></ol>
</details>



```R
cluster_info <-  c('Mono/Mø','SMC','EC','T', 'NK','B','PC','MC')
names(cluster_info) <- c('Mono/Mø','SMC','EC','T', 'NK','B','PC','MC')
```


```R
col <- c('#F6313E','#0081C9','#FFA300','#46A040','#89774A','#A65AC2','#0DB2AA','#00441B')
names(col) <- levels(immune.combined@meta.data$celltype_Brief)
```


```R
mat_heatmap <- mat_heatmap[features,]
```


```R
mat_heatmap
```


<table class="dataframe">
<caption>A data.frame: 24 × 8</caption>
<thead>
	<tr><th></th><th scope=col>Mono/Mø</th><th scope=col>SMC</th><th scope=col>EC</th><th scope=col>T</th><th scope=col>NK</th><th scope=col>B</th><th scope=col>PC</th><th scope=col>MC</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>S100A8</th><td> 2.19172032</td><td>-0.3587418</td><td>-0.346451254</td><td>-0.35178122</td><td>-0.34478465</td><td>-0.21574790</td><td>-0.19108576</td><td>-0.3831277</td></tr>
	<tr><th scope=row>LYZ</th><td> 2.61061099</td><td>-0.3900233</td><td>-0.379780474</td><td>-0.38115520</td><td>-0.38558709</td><td>-0.39956447</td><td>-0.37518891</td><td>-0.2993115</td></tr>
	<tr><th scope=row>CD68</th><td> 2.29405000</td><td>-0.3755227</td><td>-0.478591686</td><td>-0.60086329</td><td>-0.60132238</td><td>-0.54220082</td><td> 0.60521411</td><td>-0.3007632</td></tr>
	<tr><th scope=row>ACTA2</th><td>-0.32522587</td><td> 2.4978430</td><td>-0.361179484</td><td>-0.31447034</td><td>-0.33650965</td><td>-0.39475480</td><td>-0.40041725</td><td>-0.3652856</td></tr>
	<tr><th scope=row>MYL9</th><td>-0.39144236</td><td> 2.5466548</td><td> 0.005691627</td><td>-0.42715197</td><td>-0.42603367</td><td>-0.43050904</td><td>-0.45228473</td><td>-0.4249246</td></tr>
	<tr><th scope=row>RGS5</th><td>-0.34160997</td><td> 2.4085310</td><td>-0.313432207</td><td>-0.36782629</td><td>-0.35430325</td><td>-0.37667493</td><td>-0.36325635</td><td>-0.2914280</td></tr>
	<tr><th scope=row>PECAM1</th><td>-0.01321306</td><td>-0.5420676</td><td> 2.451116864</td><td>-0.52268235</td><td>-0.44358258</td><td>-0.50491252</td><td>-0.04029025</td><td>-0.3843685</td></tr>
	<tr><th scope=row>ACKR1</th><td>-0.26731847</td><td>-0.2936880</td><td> 2.017285454</td><td>-0.30260631</td><td>-0.30657200</td><td>-0.31393696</td><td>-0.31724439</td><td>-0.2159193</td></tr>
	<tr><th scope=row>VWF</th><td>-0.34492318</td><td>-0.3647904</td><td> 2.512063579</td><td>-0.36243864</td><td>-0.37362148</td><td>-0.37516293</td><td>-0.38105676</td><td>-0.3100702</td></tr>
	<tr><th scope=row>IL7R</th><td>-0.36643421</td><td>-0.4206861</td><td>-0.402364896</td><td> 2.44048773</td><td>-0.02392908</td><td>-0.35328054</td><td>-0.39981226</td><td>-0.4739806</td></tr>
	<tr><th scope=row>CD3D</th><td>-0.39289934</td><td>-0.3861210</td><td>-0.389543632</td><td> 2.59493925</td><td>-0.25575680</td><td>-0.35726755</td><td>-0.43175978</td><td>-0.3815911</td></tr>
	<tr><th scope=row>TRBC2</th><td>-0.66251729</td><td>-0.6491314</td><td>-0.659935973</td><td> 2.08549570</td><td> 0.91244893</td><td> 0.30059331</td><td>-0.68574127</td><td>-0.6412120</td></tr>
	<tr><th scope=row>XCL1</th><td>-0.40208660</td><td>-0.4031994</td><td>-0.409181196</td><td>-0.06628764</td><td> 2.50737052</td><td>-0.40549570</td><td>-0.41620561</td><td>-0.4049144</td></tr>
	<tr><th scope=row>GNLY</th><td>-0.45167809</td><td>-0.4651669</td><td>-0.451330310</td><td> 0.48971733</td><td> 2.31092648</td><td>-0.48455206</td><td>-0.47596763</td><td>-0.4719488</td></tr>
	<tr><th scope=row>KLRD1</th><td>-0.44848203</td><td>-0.4438365</td><td>-0.435403154</td><td> 0.16811804</td><td> 2.51677128</td><td>-0.44889421</td><td>-0.45287542</td><td>-0.4553980</td></tr>
	<tr><th scope=row>CD79A</th><td>-0.44914597</td><td>-0.4475425</td><td>-0.445503746</td><td>-0.44668185</td><td>-0.45003735</td><td> 2.51913644</td><td> 0.17096644</td><td>-0.4511914</td></tr>
	<tr><th scope=row>CD19</th><td>-0.31540437</td><td>-0.3222459</td><td>-0.329529870</td><td>-0.34320801</td><td>-0.33552019</td><td> 2.22737967</td><td>-0.23639034</td><td>-0.3450810</td></tr>
	<tr><th scope=row>VPREB3</th><td>-0.38156108</td><td>-0.3695394</td><td>-0.374206297</td><td>-0.38212388</td><td>-0.38088328</td><td> 2.54380214</td><td>-0.27072676</td><td>-0.3847615</td></tr>
	<tr><th scope=row>JCHAIN</th><td>-0.42041578</td><td>-0.4274181</td><td>-0.426251411</td><td>-0.41816066</td><td>-0.42340786</td><td>-0.06339236</td><td> 2.55616477</td><td>-0.3771186</td></tr>
	<tr><th scope=row>IGHG4</th><td>-0.31115968</td><td>-0.3159909</td><td>-0.320917645</td><td>-0.31830295</td><td>-0.28839881</td><td> 0.08615579</td><td> 1.79790862</td><td>-0.3292944</td></tr>
	<tr><th scope=row>LAMP5</th><td>-0.28308899</td><td>-0.3059638</td><td>-0.318418925</td><td>-0.31876273</td><td>-0.31915770</td><td>-0.28229340</td><td> 2.14684323</td><td>-0.3191577</td></tr>
	<tr><th scope=row>CPA3</th><td>-0.37065702</td><td>-0.3689011</td><td>-0.367324190</td><td>-0.37068935</td><td>-0.35695510</td><td>-0.36862066</td><td>-0.37245182</td><td> 2.5755992</td></tr>
	<tr><th scope=row>TPSD1</th><td>-0.24834909</td><td>-0.2482018</td><td>-0.246503447</td><td>-0.24947821</td><td>-0.23587129</td><td>-0.24947821</td><td>-0.23378656</td><td> 1.7116686</td></tr>
	<tr><th scope=row>TPSAB1</th><td>-0.36255106</td><td>-0.3615543</td><td>-0.359244833</td><td>-0.36276750</td><td>-0.35717963</td><td>-0.36572738</td><td>-0.36017848</td><td> 2.5292031</td></tr>
</tbody>
</table>




```R
cluster_info <- factor(cluster_info, levels = c('Mono/Mø','SMC','EC','T','NK','B','PC','MC'))
```


```R
cell_ratio <- immune.combined@meta.data %>%
group_by(celltype_Brief) %>%
summarise(count = n()) %>%
mutate(total = sum(count)) %>%
mutate(ratio = count/total)
```


```R
cell_ratio
```


<table class="dataframe">
<caption>A tibble: 8 × 4</caption>
<thead>
	<tr><th scope=col>celltype_Brief</th><th scope=col>count</th><th scope=col>total</th><th scope=col>ratio</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Mono/Mø</td><td>10001</td><td>48235</td><td>0.20733907</td></tr>
	<tr><td>SMC    </td><td>10289</td><td>48235</td><td>0.21330984</td></tr>
	<tr><td>EC     </td><td> 6312</td><td>48235</td><td>0.13085933</td></tr>
	<tr><td>T      </td><td>17286</td><td>48235</td><td>0.35837048</td></tr>
	<tr><td>NK     </td><td> 1606</td><td>48235</td><td>0.03329532</td></tr>
	<tr><td>B      </td><td> 1884</td><td>48235</td><td>0.03905877</td></tr>
	<tr><td>PC     </td><td>  311</td><td>48235</td><td>0.00644760</td></tr>
	<tr><td>MC     </td><td>  546</td><td>48235</td><td>0.01131958</td></tr>
</tbody>
</table>




```R
top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.5, col = "white")),
  cell_ratio = anno_barplot(cell_ratio$ratio,)) # 设置字体
```


```R
cluster_info
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>Mono/Mø</dt><dd>Mono/Mø</dd><dt>SMC</dt><dd>SMC</dd><dt>EC</dt><dd>EC</dd><dt>T</dt><dd>T</dd><dt>NK</dt><dd>NK</dd><dt>B</dt><dd>B</dd><dt>PC</dt><dd>PC</dd><dt>MC</dt><dd>MC</dd></dl>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Mono/Mø'</li><li>'SMC'</li><li>'EC'</li><li>'T'</li><li>'NK'</li><li>'B'</li><li>'PC'</li><li>'MC'</li></ol>
</details>



```R
p3 <- Heatmap(mat_heatmap,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = cluster_info,
         top_annotation = top_anno,
         column_title = NULL,
       heatmap_legend_param = list(
          title = "zscore"
        ))
p3 <- ggplotify::as.ggplot(p3)
p3
```

    Warning message:
    “The input is a data frame-like object, convert it to a matrix.”



    
![png](Step10.4_Exclude_Heatmap_files/Step10.4_Exclude_Heatmap_25_1.png)
    



```R
ggsave(p3, file ='Heatmap_features.pdf')
```

    Saving 6.67 x 6.67 in image
    

