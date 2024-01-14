```R
############select brown module
rm(list = ls())
options(stringsAsFactors = F)
```


```R
library(tidyverse)
library(ggpubr)
library(DESeq2)
```

    Warning message in system("timedatectl", intern = TRUE):
    “running command 'timedatectl' had status 1”
    ── [1mAttaching packages[22m ─────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
    
    [32m✔[39m [34mggplot2[39m 3.4.0     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.5     [32m✔[39m [34mdplyr  [39m 1.0.7
    [32m✔[39m [34mtidyr  [39m 1.1.4     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 2.0.2     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    
    Loading required package: S4Vectors
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    Loading required package: parallel
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:parallel’:
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    
    The following objects are masked from ‘package:dplyr’:
    
        combine, intersect, setdiff, union
    
    
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
    
    
    The following objects are masked from ‘package:dplyr’:
    
        first, rename
    
    
    The following object is masked from ‘package:tidyr’:
    
        expand
    
    
    The following object is masked from ‘package:base’:
    
        expand.grid
    
    
    Loading required package: IRanges
    
    
    Attaching package: ‘IRanges’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        collapse, desc, slice
    
    
    The following object is masked from ‘package:purrr’:
    
        reduce
    
    
    Loading required package: GenomicRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following object is masked from ‘package:dplyr’:
    
        count
    
    
    
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
    
    



```R
######select top 25% brown gene 
####normalized data 
RNA_counts <- read_csv('./RNA_counts.csv')
```

    [1mRows: [22m[34m15871[39m [1mColumns: [22m[34m17[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m  (1): gene
    [32mdbl[39m (16): H_C_1, H_C_2, H_C_3, H_C_4, H_O_1, H_O_2, H_O_3, H_O_4, L_C_1, L_C...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
exp_dat <-  read_csv('./RNA_counts.csv') %>%
 column_to_rownames(var ='gene')
```

    [1mRows: [22m[34m15871[39m [1mColumns: [22m[34m17[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m  (1): gene
    [32mdbl[39m (16): H_C_1, H_C_2, H_C_3, H_C_4, H_O_1, H_O_2, H_O_3, H_O_4, L_C_1, L_C...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
exp_dat
```


<table class="dataframe">
<caption>A data.frame: 15871 × 16</caption>
<thead>
	<tr><th></th><th scope=col>H_C_1</th><th scope=col>H_C_2</th><th scope=col>H_C_3</th><th scope=col>H_C_4</th><th scope=col>H_O_1</th><th scope=col>H_O_2</th><th scope=col>H_O_3</th><th scope=col>H_O_4</th><th scope=col>L_C_1</th><th scope=col>L_C_2</th><th scope=col>L_C_3</th><th scope=col>L_C_4</th><th scope=col>L_O_1</th><th scope=col>L_O_2</th><th scope=col>L_O_3</th><th scope=col>L_O_4</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Gpnmb</th><td>856347</td><td>880695</td><td>988240</td><td>879217</td><td>734951</td><td>688516</td><td>698805</td><td>796416</td><td>804019</td><td>841302</td><td>798781</td><td>1223851</td><td>595370</td><td>712569</td><td>628444</td><td>663962</td></tr>
	<tr><th scope=row>Psap</th><td>757491</td><td>802603</td><td>889556</td><td>788284</td><td>777083</td><td>759529</td><td>789016</td><td>847598</td><td>803199</td><td>887335</td><td>816306</td><td> 934044</td><td>756695</td><td>861610</td><td>774337</td><td>817682</td></tr>
	<tr><th scope=row>Fth1</th><td>597268</td><td>600612</td><td>678216</td><td>599582</td><td>717214</td><td>702500</td><td>722328</td><td>760409</td><td>566725</td><td>622977</td><td>574328</td><td> 882613</td><td>632192</td><td>731840</td><td>644627</td><td>693134</td></tr>
	<tr><th scope=row>Ctsd</th><td>724282</td><td>742242</td><td>823936</td><td>729498</td><td>763092</td><td>737518</td><td>721393</td><td>795971</td><td>711425</td><td>751417</td><td>719019</td><td> 586377</td><td>670851</td><td>785205</td><td>696214</td><td>729278</td></tr>
	<tr><th scope=row>Ctsb</th><td>594753</td><td>622679</td><td>705216</td><td>617065</td><td>598325</td><td>595797</td><td>587809</td><td>657341</td><td>584321</td><td>635490</td><td>590337</td><td> 652491</td><td>515100</td><td>607682</td><td>538174</td><td>578670</td></tr>
	<tr><th scope=row>mt-Co1</th><td>273737</td><td>323889</td><td>348948</td><td>330516</td><td>355158</td><td>350056</td><td>345045</td><td>399833</td><td>260165</td><td>333136</td><td>290118</td><td> 703417</td><td>307928</td><td>335701</td><td>310949</td><td>337785</td></tr>
	<tr><th scope=row>Lyz2</th><td>483866</td><td>497642</td><td>588326</td><td>502520</td><td>593584</td><td>639652</td><td>574689</td><td>627841</td><td>464809</td><td>513088</td><td>448003</td><td> 571813</td><td>491525</td><td>548792</td><td>496298</td><td>524792</td></tr>
	<tr><th scope=row>Ctss</th><td>411338</td><td>421152</td><td>500370</td><td>428558</td><td>326509</td><td>391092</td><td>341348</td><td>351372</td><td>393340</td><td>451811</td><td>383562</td><td> 397366</td><td>262944</td><td>307433</td><td>267365</td><td>278393</td></tr>
	<tr><th scope=row>Mpeg1</th><td>368774</td><td>382868</td><td>438925</td><td>385075</td><td>310058</td><td>325328</td><td>295853</td><td>320687</td><td>323447</td><td>356596</td><td>311483</td><td> 334326</td><td>242072</td><td>282749</td><td>255083</td><td>270388</td></tr>
	<tr><th scope=row>Ftl1</th><td>289869</td><td>278015</td><td>313721</td><td>267817</td><td>352256</td><td>320765</td><td>370271</td><td>389492</td><td>254153</td><td>275232</td><td>265281</td><td> 199652</td><td>289308</td><td>360117</td><td>319179</td><td>334665</td></tr>
	<tr><th scope=row>Eef1a1</th><td>293979</td><td>306588</td><td>343862</td><td>305052</td><td>276501</td><td>270577</td><td>266888</td><td>287759</td><td>288259</td><td>328001</td><td>296657</td><td> 334017</td><td>249176</td><td>281537</td><td>253434</td><td>274041</td></tr>
	<tr><th scope=row>H2-D1</th><td>283335</td><td>285259</td><td>318880</td><td>278487</td><td>244785</td><td>263733</td><td>255186</td><td>261460</td><td>265467</td><td>293524</td><td>259833</td><td> 210670</td><td>209525</td><td>249274</td><td>218007</td><td>222180</td></tr>
	<tr><th scope=row>Apoe</th><td>280662</td><td>276820</td><td>309714</td><td>271198</td><td>246997</td><td>244755</td><td>256949</td><td>259950</td><td>205193</td><td>216631</td><td>193438</td><td> 244035</td><td>156612</td><td>194578</td><td>175086</td><td>173490</td></tr>
	<tr><th scope=row>B2m</th><td>135273</td><td>137548</td><td>166652</td><td>143974</td><td>107638</td><td>129926</td><td>114980</td><td>117056</td><td>141763</td><td>169492</td><td>145298</td><td> 264207</td><td>104724</td><td>115716</td><td>101954</td><td>108905</td></tr>
	<tr><th scope=row>Lrp1</th><td>177885</td><td>187158</td><td>207119</td><td>186177</td><td>185878</td><td>171981</td><td>180490</td><td>194928</td><td>172272</td><td>173108</td><td>164602</td><td> 190435</td><td>154204</td><td>165977</td><td>163347</td><td>175098</td></tr>
	<tr><th scope=row>Lamp1</th><td>136614</td><td>143361</td><td>159547</td><td>138047</td><td>141340</td><td>134688</td><td>133434</td><td>150047</td><td>136169</td><td>148988</td><td>137154</td><td> 197678</td><td>127509</td><td>143246</td><td>128712</td><td>133807</td></tr>
	<tr><th scope=row>Lgals3</th><td>169230</td><td>162609</td><td>180106</td><td>159225</td><td>169122</td><td>156664</td><td>172572</td><td>186985</td><td>156705</td><td>171373</td><td>168771</td><td> 162836</td><td>148642</td><td>178389</td><td>159713</td><td>159082</td></tr>
	<tr><th scope=row>Cd36</th><td> 50957</td><td> 52421</td><td> 58184</td><td> 54839</td><td>174616</td><td>146584</td><td>153625</td><td>179670</td><td> 56830</td><td> 64203</td><td> 67082</td><td>  62588</td><td>149632</td><td>176869</td><td>161680</td><td>171309</td></tr>
	<tr><th scope=row>mt-Nd5</th><td> 35117</td><td> 43260</td><td> 47187</td><td> 43010</td><td> 48696</td><td> 47634</td><td> 47568</td><td> 53550</td><td> 40970</td><td> 49335</td><td> 47574</td><td> 165596</td><td> 46928</td><td> 49272</td><td> 48242</td><td> 51688</td></tr>
	<tr><th scope=row>Cd74</th><td> 63240</td><td> 69060</td><td> 84860</td><td> 72755</td><td> 34268</td><td> 60211</td><td> 41118</td><td> 36102</td><td> 95460</td><td>124119</td><td> 91737</td><td> 163395</td><td> 35744</td><td> 37573</td><td> 29869</td><td> 30152</td></tr>
	<tr><th scope=row>Anpep</th><td> 91802</td><td> 94745</td><td>107559</td><td> 96798</td><td> 81361</td><td> 74166</td><td> 79434</td><td> 88893</td><td> 89960</td><td> 93915</td><td> 86383</td><td> 146733</td><td> 70678</td><td> 82900</td><td> 73553</td><td> 78049</td></tr>
	<tr><th scope=row>H2-K1</th><td>101692</td><td>105677</td><td>116270</td><td>103865</td><td> 79104</td><td> 98372</td><td> 82201</td><td> 82486</td><td>112290</td><td>132597</td><td>114176</td><td>  89688</td><td> 78285</td><td> 91675</td><td> 77545</td><td> 79115</td></tr>
	<tr><th scope=row>mt-Nd4</th><td> 31179</td><td> 47332</td><td> 37417</td><td> 48399</td><td> 37355</td><td> 38779</td><td> 43754</td><td> 59910</td><td> 29395</td><td> 50532</td><td> 37388</td><td> 124241</td><td> 45413</td><td> 50294</td><td> 44507</td><td> 49244</td></tr>
	<tr><th scope=row>Cybb</th><td> 86367</td><td> 94489</td><td>110307</td><td> 97920</td><td> 56925</td><td> 76049</td><td> 55145</td><td> 55479</td><td>102751</td><td>119456</td><td> 99329</td><td> 119558</td><td> 54689</td><td> 59286</td><td> 49802</td><td> 51744</td></tr>
	<tr><th scope=row>Cd68</th><td> 54690</td><td> 56856</td><td> 63846</td><td> 56886</td><td> 59155</td><td> 57481</td><td> 59007</td><td> 64065</td><td> 53632</td><td> 59625</td><td> 53432</td><td> 118633</td><td> 51358</td><td> 58424</td><td> 51728</td><td> 57107</td></tr>
	<tr><th scope=row>Pld3</th><td> 69050</td><td> 72002</td><td> 81778</td><td> 71797</td><td> 62344</td><td> 61569</td><td> 59520</td><td> 66118</td><td> 71220</td><td> 76726</td><td> 70734</td><td> 112232</td><td> 56105</td><td> 65753</td><td> 57815</td><td> 59861</td></tr>
	<tr><th scope=row>Laptm5</th><td> 59340</td><td> 60428</td><td> 66591</td><td> 60883</td><td> 63666</td><td> 59452</td><td> 62266</td><td> 66805</td><td> 59532</td><td> 63812</td><td> 58023</td><td> 107099</td><td> 58683</td><td> 65141</td><td> 58723</td><td> 60982</td></tr>
	<tr><th scope=row>mt-Nd1</th><td> 36150</td><td> 39977</td><td> 47213</td><td> 41104</td><td> 52676</td><td> 51430</td><td> 50534</td><td> 58476</td><td> 33872</td><td> 40001</td><td> 38306</td><td> 106797</td><td> 42452</td><td> 47353</td><td> 45427</td><td> 49693</td></tr>
	<tr><th scope=row>Vim</th><td> 77796</td><td> 83983</td><td> 92950</td><td> 83805</td><td> 99801</td><td> 96653</td><td> 95343</td><td>100978</td><td> 83569</td><td> 94841</td><td> 90315</td><td>  89998</td><td> 92373</td><td>106449</td><td> 88642</td><td> 96764</td></tr>
	<tr><th scope=row>Grn</th><td> 80332</td><td> 83385</td><td> 94522</td><td> 83547</td><td> 85384</td><td> 87143</td><td> 85045</td><td> 90492</td><td> 79737</td><td> 89793</td><td> 81930</td><td> 104143</td><td> 74184</td><td> 83601</td><td> 74819</td><td> 79185</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>Pnoc</th><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>2</td><td>2</td></tr>
	<tr><th scope=row>Mus_musculus_newGene_4398</th><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>2</td><td>2</td><td>2</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Olfr624</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>1</td><td>1</td></tr>
	<tr><th scope=row>Hist1h2ap</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>1</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>Gm15294</th><td>0</td><td>2</td><td>2</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>Gjb1</th><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>2</td><td>0</td><td>2</td><td>0</td><td>2</td></tr>
	<tr><th scope=row>Cldn24</th><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>1</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>G6pd2</th><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>Gm10271</th><td>2</td><td>0</td><td>1</td><td>1</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Fam205a2</th><td>0</td><td>2</td><td>1</td><td>1</td><td>0</td><td>0</td><td>2</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>4930544D05Rik</th><td>0</td><td>1</td><td>0</td><td>0</td><td>2</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Gm17622</th><td>0</td><td>1</td><td>0</td><td>2</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Actl7a</th><td>2</td><td>2</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Mus_musculus_newGene_1692</th><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>2</td><td>0</td><td>2</td></tr>
	<tr><th scope=row>Nphs2</th><td>2</td><td>2</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Gm42791</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>1</td><td>2</td><td>2</td></tr>
	<tr><th scope=row>Mus_musculus_newGene_375</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>2</td><td>1</td><td>2</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>AC140267.1</th><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td></tr>
	<tr><th scope=row>AL731706.2</th><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td></tr>
	<tr><th scope=row>Hist1h2ba</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Epo</th><td>1</td><td>1</td><td>1</td><td>0</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>0</td><td>1</td><td>1</td><td>1</td><td>1</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>Kif28</th><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>0</td><td>1</td><td>0</td><td>1</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>AL607142.1</th><td>1</td><td>0</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>0</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>Hist1h2ah</th><td>1</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td><td>1</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>Gm20538</th><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Gm11084</th><td>1</td><td>1</td><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Gm20662</th><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>Gm49342</th><td>0</td><td>1</td><td>1</td><td>0</td><td>1</td><td>1</td><td>0</td><td>1</td><td>1</td><td>0</td><td>0</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td></tr>
	<tr><th scope=row>Exosc6</th><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td><td>1</td></tr>
	<tr><th scope=row>Gm28053</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>1</td><td>0</td><td>0</td><td>0</td><td>1</td><td>1</td><td>1</td><td>0</td></tr>
</tbody>
</table>




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
###select interested genes
foam_gene <- c('Gdf15','Aldoa','Creg1','Lgmn','Lilr4b','Pkm','Plin2','Taldo1','Tnfaip2')
```


```R
###exclude abnormal samples and select brown top 10%
###prepare RF_input_for_r 
RF_input <- assay(dds_norm) %>%
.[,c(1:11,13:16)] %>%
.[foam_gene,] %>% 
t() %>%
as.data.frame() %>%
mutate(group = rep(c('H_C','H_O','L_C','L_O'), c(4,4,3,4)))
head(RF_input)
###rename - by t
colnames(RF_input) <- str_replace_all(colnames(RF_input),'-','_')
```


<table class="dataframe">
<caption>A data.frame: 6 × 10</caption>
<thead>
	<tr><th></th><th scope=col>Gdf15</th><th scope=col>Aldoa</th><th scope=col>Creg1</th><th scope=col>Lgmn</th><th scope=col>Lilr4b</th><th scope=col>Pkm</th><th scope=col>Plin2</th><th scope=col>Taldo1</th><th scope=col>Tnfaip2</th><th scope=col>group</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>H_C_1</th><td>10.98536</td><td>14.93508</td><td>15.52910</td><td>16.19866</td><td>14.58548</td><td>15.44006</td><td>14.95568</td><td>14.02696</td><td>14.43839</td><td>H_C</td></tr>
	<tr><th scope=row>H_C_2</th><td>10.95770</td><td>14.88724</td><td>15.43755</td><td>16.21587</td><td>14.66305</td><td>15.32849</td><td>14.96179</td><td>13.91388</td><td>14.45043</td><td>H_C</td></tr>
	<tr><th scope=row>H_C_3</th><td>10.77258</td><td>14.85852</td><td>15.48585</td><td>16.18843</td><td>14.65579</td><td>15.37205</td><td>14.92798</td><td>13.87415</td><td>14.37119</td><td>H_C</td></tr>
	<tr><th scope=row>H_C_4</th><td>10.85834</td><td>14.86347</td><td>15.47195</td><td>16.18000</td><td>14.66837</td><td>15.36304</td><td>14.96870</td><td>13.86905</td><td>14.46267</td><td>H_C</td></tr>
	<tr><th scope=row>H_O_1</th><td>11.66893</td><td>15.03698</td><td>16.01711</td><td>16.34234</td><td>14.96809</td><td>15.91897</td><td>15.76242</td><td>14.45946</td><td>14.58633</td><td>H_O</td></tr>
	<tr><th scope=row>H_O_2</th><td>11.51231</td><td>15.03044</td><td>16.06649</td><td>16.37537</td><td>14.98776</td><td>15.98559</td><td>15.66953</td><td>14.35520</td><td>14.58289</td><td>H_O</td></tr>
</tbody>
</table>




```R
plot_gene <- function(x){
p <- RF_input[,c(x,'group')] %>%
mutate(group = factor(group, levels = c('L_C','L_O','H_C','H_O'))) %>%
magrittr::set_colnames(c('gene','group')) %>%
ggboxplot(title = x, ylab = '', xlab = '',
          x = 'group', y = 'gene',fill  = 'group', add = 'jitter', palette = 'lancet',legend = "none") +
stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('L_O', 'H_O'),c('H_C','H_O')),label = 'p.signif') +
  theme(plot.title=element_text(hjust=0.5, size = 8),
         axis.text.x = element_text(angle=0,hjust=1,size=6),
         axis.text.y = element_text(size=6),
         axis.title.y = element_text(size=6),
          legend.title =  element_text(size = 6),
           legend.text = element_text(size = 6))
return(p)}
```


```R
colnames(RF_input)
p1 <- plot_gene(colnames(RF_input)[1])
p2 <- plot_gene(colnames(RF_input)[2])
p3 <- plot_gene(colnames(RF_input)[3])
p4 <- plot_gene(colnames(RF_input)[4])
p5 <- plot_gene(colnames(RF_input)[5])
p6 <- plot_gene(colnames(RF_input)[6])
p7 <- plot_gene(colnames(RF_input)[7])
p8 <- plot_gene(colnames(RF_input)[8])
p9 <- plot_gene(colnames(RF_input)[9])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Gdf15'</li><li>'Aldoa'</li><li>'Creg1'</li><li>'Lgmn'</li><li>'Lilr4b'</li><li>'Pkm'</li><li>'Plin2'</li><li>'Taldo1'</li><li>'Tnfaip2'</li><li>'group'</li></ol>




```R
library(patchwork)
p0 <- p2 + p3 + p4 + p6+ p5 +p1 + p7 + p8 + p9 + plot_layout(ncol = 3)
p0
ggsave(p0, file = './Final_Results/Figure3/Selected_gene.pdf')
```

    [1m[22mSaving 6.67 x 6.67 in image



    
![png](Step9_SelectGene_Group_Compare_files/Step9_SelectGene_Group_Compare_13_1.png)
    



```R
RF_input %>%
pivot_longer(-group, names_to = 'gene', values_to = 'count') %>%
mutate(group = factor(group, levels = c('L_C','L_O','H_C','H_O'))) %>%
ggboxplot( ylab = 'count', xlab = '', x = 'group', y = 'count', 
          fill  = 'group', add = 'jitter', palette = 'lancet', 
          legend = "right")+
stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('L_O', 'H_O'),c('H_C','H_O')),label = 'p.signif') +
 facet_wrap(facets = 'gene', scales = "free", ncol = 5)
```


    
![png](Step9_SelectGene_Group_Compare_files/Step9_SelectGene_Group_Compare_14_0.png)
    

