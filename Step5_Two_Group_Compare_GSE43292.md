```R
library(tidyverse)
library(GSVA)
library(GEOquery) 
library(limma)
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
    
    Loading required package: Biobase
    
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
    
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Setting options('download.file.method.GEOquery'='auto')
    
    Setting options('GEOquery.inmemory.gpl'=FALSE)
    
    
    Attaching package: ‘limma’
    
    
    The following object is masked from ‘package:BiocGenerics’:
    
        plotMA
    
    



```R
gset  <-  GEOquery::getGEO('GSE43292',getGPL = F)
```

    Found 1 file(s)
    
    GSE43292_series_matrix.txt.gz
    
    [1mRows: [22m[34m33297[39m [1mColumns: [22m[34m65[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m "\t"
    [32mdbl[39m (65): ID_REF, GSM1060117, GSM1060118, GSM1060119, GSM1060120, GSM1060121...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
normalized_gset  <- gset$GSE43292_series_matrix.txt.gz@assayData$exprs
```


```R
range(normalized_gset)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>1.70657</li><li>14.3077</li></ol>




```R
normalized_gset
```


<table class="dataframe">
<caption>A matrix: 33297 × 64 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>GSM1060117</th><th scope=col>GSM1060118</th><th scope=col>GSM1060119</th><th scope=col>GSM1060120</th><th scope=col>GSM1060121</th><th scope=col>GSM1060122</th><th scope=col>GSM1060123</th><th scope=col>GSM1060124</th><th scope=col>GSM1060125</th><th scope=col>GSM1060126</th><th scope=col>⋯</th><th scope=col>GSM1060171</th><th scope=col>GSM1060172</th><th scope=col>GSM1060173</th><th scope=col>GSM1060174</th><th scope=col>GSM1060175</th><th scope=col>GSM1060176</th><th scope=col>GSM1060177</th><th scope=col>GSM1060178</th><th scope=col>GSM1060179</th><th scope=col>GSM1060180</th></tr>
</thead>
<tbody>
	<tr><th scope=row>7892501</th><td> 3.78120</td><td> 3.14761</td><td> 4.45896</td><td> 4.02410</td><td> 4.86666</td><td> 5.52737</td><td> 5.13115</td><td> 4.68490</td><td> 5.33053</td><td> 5.21742</td><td>⋯</td><td> 3.99150</td><td> 4.99254</td><td> 6.80393</td><td> 4.07402</td><td> 5.86497</td><td> 5.91446</td><td> 5.35851</td><td> 4.94704</td><td> 6.05089</td><td> 3.09624</td></tr>
	<tr><th scope=row>7892502</th><td> 4.44462</td><td> 4.32176</td><td> 5.30920</td><td> 5.40423</td><td> 5.15199</td><td> 5.16896</td><td> 5.19199</td><td> 5.52297</td><td> 4.50558</td><td> 4.65962</td><td>⋯</td><td> 5.37767</td><td> 5.34172</td><td> 5.10396</td><td> 5.29750</td><td> 4.54258</td><td> 4.67422</td><td> 5.81923</td><td> 4.49336</td><td> 4.65179</td><td> 4.62726</td></tr>
	<tr><th scope=row>7892503</th><td> 3.35557</td><td> 2.99685</td><td> 4.28775</td><td> 3.96073</td><td> 4.36725</td><td> 3.51501</td><td> 4.23388</td><td> 4.23024</td><td> 3.25248</td><td> 4.60459</td><td>⋯</td><td> 3.65493</td><td> 4.42915</td><td> 5.12905</td><td> 4.35436</td><td> 3.90701</td><td> 4.52855</td><td> 3.68901</td><td> 4.85516</td><td> 4.45861</td><td> 3.11494</td></tr>
	<tr><th scope=row>7892504</th><td>10.77280</td><td>10.39450</td><td> 9.92374</td><td> 9.90876</td><td>10.64390</td><td>10.32270</td><td>10.65320</td><td>10.19200</td><td>10.47110</td><td>10.11110</td><td>⋯</td><td> 9.56718</td><td> 9.29764</td><td> 9.77111</td><td> 9.41572</td><td> 9.78930</td><td> 9.85026</td><td> 9.79740</td><td> 9.58430</td><td> 9.93696</td><td> 9.85777</td></tr>
	<tr><th scope=row>7892505</th><td> 3.52998</td><td> 3.69160</td><td> 2.93842</td><td> 2.58252</td><td> 3.38546</td><td> 3.43945</td><td> 2.80291</td><td> 3.06008</td><td> 3.36183</td><td> 3.68809</td><td>⋯</td><td> 3.14876</td><td> 3.54822</td><td> 3.36739</td><td> 2.77438</td><td> 3.11674</td><td> 3.33899</td><td> 3.11424</td><td> 3.31875</td><td> 3.10718</td><td> 3.78139</td></tr>
	<tr><th scope=row>7892506</th><td> 3.53808</td><td> 3.63492</td><td> 4.47884</td><td> 4.21690</td><td> 3.36820</td><td> 4.02052</td><td> 3.76121</td><td> 4.51591</td><td> 3.31265</td><td> 4.10651</td><td>⋯</td><td> 3.93250</td><td> 3.73839</td><td> 4.60891</td><td> 4.63011</td><td> 4.06505</td><td> 4.01099</td><td> 3.79005</td><td> 3.88709</td><td> 3.73716</td><td> 3.76657</td></tr>
	<tr><th scope=row>7892507</th><td> 4.90025</td><td> 3.88422</td><td> 4.68580</td><td> 4.92572</td><td> 4.89072</td><td> 4.67694</td><td> 4.37426</td><td> 4.06716</td><td> 4.17041</td><td> 4.18189</td><td>⋯</td><td> 4.63240</td><td> 4.74839</td><td> 4.86914</td><td> 4.47295</td><td> 4.85909</td><td> 5.21971</td><td> 4.33575</td><td> 4.97125</td><td> 5.42417</td><td> 5.59303</td></tr>
	<tr><th scope=row>7892508</th><td> 3.12491</td><td> 4.60145</td><td> 6.45946</td><td> 6.18453</td><td> 5.71111</td><td> 5.88986</td><td> 4.99173</td><td> 6.46703</td><td> 4.51229</td><td> 5.46217</td><td>⋯</td><td> 6.16781</td><td> 5.56687</td><td> 4.87189</td><td> 5.21983</td><td> 5.45599</td><td> 5.59809</td><td> 5.17229</td><td> 4.81833</td><td> 5.14684</td><td> 4.94791</td></tr>
	<tr><th scope=row>7892509</th><td>10.98500</td><td>11.07720</td><td>11.09760</td><td>11.21710</td><td>11.45750</td><td>11.23830</td><td>11.51570</td><td>11.14390</td><td>11.21840</td><td>11.38350</td><td>⋯</td><td>11.05730</td><td>11.00190</td><td>11.49150</td><td>11.20460</td><td>11.51930</td><td>11.23650</td><td>11.28410</td><td>10.77650</td><td>11.26450</td><td>10.55250</td></tr>
	<tr><th scope=row>7892510</th><td> 4.59634</td><td> 4.47389</td><td> 5.70508</td><td> 5.33982</td><td> 4.10943</td><td> 5.24095</td><td> 4.56018</td><td> 5.67107</td><td> 4.09623</td><td> 4.93858</td><td>⋯</td><td> 5.78777</td><td> 3.70906</td><td> 3.61431</td><td> 3.58086</td><td> 4.49461</td><td> 4.24192</td><td> 4.97155</td><td> 5.32210</td><td> 3.82910</td><td> 3.80811</td></tr>
	<tr><th scope=row>7892511</th><td> 3.58720</td><td> 3.32765</td><td> 3.34008</td><td> 3.77569</td><td> 3.42512</td><td> 3.54140</td><td> 3.34277</td><td> 3.02534</td><td> 2.70667</td><td> 2.92360</td><td>⋯</td><td> 3.56567</td><td> 3.95456</td><td> 3.39193</td><td> 3.42260</td><td> 2.69664</td><td> 2.96452</td><td> 3.29143</td><td> 3.26436</td><td> 2.95017</td><td> 3.32259</td></tr>
	<tr><th scope=row>7892512</th><td> 7.32139</td><td> 6.76761</td><td> 7.25412</td><td> 7.44505</td><td> 7.38296</td><td> 7.33182</td><td> 7.49915</td><td> 7.38601</td><td> 7.25348</td><td> 7.41457</td><td>⋯</td><td> 7.36474</td><td> 7.44655</td><td> 7.20203</td><td> 7.02896</td><td> 7.62207</td><td> 7.04899</td><td> 6.71407</td><td> 6.84409</td><td> 7.72382</td><td> 7.66248</td></tr>
	<tr><th scope=row>7892513</th><td> 3.98102</td><td> 3.75270</td><td> 4.38208</td><td> 4.62070</td><td> 4.10084</td><td> 3.50630</td><td> 3.38108</td><td> 3.73943</td><td> 3.99348</td><td> 3.88434</td><td>⋯</td><td> 4.05522</td><td> 4.50165</td><td> 4.39227</td><td> 5.04413</td><td> 4.21715</td><td> 3.48366</td><td> 4.84921</td><td> 3.85346</td><td> 3.63921</td><td> 3.66646</td></tr>
	<tr><th scope=row>7892514</th><td>10.56670</td><td>10.32720</td><td>10.26690</td><td>10.16720</td><td> 9.46315</td><td>10.45450</td><td> 9.93993</td><td>10.32290</td><td>10.10150</td><td>10.46820</td><td>⋯</td><td>11.00670</td><td>11.21430</td><td>11.19890</td><td>10.83750</td><td>10.99960</td><td>10.93890</td><td>10.96410</td><td>10.82480</td><td>10.87610</td><td>10.94650</td></tr>
	<tr><th scope=row>7892515</th><td> 9.92836</td><td> 9.41894</td><td> 9.32436</td><td> 9.41651</td><td> 9.30495</td><td> 9.56982</td><td> 9.48037</td><td> 9.33006</td><td> 9.28653</td><td> 9.59452</td><td>⋯</td><td> 9.12917</td><td> 8.57535</td><td> 9.14983</td><td> 9.14164</td><td> 9.07543</td><td> 9.26962</td><td> 9.12941</td><td> 9.19487</td><td> 8.92107</td><td> 9.01572</td></tr>
	<tr><th scope=row>7892516</th><td> 4.04393</td><td> 3.62371</td><td> 4.58426</td><td> 3.53703</td><td> 3.54068</td><td> 3.93546</td><td> 3.73194</td><td> 4.55927</td><td> 3.76566</td><td> 4.45994</td><td>⋯</td><td> 6.47003</td><td> 6.05982</td><td> 5.88598</td><td> 5.44331</td><td> 4.73576</td><td> 5.54231</td><td> 5.61045</td><td> 4.76609</td><td> 5.42525</td><td> 5.29021</td></tr>
	<tr><th scope=row>7892517</th><td> 4.96047</td><td> 5.64362</td><td> 6.99294</td><td> 6.90937</td><td> 6.15123</td><td> 6.30916</td><td> 5.73223</td><td> 6.33918</td><td> 6.11746</td><td> 6.19448</td><td>⋯</td><td> 6.87188</td><td> 6.20869</td><td> 6.52633</td><td> 6.39396</td><td> 5.72065</td><td> 6.23042</td><td> 5.68169</td><td> 7.02973</td><td> 6.34240</td><td> 6.82318</td></tr>
	<tr><th scope=row>7892518</th><td> 2.91881</td><td> 3.08391</td><td> 3.03377</td><td> 2.46635</td><td> 3.08841</td><td> 3.05242</td><td> 3.20814</td><td> 2.70778</td><td> 2.94227</td><td> 2.74480</td><td>⋯</td><td> 2.81872</td><td> 3.11844</td><td> 2.93278</td><td> 2.67548</td><td> 2.78379</td><td> 3.31550</td><td> 3.07297</td><td> 2.85913</td><td> 3.50342</td><td> 3.88708</td></tr>
	<tr><th scope=row>7892519</th><td> 4.81200</td><td> 4.75262</td><td> 4.72547</td><td> 4.87091</td><td> 4.77413</td><td> 5.12273</td><td> 5.01338</td><td> 4.69190</td><td> 4.51101</td><td> 4.43525</td><td>⋯</td><td> 4.79111</td><td> 4.62661</td><td> 3.71532</td><td> 4.21932</td><td> 4.78999</td><td> 4.72737</td><td> 4.30614</td><td> 4.60119</td><td> 4.74467</td><td> 5.19458</td></tr>
	<tr><th scope=row>7892520</th><td> 9.16042</td><td> 9.34007</td><td> 9.59116</td><td> 9.66963</td><td> 9.46023</td><td> 9.63462</td><td> 9.54866</td><td> 9.61916</td><td> 9.41118</td><td> 9.54226</td><td>⋯</td><td> 9.75179</td><td> 9.49257</td><td> 9.80800</td><td> 9.42165</td><td> 9.49875</td><td> 9.68014</td><td> 9.08274</td><td> 9.71196</td><td> 9.69190</td><td> 9.44335</td></tr>
	<tr><th scope=row>7892521</th><td> 6.26117</td><td> 6.21798</td><td> 5.55011</td><td> 6.34510</td><td> 6.12019</td><td> 6.65748</td><td> 5.93410</td><td> 6.57036</td><td> 5.77405</td><td> 6.28758</td><td>⋯</td><td> 6.68947</td><td> 6.41834</td><td> 6.30424</td><td> 6.65764</td><td> 6.48381</td><td> 6.85972</td><td> 6.84596</td><td> 6.28956</td><td> 6.48117</td><td> 6.78590</td></tr>
	<tr><th scope=row>7892522</th><td> 7.24132</td><td> 6.50037</td><td> 6.70237</td><td> 7.09848</td><td> 6.54655</td><td> 6.70077</td><td> 6.32955</td><td> 6.97728</td><td> 6.54321</td><td> 6.48605</td><td>⋯</td><td> 6.49853</td><td> 6.98439</td><td> 6.84969</td><td> 6.44543</td><td> 7.01268</td><td> 6.90304</td><td> 6.64695</td><td> 6.74607</td><td> 6.41080</td><td> 6.71614</td></tr>
	<tr><th scope=row>7892524</th><td> 6.50475</td><td> 6.36394</td><td> 6.25285</td><td> 6.43161</td><td> 7.06881</td><td> 7.35106</td><td> 6.97076</td><td> 7.24345</td><td> 6.47792</td><td> 7.07710</td><td>⋯</td><td> 7.40984</td><td> 7.00005</td><td> 6.96628</td><td> 6.78054</td><td> 7.31059</td><td> 6.90798</td><td> 7.21225</td><td> 6.60862</td><td> 6.78976</td><td> 7.19713</td></tr>
	<tr><th scope=row>7892525</th><td> 5.83484</td><td> 5.64223</td><td> 6.02651</td><td> 5.84238</td><td> 5.51708</td><td> 6.37616</td><td> 6.10214</td><td> 6.51659</td><td> 6.41889</td><td> 6.12325</td><td>⋯</td><td> 6.34923</td><td> 6.30351</td><td> 5.99033</td><td> 6.19882</td><td> 6.01760</td><td> 6.11311</td><td> 5.58630</td><td> 6.20321</td><td> 6.18283</td><td> 6.46917</td></tr>
	<tr><th scope=row>7892526</th><td> 5.16558</td><td> 4.55290</td><td> 5.05065</td><td> 4.94202</td><td> 4.71531</td><td> 5.52206</td><td> 4.63619</td><td> 5.08166</td><td> 4.43056</td><td> 4.87026</td><td>⋯</td><td> 4.95583</td><td> 5.10848</td><td> 5.71482</td><td> 5.25144</td><td> 5.34415</td><td> 5.57969</td><td> 5.31767</td><td> 4.73032</td><td> 5.25397</td><td> 5.10265</td></tr>
	<tr><th scope=row>7892527</th><td> 8.19672</td><td> 8.28713</td><td> 7.97270</td><td> 8.41051</td><td> 8.34531</td><td> 8.61247</td><td> 8.34169</td><td> 8.89471</td><td> 8.07311</td><td> 8.66033</td><td>⋯</td><td> 8.94087</td><td> 8.97108</td><td> 8.32290</td><td> 8.62013</td><td> 8.31679</td><td> 8.73433</td><td> 8.69310</td><td> 8.25254</td><td> 8.44265</td><td> 8.66906</td></tr>
	<tr><th scope=row>7892528</th><td> 2.64652</td><td> 2.44573</td><td> 2.61131</td><td> 2.72695</td><td> 2.66193</td><td> 2.63984</td><td> 2.48122</td><td> 2.39733</td><td> 2.52096</td><td> 2.87177</td><td>⋯</td><td> 2.11999</td><td> 2.60319</td><td> 2.61411</td><td> 2.64521</td><td> 2.71584</td><td> 2.16796</td><td> 2.52625</td><td> 2.27893</td><td> 2.71995</td><td> 2.66174</td></tr>
	<tr><th scope=row>7892529</th><td> 6.50272</td><td> 5.91494</td><td> 6.33529</td><td> 6.47753</td><td> 5.99535</td><td> 6.36200</td><td> 6.45720</td><td> 6.07984</td><td> 5.67637</td><td> 6.62035</td><td>⋯</td><td> 6.81971</td><td> 5.63092</td><td> 6.24533</td><td> 6.42928</td><td> 6.43303</td><td> 6.32516</td><td> 5.39851</td><td> 5.94857</td><td> 5.86997</td><td> 6.36188</td></tr>
	<tr><th scope=row>7892530</th><td> 8.32169</td><td> 8.51596</td><td> 8.94346</td><td> 8.50791</td><td> 9.03597</td><td> 9.02741</td><td> 9.25164</td><td> 8.71113</td><td> 8.83785</td><td> 8.93800</td><td>⋯</td><td>10.10480</td><td> 9.17296</td><td>10.12940</td><td> 9.63970</td><td> 9.86072</td><td> 9.61936</td><td>10.49580</td><td> 9.16675</td><td> 9.57942</td><td> 9.23620</td></tr>
	<tr><th scope=row>7892531</th><td> 6.70481</td><td> 6.63265</td><td> 7.23591</td><td> 7.01465</td><td> 6.14176</td><td> 6.70057</td><td> 7.02832</td><td> 6.09762</td><td> 6.06085</td><td> 6.53617</td><td>⋯</td><td> 7.50582</td><td> 6.88529</td><td> 7.56922</td><td> 7.52210</td><td> 7.40056</td><td> 7.40839</td><td> 7.33634</td><td> 6.99176</td><td> 7.06963</td><td> 7.27906</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>8180389</th><td>10.25940</td><td>10.48780</td><td> 9.50512</td><td> 9.52102</td><td> 9.75484</td><td> 9.97529</td><td> 9.52479</td><td>10.10150</td><td>10.32440</td><td>10.04660</td><td>⋯</td><td>10.11880</td><td>10.25080</td><td>10.21280</td><td>10.16020</td><td> 9.99636</td><td> 9.99221</td><td>10.64210</td><td>10.23000</td><td> 9.93136</td><td> 9.98756</td></tr>
	<tr><th scope=row>8180390</th><td> 6.55343</td><td> 6.94296</td><td> 7.25211</td><td> 7.51746</td><td> 7.33980</td><td> 7.42288</td><td> 7.39122</td><td> 7.37606</td><td> 7.24983</td><td> 7.36211</td><td>⋯</td><td> 7.18293</td><td> 7.30927</td><td> 7.06765</td><td> 6.91456</td><td> 7.39217</td><td> 7.27983</td><td> 6.60256</td><td> 6.91987</td><td> 7.16131</td><td> 7.02385</td></tr>
	<tr><th scope=row>8180391</th><td> 6.27752</td><td> 6.41082</td><td> 6.10720</td><td> 6.32630</td><td> 6.16498</td><td> 6.22538</td><td> 6.24410</td><td> 6.19746</td><td> 6.24531</td><td> 6.08057</td><td>⋯</td><td> 6.11179</td><td> 5.96711</td><td> 5.92602</td><td> 6.00924</td><td> 6.04623</td><td> 6.11470</td><td> 5.95033</td><td> 6.34734</td><td> 6.37751</td><td> 6.01595</td></tr>
	<tr><th scope=row>8180392</th><td> 7.47708</td><td> 7.71490</td><td> 7.28661</td><td> 7.10324</td><td> 7.33411</td><td> 7.19557</td><td> 7.24670</td><td> 7.31382</td><td> 7.42553</td><td> 7.49212</td><td>⋯</td><td> 7.96154</td><td> 7.58768</td><td> 7.64246</td><td> 7.59070</td><td> 7.37142</td><td> 7.41636</td><td> 8.42176</td><td> 7.53810</td><td> 7.18542</td><td> 7.47624</td></tr>
	<tr><th scope=row>8180393</th><td> 9.62144</td><td> 9.45111</td><td>10.00760</td><td> 9.89795</td><td>10.17240</td><td> 9.76127</td><td>10.13490</td><td> 9.57931</td><td> 9.56152</td><td> 9.66254</td><td>⋯</td><td> 9.32376</td><td> 9.25628</td><td> 9.93554</td><td> 9.76309</td><td> 9.92998</td><td> 9.61374</td><td> 9.85508</td><td> 8.80390</td><td> 9.62855</td><td> 9.34498</td></tr>
	<tr><th scope=row>8180394</th><td> 9.68172</td><td> 9.53056</td><td>10.09810</td><td> 9.98183</td><td>10.08600</td><td> 9.74550</td><td>10.13510</td><td> 9.59952</td><td> 9.54534</td><td> 9.73178</td><td>⋯</td><td> 9.36752</td><td> 9.22824</td><td> 9.83735</td><td> 9.71157</td><td> 9.87270</td><td> 9.62582</td><td> 9.85766</td><td> 8.68072</td><td> 9.55136</td><td> 9.33808</td></tr>
	<tr><th scope=row>8180395</th><td> 9.39315</td><td> 9.27998</td><td> 9.73331</td><td> 9.59819</td><td> 9.96867</td><td> 9.52534</td><td> 9.88094</td><td> 9.23278</td><td> 9.39605</td><td> 9.44490</td><td>⋯</td><td> 9.10917</td><td> 9.08255</td><td> 9.79733</td><td> 9.59374</td><td> 9.77394</td><td> 9.42014</td><td> 9.66356</td><td> 9.05709</td><td> 9.48309</td><td> 9.25779</td></tr>
	<tr><th scope=row>8180396</th><td> 7.29822</td><td> 7.32523</td><td> 7.52480</td><td> 7.50253</td><td> 7.39994</td><td> 7.50759</td><td> 7.41806</td><td> 7.37714</td><td> 7.26622</td><td> 7.42727</td><td>⋯</td><td> 7.76933</td><td> 7.71499</td><td> 7.72029</td><td> 7.63131</td><td> 7.57446</td><td> 7.88152</td><td> 8.04252</td><td> 7.25821</td><td> 7.51894</td><td> 7.92866</td></tr>
	<tr><th scope=row>8180397</th><td> 7.11367</td><td> 7.03657</td><td> 7.35304</td><td> 7.31792</td><td> 7.21629</td><td> 7.33181</td><td> 7.19340</td><td> 7.19468</td><td> 7.07087</td><td> 7.35402</td><td>⋯</td><td> 7.58476</td><td> 7.52704</td><td> 7.53266</td><td> 7.48995</td><td> 7.32333</td><td> 7.69847</td><td> 7.85365</td><td> 6.85940</td><td> 7.30787</td><td> 7.91634</td></tr>
	<tr><th scope=row>8180398</th><td>10.97400</td><td>11.03890</td><td>10.30640</td><td>10.19620</td><td>10.83920</td><td>11.05650</td><td>10.42290</td><td>10.80710</td><td>10.94450</td><td>10.72370</td><td>⋯</td><td>10.71350</td><td>10.59840</td><td>10.75180</td><td>10.61090</td><td>10.79080</td><td>10.73600</td><td>10.98840</td><td>10.69370</td><td>10.64030</td><td>10.64270</td></tr>
	<tr><th scope=row>8180399</th><td> 4.25999</td><td> 4.17867</td><td> 4.31537</td><td> 4.44603</td><td> 4.06956</td><td> 4.40342</td><td> 4.15644</td><td> 4.20676</td><td> 3.93901</td><td> 4.43557</td><td>⋯</td><td> 4.22084</td><td> 4.04266</td><td> 3.95714</td><td> 4.02875</td><td> 4.29990</td><td> 4.08211</td><td> 3.94899</td><td> 4.25583</td><td> 4.05437</td><td> 4.17966</td></tr>
	<tr><th scope=row>8180400</th><td> 3.65236</td><td> 4.06900</td><td> 4.03755</td><td> 3.21823</td><td> 3.88289</td><td> 4.01905</td><td> 3.51342</td><td> 3.51326</td><td> 4.05536</td><td> 3.21948</td><td>⋯</td><td> 3.41229</td><td> 4.23570</td><td> 3.80310</td><td> 3.49330</td><td> 3.29536</td><td> 3.21670</td><td> 3.61762</td><td> 3.32961</td><td> 3.56860</td><td> 3.68296</td></tr>
	<tr><th scope=row>8180401</th><td> 4.83523</td><td> 4.88051</td><td> 4.94395</td><td> 4.82221</td><td> 4.89523</td><td> 4.70820</td><td> 5.32654</td><td> 4.77466</td><td> 4.83681</td><td> 4.83352</td><td>⋯</td><td> 4.73938</td><td> 4.76251</td><td> 4.68211</td><td> 4.69493</td><td> 5.03794</td><td> 4.83336</td><td> 4.71328</td><td> 4.64660</td><td> 4.80849</td><td> 4.74824</td></tr>
	<tr><th scope=row>8180402</th><td> 9.70755</td><td> 9.86928</td><td>10.11720</td><td>10.20820</td><td> 9.75277</td><td> 9.89166</td><td> 9.53704</td><td> 9.78637</td><td> 9.76651</td><td> 9.75701</td><td>⋯</td><td> 9.88208</td><td>10.03230</td><td> 9.68623</td><td> 9.60079</td><td> 9.55866</td><td> 9.72917</td><td> 9.64189</td><td> 9.87875</td><td> 9.78194</td><td>10.02910</td></tr>
	<tr><th scope=row>8180403</th><td> 8.57627</td><td> 8.44579</td><td> 8.91200</td><td> 8.91013</td><td> 9.37731</td><td> 8.67250</td><td> 9.35507</td><td> 8.69271</td><td> 8.80034</td><td> 8.90616</td><td>⋯</td><td> 8.33870</td><td> 8.41106</td><td> 8.75546</td><td> 8.88114</td><td> 8.91929</td><td> 8.49891</td><td> 8.29798</td><td> 8.33779</td><td> 8.56979</td><td> 8.29840</td></tr>
	<tr><th scope=row>8180404</th><td>11.03150</td><td>11.05740</td><td>11.00030</td><td>10.79400</td><td>10.90360</td><td>11.03720</td><td>10.58610</td><td>10.76890</td><td>10.92230</td><td>10.71060</td><td>⋯</td><td>10.67790</td><td>10.54910</td><td>10.86370</td><td>10.81430</td><td>10.85990</td><td>10.68160</td><td>11.11630</td><td>10.91790</td><td>10.83590</td><td>10.78740</td></tr>
	<tr><th scope=row>8180405</th><td> 7.06736</td><td> 7.06660</td><td> 7.06434</td><td> 7.05768</td><td> 6.92040</td><td> 7.32041</td><td> 6.69195</td><td> 7.07610</td><td> 6.34975</td><td> 6.91917</td><td>⋯</td><td> 7.24386</td><td> 7.25713</td><td> 6.76863</td><td> 6.59446</td><td> 6.75317</td><td> 6.88795</td><td> 7.07286</td><td> 7.01567</td><td> 6.44596</td><td> 7.00449</td></tr>
	<tr><th scope=row>8180406</th><td> 6.88872</td><td> 6.23611</td><td> 6.21732</td><td> 6.31469</td><td> 6.33552</td><td> 6.63640</td><td> 6.38918</td><td> 6.42896</td><td> 6.69006</td><td> 6.37840</td><td>⋯</td><td> 6.28070</td><td> 6.40909</td><td> 6.31329</td><td> 6.43169</td><td> 6.27090</td><td> 6.43424</td><td> 6.62205</td><td> 6.94218</td><td> 6.57735</td><td> 6.81541</td></tr>
	<tr><th scope=row>8180407</th><td>10.96870</td><td>11.04100</td><td>11.00200</td><td>10.92260</td><td>11.30090</td><td>11.33560</td><td>11.38710</td><td>11.17690</td><td>11.07330</td><td>11.17450</td><td>⋯</td><td>10.99780</td><td>11.00220</td><td>11.25560</td><td>11.05570</td><td>11.18210</td><td>11.13550</td><td>11.14100</td><td>10.82810</td><td>11.15740</td><td>10.88130</td></tr>
	<tr><th scope=row>8180408</th><td>11.05830</td><td>11.09920</td><td>11.08440</td><td>10.99330</td><td>11.34160</td><td>11.39000</td><td>11.44640</td><td>11.23050</td><td>11.14290</td><td>11.23050</td><td>⋯</td><td>11.09070</td><td>11.05300</td><td>11.27280</td><td>11.11970</td><td>11.22700</td><td>11.18270</td><td>11.18610</td><td>10.87430</td><td>11.20130</td><td>10.92530</td></tr>
	<tr><th scope=row>8180409</th><td>11.38070</td><td>11.43320</td><td>11.35840</td><td>11.26210</td><td>11.61950</td><td>11.65860</td><td>11.68930</td><td>11.49070</td><td>11.40840</td><td>11.49980</td><td>⋯</td><td>11.24730</td><td>11.21450</td><td>11.39920</td><td>11.26170</td><td>11.40340</td><td>11.37870</td><td>11.35600</td><td>11.04900</td><td>11.31690</td><td>11.09340</td></tr>
	<tr><th scope=row>8180410</th><td> 9.86276</td><td> 9.77370</td><td> 9.85765</td><td> 9.76752</td><td> 9.96494</td><td> 9.87393</td><td> 9.78035</td><td> 9.77831</td><td> 9.77569</td><td> 9.75653</td><td>⋯</td><td> 9.60899</td><td> 9.72797</td><td>10.02570</td><td>10.04280</td><td>10.07610</td><td> 9.99215</td><td> 9.84565</td><td> 9.89990</td><td>10.13640</td><td> 9.89509</td></tr>
	<tr><th scope=row>8180411</th><td> 8.86882</td><td> 8.89400</td><td> 8.80103</td><td> 8.67483</td><td> 7.01018</td><td> 8.33998</td><td> 6.94540</td><td> 8.59880</td><td> 8.16177</td><td> 8.38775</td><td>⋯</td><td> 8.50608</td><td> 8.68022</td><td> 7.94113</td><td> 8.06512</td><td> 7.53432</td><td> 8.37084</td><td> 8.17369</td><td> 8.70992</td><td> 8.09877</td><td> 9.38929</td></tr>
	<tr><th scope=row>8180412</th><td> 7.33219</td><td> 7.55893</td><td> 6.78430</td><td> 6.84587</td><td> 7.18445</td><td> 7.23077</td><td> 7.01315</td><td> 7.22381</td><td> 7.37274</td><td> 7.53646</td><td>⋯</td><td> 7.03117</td><td> 7.09762</td><td> 7.14922</td><td> 7.19199</td><td> 7.18049</td><td> 7.16591</td><td> 7.25525</td><td> 7.45072</td><td> 7.10824</td><td> 7.12743</td></tr>
	<tr><th scope=row>8180413</th><td> 7.14770</td><td> 7.36702</td><td> 6.58836</td><td> 6.64284</td><td> 7.03858</td><td> 7.00432</td><td> 6.73689</td><td> 7.03558</td><td> 7.17611</td><td> 7.24642</td><td>⋯</td><td> 6.80760</td><td> 6.83930</td><td> 6.95357</td><td> 7.01256</td><td> 6.98837</td><td> 6.87821</td><td> 7.17545</td><td> 7.30428</td><td> 6.91156</td><td> 6.92976</td></tr>
	<tr><th scope=row>8180414</th><td> 5.83472</td><td> 5.61697</td><td> 5.27139</td><td> 5.67863</td><td> 5.24630</td><td> 5.38700</td><td> 5.61434</td><td> 5.60667</td><td> 5.24276</td><td> 5.36182</td><td>⋯</td><td> 5.48344</td><td> 5.52821</td><td> 5.48791</td><td> 5.28077</td><td> 5.54975</td><td> 5.83203</td><td> 5.34092</td><td> 5.61503</td><td> 5.60587</td><td> 5.41597</td></tr>
	<tr><th scope=row>8180415</th><td> 7.72986</td><td> 7.64719</td><td> 7.82310</td><td> 7.51615</td><td> 7.78083</td><td> 7.73307</td><td> 7.95015</td><td> 7.82931</td><td> 7.90996</td><td> 7.99317</td><td>⋯</td><td> 7.82651</td><td> 7.94766</td><td> 7.98722</td><td> 7.94047</td><td> 8.14802</td><td> 7.81312</td><td> 7.96662</td><td> 7.75461</td><td> 7.85791</td><td> 7.36329</td></tr>
	<tr><th scope=row>8180416</th><td> 5.40042</td><td> 5.10538</td><td> 4.85939</td><td> 4.95673</td><td> 5.37243</td><td> 5.18263</td><td> 5.21231</td><td> 4.99286</td><td> 4.79942</td><td> 5.13801</td><td>⋯</td><td> 4.90392</td><td> 5.24170</td><td> 4.98484</td><td> 5.05982</td><td> 5.06662</td><td> 5.13604</td><td> 4.83894</td><td> 5.10050</td><td> 5.11261</td><td> 5.21177</td></tr>
	<tr><th scope=row>8180417</th><td> 8.48381</td><td> 8.54317</td><td> 8.69660</td><td> 8.48338</td><td> 8.67200</td><td> 8.50409</td><td> 8.56539</td><td> 8.50912</td><td> 8.46662</td><td> 8.53186</td><td>⋯</td><td> 8.44144</td><td> 8.47813</td><td> 8.77304</td><td> 8.51178</td><td> 8.54516</td><td> 8.57926</td><td> 8.59788</td><td> 8.47102</td><td> 8.61478</td><td> 8.37552</td></tr>
	<tr><th scope=row>8180418</th><td> 8.61794</td><td> 8.52306</td><td> 8.15599</td><td> 7.66275</td><td> 8.29000</td><td> 8.15059</td><td> 7.86556</td><td> 8.43660</td><td> 8.25621</td><td> 7.73094</td><td>⋯</td><td> 8.28558</td><td> 8.33585</td><td> 8.44541</td><td> 8.15321</td><td> 8.49835</td><td> 8.48145</td><td> 8.52579</td><td> 8.57159</td><td> 8.29087</td><td> 8.01442</td></tr>
</tbody>
</table>




```R
ids <- AnnoProbe::idmap('GPL6244',type = 'soft')
```

    file downloaded in /home/jfckkiu/AS_HG/Final_Results/Figure2
    



```R
ids   <- ids  %>% 
           mutate(ID = as.character(ID))
```


```R
ids
```


<table class="dataframe">
<caption>A data.frame: 33297 × 2</caption>
<thead>
	<tr><th scope=col>ID</th><th scope=col>symbol</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>7896736</td><td>NA          </td></tr>
	<tr><td>7896738</td><td>OR4G2P      </td></tr>
	<tr><td>7896740</td><td>OR4F4       </td></tr>
	<tr><td>7896742</td><td>LOC728323   </td></tr>
	<tr><td>7896744</td><td>OR4F29      </td></tr>
	<tr><td>7896746</td><td>MT-TM       </td></tr>
	<tr><td>7896748</td><td>MT-TW       </td></tr>
	<tr><td>7896750</td><td>MT-TD       </td></tr>
	<tr><td>7896752</td><td>MT-TK       </td></tr>
	<tr><td>7896754</td><td>LOC100287497</td></tr>
	<tr><td>7896756</td><td>FAM87B      </td></tr>
	<tr><td>7896759</td><td>LINC01128   </td></tr>
	<tr><td>7896761</td><td>SAMD11      </td></tr>
	<tr><td>7896779</td><td>KLHL17      </td></tr>
	<tr><td>7896798</td><td>PLEKHN1     </td></tr>
	<tr><td>7896817</td><td>ISG15       </td></tr>
	<tr><td>7896822</td><td>AGRN        </td></tr>
	<tr><td>7896859</td><td>MIR200B     </td></tr>
	<tr><td>7896861</td><td>MIR200A     </td></tr>
	<tr><td>7896863</td><td>MIR429      </td></tr>
	<tr><td>7896865</td><td>TTLL10      </td></tr>
	<tr><td>7896878</td><td>B3GALT6     </td></tr>
	<tr><td>7896882</td><td>SCNN1D      </td></tr>
	<tr><td>7896908</td><td>PUSL1       </td></tr>
	<tr><td>7896917</td><td>CPTP        </td></tr>
	<tr><td>7896921</td><td>TAS1R3      </td></tr>
	<tr><td>7896929</td><td>VWA1        </td></tr>
	<tr><td>7896937</td><td>ATAD3C      </td></tr>
	<tr><td>7896952</td><td>ATAD3A      </td></tr>
	<tr><td>7896961</td><td>ATAD3B      </td></tr>
	<tr><td>⋮</td><td>⋮</td></tr>
	<tr><td>7896701</td><td>NA</td></tr>
	<tr><td>7896702</td><td>NA</td></tr>
	<tr><td>7896703</td><td>NA</td></tr>
	<tr><td>7896704</td><td>NA</td></tr>
	<tr><td>7896705</td><td>NA</td></tr>
	<tr><td>7896706</td><td>NA</td></tr>
	<tr><td>7896707</td><td>NA</td></tr>
	<tr><td>7896708</td><td>NA</td></tr>
	<tr><td>7896709</td><td>NA</td></tr>
	<tr><td>7896710</td><td>NA</td></tr>
	<tr><td>7896711</td><td>NA</td></tr>
	<tr><td>7896712</td><td>NA</td></tr>
	<tr><td>7896713</td><td>NA</td></tr>
	<tr><td>7896714</td><td>NA</td></tr>
	<tr><td>7896715</td><td>NA</td></tr>
	<tr><td>7896716</td><td>NA</td></tr>
	<tr><td>7896717</td><td>NA</td></tr>
	<tr><td>7896718</td><td>NA</td></tr>
	<tr><td>7896719</td><td>NA</td></tr>
	<tr><td>7896720</td><td>NA</td></tr>
	<tr><td>7896721</td><td>NA</td></tr>
	<tr><td>7896722</td><td>NA</td></tr>
	<tr><td>7896723</td><td>NA</td></tr>
	<tr><td>7896724</td><td>NA</td></tr>
	<tr><td>7896725</td><td>NA</td></tr>
	<tr><td>7896726</td><td>NA</td></tr>
	<tr><td>7896727</td><td>NA</td></tr>
	<tr><td>7896728</td><td>NA</td></tr>
	<tr><td>7896729</td><td>NA</td></tr>
	<tr><td>7896730</td><td>NA</td></tr>
</tbody>
</table>




```R
normalized_gset %>% 
 as.data.frame() %>% 
  rownames_to_column(var = 'Probe') %>%  
   left_join(ids, by = c('Probe'= 'ID')) %>% 
    filter(Probe %in% '7896740')
```


<table class="dataframe">
<caption>A data.frame: 1 × 66</caption>
<thead>
	<tr><th scope=col>Probe</th><th scope=col>GSM1060117</th><th scope=col>GSM1060118</th><th scope=col>GSM1060119</th><th scope=col>GSM1060120</th><th scope=col>GSM1060121</th><th scope=col>GSM1060122</th><th scope=col>GSM1060123</th><th scope=col>GSM1060124</th><th scope=col>GSM1060125</th><th scope=col>⋯</th><th scope=col>GSM1060172</th><th scope=col>GSM1060173</th><th scope=col>GSM1060174</th><th scope=col>GSM1060175</th><th scope=col>GSM1060176</th><th scope=col>GSM1060177</th><th scope=col>GSM1060178</th><th scope=col>GSM1060179</th><th scope=col>GSM1060180</th><th scope=col>symbol</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>7896740</td><td>3.51283</td><td>2.93909</td><td>3.51641</td><td>3.19208</td><td>3.29371</td><td>3.20059</td><td>3.27502</td><td>3.14589</td><td>3.24979</td><td>⋯</td><td>2.9636</td><td>3.01715</td><td>3.18444</td><td>2.91397</td><td>3.07929</td><td>3.0356</td><td>3.04976</td><td>3.17131</td><td>3.3187</td><td>OR4F4</td></tr>
</tbody>
</table>




```R
sample_name <- colnames(normalized_gset %>% 
 as.data.frame()  %>% 
  rownames_to_column(var = 'Probe')  %>%  
   left_join(ids, by = c('Probe'= 'ID')) %>% 
    filter(!is.na(symbol))  %>% 
     .[,2:(ncol(.)-1)]) 
sample_name    
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'GSM1060117'</li><li>'GSM1060118'</li><li>'GSM1060119'</li><li>'GSM1060120'</li><li>'GSM1060121'</li><li>'GSM1060122'</li><li>'GSM1060123'</li><li>'GSM1060124'</li><li>'GSM1060125'</li><li>'GSM1060126'</li><li>'GSM1060127'</li><li>'GSM1060128'</li><li>'GSM1060129'</li><li>'GSM1060130'</li><li>'GSM1060131'</li><li>'GSM1060132'</li><li>'GSM1060133'</li><li>'GSM1060134'</li><li>'GSM1060135'</li><li>'GSM1060136'</li><li>'GSM1060137'</li><li>'GSM1060138'</li><li>'GSM1060139'</li><li>'GSM1060140'</li><li>'GSM1060141'</li><li>'GSM1060142'</li><li>'GSM1060143'</li><li>'GSM1060144'</li><li>'GSM1060145'</li><li>'GSM1060146'</li><li>'GSM1060147'</li><li>'GSM1060148'</li><li>'GSM1060149'</li><li>'GSM1060150'</li><li>'GSM1060151'</li><li>'GSM1060152'</li><li>'GSM1060153'</li><li>'GSM1060154'</li><li>'GSM1060155'</li><li>'GSM1060156'</li><li>'GSM1060157'</li><li>'GSM1060158'</li><li>'GSM1060159'</li><li>'GSM1060160'</li><li>'GSM1060161'</li><li>'GSM1060162'</li><li>'GSM1060163'</li><li>'GSM1060164'</li><li>'GSM1060165'</li><li>'GSM1060166'</li><li>'GSM1060167'</li><li>'GSM1060168'</li><li>'GSM1060169'</li><li>'GSM1060170'</li><li>'GSM1060171'</li><li>'GSM1060172'</li><li>'GSM1060173'</li><li>'GSM1060174'</li><li>'GSM1060175'</li><li>'GSM1060176'</li><li>'GSM1060177'</li><li>'GSM1060178'</li><li>'GSM1060179'</li><li>'GSM1060180'</li></ol>




```R
forum_list  <- vector('list',length(sample_name))
```


```R
for (i in 1:length(sample_name)){
           forum_list[[i]]  <-  paste0(sample_name[i],'=mean(',sample_name[i],'),')
}
```


```R
forum  <- c()
for (i in 1:length(forum_list)){
   forum  <- str_c(forum, forum_list[[i]])
}
```


```R
forum
```


'GSM1060117=mean(GSM1060117),GSM1060118=mean(GSM1060118),GSM1060119=mean(GSM1060119),GSM1060120=mean(GSM1060120),GSM1060121=mean(GSM1060121),GSM1060122=mean(GSM1060122),GSM1060123=mean(GSM1060123),GSM1060124=mean(GSM1060124),GSM1060125=mean(GSM1060125),GSM1060126=mean(GSM1060126),GSM1060127=mean(GSM1060127),GSM1060128=mean(GSM1060128),GSM1060129=mean(GSM1060129),GSM1060130=mean(GSM1060130),GSM1060131=mean(GSM1060131),GSM1060132=mean(GSM1060132),GSM1060133=mean(GSM1060133),GSM1060134=mean(GSM1060134),GSM1060135=mean(GSM1060135),GSM1060136=mean(GSM1060136),GSM1060137=mean(GSM1060137),GSM1060138=mean(GSM1060138),GSM1060139=mean(GSM1060139),GSM1060140=mean(GSM1060140),GSM1060141=mean(GSM1060141),GSM1060142=mean(GSM1060142),GSM1060143=mean(GSM1060143),GSM1060144=mean(GSM1060144),GSM1060145=mean(GSM1060145),GSM1060146=mean(GSM1060146),GSM1060147=mean(GSM1060147),GSM1060148=mean(GSM1060148),GSM1060149=mean(GSM1060149),GSM1060150=mean(GSM1060150),GSM1060151=mean(GSM1060151),GSM1060152=mean(GSM1060152),GSM1060153=mean(GSM1060153),GSM1060154=mean(GSM1060154),GSM1060155=mean(GSM1060155),GSM1060156=mean(GSM1060156),GSM1060157=mean(GSM1060157),GSM1060158=mean(GSM1060158),GSM1060159=mean(GSM1060159),GSM1060160=mean(GSM1060160),GSM1060161=mean(GSM1060161),GSM1060162=mean(GSM1060162),GSM1060163=mean(GSM1060163),GSM1060164=mean(GSM1060164),GSM1060165=mean(GSM1060165),GSM1060166=mean(GSM1060166),GSM1060167=mean(GSM1060167),GSM1060168=mean(GSM1060168),GSM1060169=mean(GSM1060169),GSM1060170=mean(GSM1060170),GSM1060171=mean(GSM1060171),GSM1060172=mean(GSM1060172),GSM1060173=mean(GSM1060173),GSM1060174=mean(GSM1060174),GSM1060175=mean(GSM1060175),GSM1060176=mean(GSM1060176),GSM1060177=mean(GSM1060177),GSM1060178=mean(GSM1060178),GSM1060179=mean(GSM1060179),GSM1060180=mean(GSM1060180),'



```R
normalized_gset  %>% 
 as.data.frame()  %>% 
  rownames_to_column(var = 'Probe')  %>%  
   left_join(ids, by = c('Probe'= 'ID')) %>% 
    filter(!is.na(symbol))   %>% 
     .[,2:ncol(.)]  %>% 
       head()
```


<table class="dataframe">
<caption>A data.frame: 6 × 65</caption>
<thead>
	<tr><th></th><th scope=col>GSM1060117</th><th scope=col>GSM1060118</th><th scope=col>GSM1060119</th><th scope=col>GSM1060120</th><th scope=col>GSM1060121</th><th scope=col>GSM1060122</th><th scope=col>GSM1060123</th><th scope=col>GSM1060124</th><th scope=col>GSM1060125</th><th scope=col>GSM1060126</th><th scope=col>⋯</th><th scope=col>GSM1060172</th><th scope=col>GSM1060173</th><th scope=col>GSM1060174</th><th scope=col>GSM1060175</th><th scope=col>GSM1060176</th><th scope=col>GSM1060177</th><th scope=col>GSM1060178</th><th scope=col>GSM1060179</th><th scope=col>GSM1060180</th><th scope=col>symbol</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>2.80540</td><td>2.97946</td><td>2.99950</td><td>2.97630</td><td>2.90421</td><td>2.78402</td><td>2.88611</td><td>2.88375</td><td>3.07667</td><td>2.83233</td><td>⋯</td><td>2.85150</td><td>3.00979</td><td>2.93718</td><td>2.75821</td><td>2.94500</td><td>2.99704</td><td>3.06945</td><td>2.89322</td><td>2.81862</td><td>OR4G2P   </td></tr>
	<tr><th scope=row>2</th><td>3.51283</td><td>2.93909</td><td>3.51641</td><td>3.19208</td><td>3.29371</td><td>3.20059</td><td>3.27502</td><td>3.14589</td><td>3.24979</td><td>2.96047</td><td>⋯</td><td>2.96360</td><td>3.01715</td><td>3.18444</td><td>2.91397</td><td>3.07929</td><td>3.03560</td><td>3.04976</td><td>3.17131</td><td>3.31870</td><td>OR4F4    </td></tr>
	<tr><th scope=row>3</th><td>6.88099</td><td>7.06372</td><td>7.74086</td><td>7.80900</td><td>7.64122</td><td>7.73806</td><td>7.90782</td><td>7.40955</td><td>7.00492</td><td>7.26310</td><td>⋯</td><td>7.16201</td><td>7.50545</td><td>7.18688</td><td>7.29475</td><td>7.30564</td><td>7.66944</td><td>6.81594</td><td>6.83839</td><td>6.92016</td><td>LOC728323</td></tr>
	<tr><th scope=row>4</th><td>5.37283</td><td>5.03621</td><td>5.46020</td><td>5.84442</td><td>5.74767</td><td>5.94580</td><td>5.68829</td><td>5.94922</td><td>5.39755</td><td>5.76010</td><td>⋯</td><td>5.49711</td><td>6.04228</td><td>6.05929</td><td>5.69216</td><td>6.28010</td><td>5.59950</td><td>5.88187</td><td>5.51909</td><td>5.84261</td><td>OR4F29   </td></tr>
	<tr><th scope=row>5</th><td>8.39527</td><td>8.42700</td><td>8.04820</td><td>8.01527</td><td>8.44434</td><td>8.32162</td><td>8.49673</td><td>8.25884</td><td>8.82199</td><td>8.38388</td><td>⋯</td><td>8.24847</td><td>8.58651</td><td>8.28113</td><td>8.40412</td><td>8.91665</td><td>8.89946</td><td>9.04654</td><td>8.35727</td><td>8.32857</td><td>MT-TM    </td></tr>
	<tr><th scope=row>6</th><td>9.09767</td><td>8.91798</td><td>9.18361</td><td>9.05934</td><td>9.10794</td><td>8.85285</td><td>9.01828</td><td>9.14659</td><td>9.27742</td><td>9.26185</td><td>⋯</td><td>7.66310</td><td>8.19807</td><td>7.84545</td><td>8.12853</td><td>8.21389</td><td>8.55925</td><td>8.22249</td><td>7.74900</td><td>7.33013</td><td>MT-TW    </td></tr>
</tbody>
</table>




```R
normalized_gset_mean <- normalized_gset  %>% 
 as.data.frame()  %>% 
  rownames_to_column(var = 'Probe')  %>%  
   left_join(ids, by = c('Probe'= 'ID')) %>% 
    filter(!is.na(symbol))   %>% 
     .[,2:ncol(.)]   %>% 
       group_by(symbol)   %>% 
         summarise(GSM1060117=mean(GSM1060117),GSM1060118=mean(GSM1060118),GSM1060119=mean(GSM1060119),GSM1060120=mean(GSM1060120),GSM1060121=mean(GSM1060121),GSM1060122=mean(GSM1060122),GSM1060123=mean(GSM1060123),GSM1060124=mean(GSM1060124),GSM1060125=mean(GSM1060125),GSM1060126=mean(GSM1060126),GSM1060127=mean(GSM1060127),GSM1060128=mean(GSM1060128),GSM1060129=mean(GSM1060129),GSM1060130=mean(GSM1060130),GSM1060131=mean(GSM1060131),GSM1060132=mean(GSM1060132),GSM1060133=mean(GSM1060133),GSM1060134=mean(GSM1060134),GSM1060135=mean(GSM1060135),GSM1060136=mean(GSM1060136),GSM1060137=mean(GSM1060137),GSM1060138=mean(GSM1060138),GSM1060139=mean(GSM1060139),GSM1060140=mean(GSM1060140),GSM1060141=mean(GSM1060141),GSM1060142=mean(GSM1060142),GSM1060143=mean(GSM1060143),GSM1060144=mean(GSM1060144),GSM1060145=mean(GSM1060145),GSM1060146=mean(GSM1060146),GSM1060147=mean(GSM1060147),GSM1060148=mean(GSM1060148),GSM1060149=mean(GSM1060149),GSM1060150=mean(GSM1060150),GSM1060151=mean(GSM1060151),GSM1060152=mean(GSM1060152),GSM1060153=mean(GSM1060153),GSM1060154=mean(GSM1060154),GSM1060155=mean(GSM1060155),GSM1060156=mean(GSM1060156),GSM1060157=mean(GSM1060157),GSM1060158=mean(GSM1060158),GSM1060159=mean(GSM1060159),GSM1060160=mean(GSM1060160),GSM1060161=mean(GSM1060161),GSM1060162=mean(GSM1060162),GSM1060163=mean(GSM1060163),GSM1060164=mean(GSM1060164),GSM1060165=mean(GSM1060165),GSM1060166=mean(GSM1060166),GSM1060167=mean(GSM1060167),GSM1060168=mean(GSM1060168),GSM1060169=mean(GSM1060169),GSM1060170=mean(GSM1060170),GSM1060171=mean(GSM1060171),GSM1060172=mean(GSM1060172),GSM1060173=mean(GSM1060173),GSM1060174=mean(GSM1060174),GSM1060175=mean(GSM1060175),GSM1060176=mean(GSM1060176),GSM1060177=mean(GSM1060177),GSM1060178=mean(GSM1060178),GSM1060179=mean(GSM1060179),GSM1060180=mean(GSM1060180))%>% 
         column_to_rownames(var ='symbol')
```


```R
gset$GSE43292_series_matrix.txt.gz@phenoData@data  %>%  head()
```


<table class="dataframe">
<caption>A data.frame: 6 × 37</caption>
<thead>
	<tr><th></th><th scope=col>title</th><th scope=col>geo_accession</th><th scope=col>status</th><th scope=col>submission_date</th><th scope=col>last_update_date</th><th scope=col>type</th><th scope=col>channel_count</th><th scope=col>source_name_ch1</th><th scope=col>organism_ch1</th><th scope=col>characteristics_ch1</th><th scope=col>⋯</th><th scope=col>contact_address</th><th scope=col>contact_city</th><th scope=col>contact_zip/postal_code</th><th scope=col>contact_country</th><th scope=col>supplementary_file</th><th scope=col>supplementary_file.1</th><th scope=col>data_row_count</th><th scope=col>disease state:ch1</th><th scope=col>patient:ch1</th><th scope=col>tissue:ch1</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM1060117</th><td>Patient 1, Macroscopically intact tissue</td><td>GSM1060117</td><td>Public on Apr 01 2013</td><td>Jan 04 2013</td><td>Apr 01 2013</td><td>RNA</td><td>1</td><td>Macroscopically intact carotid tissue adjacent to the atheroma plaque</td><td>Homo sapiens</td><td>patient: 1</td><td>⋯</td><td>8 avenue Rockefeller</td><td>Lyon</td><td>69373</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060117/suppl/GSM1060117_JLI_1_HuGene.CEL.gz</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060117/suppl/GSM1060117_JLI_1_HuGene.rma-gene-default.chp.gz</td><td>33297</td><td>hypertension</td><td>1</td><td>Macroscopically intact tissue</td></tr>
	<tr><th scope=row>GSM1060118</th><td>Patient 1, Atheroma plaque              </td><td>GSM1060118</td><td>Public on Apr 01 2013</td><td>Jan 04 2013</td><td>Apr 01 2013</td><td>RNA</td><td>1</td><td>Carotid atheroma plaque obtained from endarterectomy                 </td><td>Homo sapiens</td><td>patient: 1</td><td>⋯</td><td>8 avenue Rockefeller</td><td>Lyon</td><td>69373</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060118/suppl/GSM1060118_JLI_2_HuGene.CEL.gz</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060118/suppl/GSM1060118_JLI_2_HuGene.rma-gene-default.chp.gz</td><td>33297</td><td>hypertension</td><td>1</td><td>Atheroma plaque              </td></tr>
	<tr><th scope=row>GSM1060119</th><td>Patient 2, Macroscopically intact tissue</td><td>GSM1060119</td><td>Public on Apr 01 2013</td><td>Jan 04 2013</td><td>Apr 01 2013</td><td>RNA</td><td>1</td><td>Macroscopically intact carotid tissue adjacent to the atheroma plaque</td><td>Homo sapiens</td><td>patient: 2</td><td>⋯</td><td>8 avenue Rockefeller</td><td>Lyon</td><td>69373</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060119/suppl/GSM1060119_JLI_3_HuGene.CEL.gz</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060119/suppl/GSM1060119_JLI_3_HuGene.rma-gene-default.chp.gz</td><td>33297</td><td>hypertension</td><td>2</td><td>Macroscopically intact tissue</td></tr>
	<tr><th scope=row>GSM1060120</th><td>Patient 2, Atheroma plaque              </td><td>GSM1060120</td><td>Public on Apr 01 2013</td><td>Jan 04 2013</td><td>Apr 01 2013</td><td>RNA</td><td>1</td><td>Carotid atheroma plaque obtained from endarterectomy                 </td><td>Homo sapiens</td><td>patient: 2</td><td>⋯</td><td>8 avenue Rockefeller</td><td>Lyon</td><td>69373</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060120/suppl/GSM1060120_JLI_4_HuGene.CEL.gz</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060120/suppl/GSM1060120_JLI_4_HuGene.rma-gene-default.chp.gz</td><td>33297</td><td>hypertension</td><td>2</td><td>Atheroma plaque              </td></tr>
	<tr><th scope=row>GSM1060121</th><td>Patient 3, Macroscopically intact tissue</td><td>GSM1060121</td><td>Public on Apr 01 2013</td><td>Jan 04 2013</td><td>Apr 01 2013</td><td>RNA</td><td>1</td><td>Macroscopically intact carotid tissue adjacent to the atheroma plaque</td><td>Homo sapiens</td><td>patient: 3</td><td>⋯</td><td>8 avenue Rockefeller</td><td>Lyon</td><td>69373</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060121/suppl/GSM1060121_JLI_7_HuGene.CEL.gz</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060121/suppl/GSM1060121_JLI_7_HuGene.rma-gene-default.chp.gz</td><td>33297</td><td>hypertension</td><td>3</td><td>Macroscopically intact tissue</td></tr>
	<tr><th scope=row>GSM1060122</th><td>Patient 3, Atheroma plaque              </td><td>GSM1060122</td><td>Public on Apr 01 2013</td><td>Jan 04 2013</td><td>Apr 01 2013</td><td>RNA</td><td>1</td><td>Carotid atheroma plaque obtained from endarterectomy                 </td><td>Homo sapiens</td><td>patient: 3</td><td>⋯</td><td>8 avenue Rockefeller</td><td>Lyon</td><td>69373</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060122/suppl/GSM1060122_JLI_8_HuGene.CEL.gz</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1060nnn/GSM1060122/suppl/GSM1060122_JLI_8_HuGene.rma-gene-default.chp.gz</td><td>33297</td><td>hypertension</td><td>3</td><td>Atheroma plaque              </td></tr>
</tbody>
</table>




```R
group_list  <- factor(gset$GSE43292_series_matrix.txt.gz@phenoData@data[,c('source_name_ch1')]) 
```


```R
levels(group_list)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Carotid atheroma plaque obtained from endarterectomy'</li><li>'Macroscopically intact carotid tissue adjacent to the atheroma plaque'</li></ol>




```R
group_list  <- fct_recode(group_list, 'MIT' = 'Macroscopically intact carotid tissue adjacent to the atheroma plaque',
                                      'ATH' = 'Carotid atheroma plaque obtained from endarterectomy')
```


```R
group_list <- relevel(group_list, ref="MIT")
```


```R
group_list
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li><li>MIT</li><li>ATH</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'MIT'</li><li>'ATH'</li></ol>
</details>



```R
design <- model.matrix(~0 + group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(normalized_gset_mean)
```


```R
design
```


<table class="dataframe">
<caption>A matrix: 64 × 2 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>MIT</th><th scope=col>ATH</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM1060117</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060118</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060119</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060120</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060121</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060122</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060123</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060124</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060125</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060126</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060127</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060128</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060129</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060130</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060131</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060132</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060133</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060134</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060135</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060136</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060137</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060138</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060139</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060140</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060141</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060142</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060143</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060144</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060145</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060146</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>GSM1060151</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060152</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060153</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060154</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060155</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060156</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060157</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060158</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060159</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060160</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060161</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060162</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060163</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060164</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060165</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060166</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060167</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060168</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060169</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060170</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060171</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060172</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060173</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060174</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060175</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060176</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060177</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060178</th><td>0</td><td>1</td></tr>
	<tr><th scope=row>GSM1060179</th><td>1</td><td>0</td></tr>
	<tr><th scope=row>GSM1060180</th><td>0</td><td>1</td></tr>
</tbody>
</table>




```R
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix
```


<table class="dataframe">
<caption>A matrix: 2 × 1 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>MIT-ATH</th></tr>
</thead>
<tbody>
	<tr><th scope=row>MIT</th><td> 1</td></tr>
	<tr><th scope=row>ATH</th><td>-1</td></tr>
</tbody>
</table>




```R
##step1
fit <- lmFit(normalized_gset_mean,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
nrDEG <- nrDEG %>% 
mutate(logFC = -logFC) %>% 
arrange(desc(logFC)) %>%
rownames_to_column(var = 'Gene_Symbol') 
##export
```


```R
normalized_gset_mean
```


<table class="dataframe">
<caption>A data.frame: 23307 × 64</caption>
<thead>
	<tr><th></th><th scope=col>GSM1060117</th><th scope=col>GSM1060118</th><th scope=col>GSM1060119</th><th scope=col>GSM1060120</th><th scope=col>GSM1060121</th><th scope=col>GSM1060122</th><th scope=col>GSM1060123</th><th scope=col>GSM1060124</th><th scope=col>GSM1060125</th><th scope=col>GSM1060126</th><th scope=col>⋯</th><th scope=col>GSM1060171</th><th scope=col>GSM1060172</th><th scope=col>GSM1060173</th><th scope=col>GSM1060174</th><th scope=col>GSM1060175</th><th scope=col>GSM1060176</th><th scope=col>GSM1060177</th><th scope=col>GSM1060178</th><th scope=col>GSM1060179</th><th scope=col>GSM1060180</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>A1BG</th><td> 6.33413</td><td> 6.23136</td><td> 6.17699</td><td> 6.28800</td><td> 6.205660</td><td> 5.995310</td><td> 6.162800</td><td> 5.914780</td><td> 6.172340</td><td> 6.187750</td><td>⋯</td><td> 5.884050</td><td> 6.097230</td><td> 5.809070</td><td> 5.962630</td><td> 5.929280</td><td> 6.054380</td><td> 5.831460</td><td> 6.352940</td><td> 6.144100</td><td> 6.11485</td></tr>
	<tr><th scope=row>A1CF</th><td> 4.10072</td><td> 4.01097</td><td> 3.96745</td><td> 4.02464</td><td> 4.006060</td><td> 4.014880</td><td> 4.276880</td><td> 3.950280</td><td> 4.148390</td><td> 3.829580</td><td>⋯</td><td> 3.849670</td><td> 3.874750</td><td> 3.768630</td><td> 3.833250</td><td> 4.138480</td><td> 3.845660</td><td> 3.928350</td><td> 3.939520</td><td> 4.118220</td><td> 3.96717</td></tr>
	<tr><th scope=row>A2M</th><td>12.49600</td><td>12.26350</td><td>12.92450</td><td>12.97690</td><td>12.933100</td><td>12.958000</td><td>12.692200</td><td>12.813700</td><td>12.734900</td><td>12.929600</td><td>⋯</td><td>12.636900</td><td>12.602800</td><td>12.873400</td><td>12.867000</td><td>12.930700</td><td>12.865700</td><td>12.335300</td><td>12.131500</td><td>12.977600</td><td>12.78230</td></tr>
	<tr><th scope=row>A2ML1</th><td> 4.03661</td><td> 4.10083</td><td> 3.81421</td><td> 4.16249</td><td> 4.030510</td><td> 4.130190</td><td> 4.105710</td><td> 3.739370</td><td> 3.740730</td><td> 3.956200</td><td>⋯</td><td> 4.140920</td><td> 3.956720</td><td> 3.950060</td><td> 4.059610</td><td> 4.022620</td><td> 3.899440</td><td> 4.085170</td><td> 4.255340</td><td> 4.189320</td><td> 3.95617</td></tr>
	<tr><th scope=row>A3GALT2</th><td> 4.83119</td><td> 4.82828</td><td> 4.95262</td><td> 4.93172</td><td> 4.727280</td><td> 4.836350</td><td> 4.982820</td><td> 4.730410</td><td> 4.607110</td><td> 4.670480</td><td>⋯</td><td> 4.481430</td><td> 4.880340</td><td> 4.651580</td><td> 4.920190</td><td> 5.253800</td><td> 5.302240</td><td> 4.900220</td><td> 5.061330</td><td> 5.091100</td><td> 4.88142</td></tr>
	<tr><th scope=row>A4GALT</th><td> 7.33187</td><td> 7.26400</td><td> 6.81011</td><td> 6.81327</td><td> 6.873140</td><td> 6.717460</td><td> 6.654660</td><td> 7.011190</td><td> 7.038280</td><td> 6.839140</td><td>⋯</td><td> 6.578780</td><td> 6.307620</td><td> 6.593280</td><td> 6.963220</td><td> 6.698600</td><td> 6.508300</td><td> 7.005950</td><td> 6.676420</td><td> 6.846110</td><td> 6.55254</td></tr>
	<tr><th scope=row>A4GNT</th><td> 3.57305</td><td> 3.74132</td><td> 3.44454</td><td> 3.56300</td><td> 3.536210</td><td> 3.670310</td><td> 3.520020</td><td> 3.640490</td><td> 3.651710</td><td> 3.493550</td><td>⋯</td><td> 3.366450</td><td> 3.507460</td><td> 3.458870</td><td> 3.171950</td><td> 3.447070</td><td> 3.573210</td><td> 3.406360</td><td> 3.744170</td><td> 3.683410</td><td> 3.58425</td></tr>
	<tr><th scope=row>AAAS</th><td> 7.25019</td><td> 7.73544</td><td> 7.84654</td><td> 7.83944</td><td> 8.057220</td><td> 7.962680</td><td> 7.852890</td><td> 8.029380</td><td> 7.854730</td><td> 8.011330</td><td>⋯</td><td> 7.899490</td><td> 7.900170</td><td> 7.866410</td><td> 7.802050</td><td> 7.939140</td><td> 7.964970</td><td> 7.880880</td><td> 7.863170</td><td> 8.069680</td><td> 7.78399</td></tr>
	<tr><th scope=row>AACS</th><td> 6.52156</td><td> 6.56949</td><td> 6.50518</td><td> 6.16380</td><td> 6.794600</td><td> 6.630790</td><td> 6.713580</td><td> 6.689270</td><td> 6.841120</td><td> 6.664670</td><td>⋯</td><td> 6.558760</td><td> 6.312700</td><td> 6.389090</td><td> 6.908250</td><td> 6.357360</td><td> 6.694810</td><td> 6.294770</td><td> 6.590790</td><td> 6.427850</td><td> 6.32270</td></tr>
	<tr><th scope=row>AADAC</th><td> 2.84985</td><td> 2.82268</td><td> 2.79593</td><td> 2.78707</td><td> 2.768170</td><td> 2.703740</td><td> 3.043500</td><td> 3.035650</td><td> 2.780560</td><td> 2.614740</td><td>⋯</td><td> 2.779440</td><td> 2.899550</td><td> 2.668780</td><td> 2.890830</td><td> 2.848480</td><td> 2.701900</td><td> 2.746710</td><td> 2.709050</td><td> 2.676910</td><td> 2.95428</td></tr>
	<tr><th scope=row>AADACL2</th><td> 3.30588</td><td> 3.25592</td><td> 3.30354</td><td> 3.27986</td><td> 3.350960</td><td> 3.078130</td><td> 3.152720</td><td> 3.158270</td><td> 3.101850</td><td> 3.129590</td><td>⋯</td><td> 3.142570</td><td> 3.212370</td><td> 3.140080</td><td> 3.102040</td><td> 2.942280</td><td> 3.071950</td><td> 3.365670</td><td> 3.292630</td><td> 3.110640</td><td> 3.19635</td></tr>
	<tr><th scope=row>AADACL3</th><td> 3.90194</td><td> 4.28039</td><td> 4.01338</td><td> 4.37966</td><td> 4.308860</td><td> 4.302990</td><td> 4.719950</td><td> 4.420110</td><td> 4.291320</td><td> 4.082200</td><td>⋯</td><td> 3.989950</td><td> 4.300350</td><td> 4.447350</td><td> 4.065000</td><td> 4.244920</td><td> 4.321740</td><td> 4.029460</td><td> 4.611170</td><td> 4.362510</td><td> 4.16190</td></tr>
	<tr><th scope=row>AADACL4</th><td> 4.05240</td><td> 3.87569</td><td> 3.76920</td><td> 3.91195</td><td> 3.835640</td><td> 4.040580</td><td> 4.054060</td><td> 3.703250</td><td> 3.915150</td><td> 3.958500</td><td>⋯</td><td> 3.938130</td><td> 4.060640</td><td> 3.776560</td><td> 4.054650</td><td> 3.997450</td><td> 3.833180</td><td> 3.870620</td><td> 3.994450</td><td> 4.283060</td><td> 3.98729</td></tr>
	<tr><th scope=row>AADAT</th><td> 5.34443</td><td> 5.18716</td><td> 6.42704</td><td> 6.36462</td><td> 6.586680</td><td> 5.647800</td><td> 7.056810</td><td> 5.847700</td><td> 5.805940</td><td> 6.033930</td><td>⋯</td><td> 5.761140</td><td> 5.258240</td><td> 6.527180</td><td> 6.496410</td><td> 6.367350</td><td> 5.929410</td><td> 6.246180</td><td> 5.152180</td><td> 5.931100</td><td> 5.23681</td></tr>
	<tr><th scope=row>AAED1</th><td> 7.77722</td><td> 7.87651</td><td> 8.05750</td><td> 7.77185</td><td> 7.627100</td><td> 7.727870</td><td> 7.479700</td><td> 7.742790</td><td> 7.911490</td><td> 7.935770</td><td>⋯</td><td> 7.826170</td><td> 7.691350</td><td> 7.861950</td><td> 7.514610</td><td> 7.477840</td><td> 7.398790</td><td> 7.851200</td><td> 7.622080</td><td> 7.656300</td><td> 7.34300</td></tr>
	<tr><th scope=row>AAGAB</th><td> 7.85477</td><td> 7.70956</td><td> 7.19151</td><td> 7.22967</td><td> 7.488280</td><td> 7.314570</td><td> 7.575450</td><td> 7.494540</td><td> 7.642140</td><td> 7.432660</td><td>⋯</td><td> 7.758110</td><td> 7.923780</td><td> 7.718050</td><td> 7.748940</td><td> 7.570000</td><td> 7.654970</td><td> 8.127740</td><td> 7.788350</td><td> 7.501440</td><td> 7.84698</td></tr>
	<tr><th scope=row>AAK1</th><td> 8.22717</td><td> 8.31815</td><td> 8.74060</td><td> 8.66952</td><td> 8.095265</td><td> 8.313355</td><td> 8.128905</td><td> 8.447590</td><td> 8.398370</td><td> 8.440045</td><td>⋯</td><td> 8.375805</td><td> 8.456655</td><td> 8.012455</td><td> 8.164350</td><td> 7.946180</td><td> 8.132335</td><td> 7.754315</td><td> 8.327455</td><td> 7.955945</td><td> 8.44443</td></tr>
	<tr><th scope=row>AAMDC</th><td> 5.65657</td><td> 5.57825</td><td> 6.05968</td><td> 6.35897</td><td> 6.521880</td><td> 6.125580</td><td> 6.337650</td><td> 6.076220</td><td> 6.213170</td><td> 6.116430</td><td>⋯</td><td> 6.518710</td><td> 6.249870</td><td> 6.586000</td><td> 6.512410</td><td> 6.691830</td><td> 6.426220</td><td> 6.431450</td><td> 5.910120</td><td> 6.358170</td><td> 6.15294</td></tr>
	<tr><th scope=row>AAMP</th><td> 8.18316</td><td> 8.33303</td><td> 7.66461</td><td> 7.71843</td><td> 8.297810</td><td> 7.988820</td><td> 8.147270</td><td> 8.171650</td><td> 8.147440</td><td> 8.082130</td><td>⋯</td><td> 7.850910</td><td> 7.861290</td><td> 8.165260</td><td> 7.989590</td><td> 7.910050</td><td> 7.907020</td><td> 8.024960</td><td> 7.811020</td><td> 8.089790</td><td> 7.70452</td></tr>
	<tr><th scope=row>AANAT</th><td> 5.63924</td><td> 5.65591</td><td> 5.41210</td><td> 5.43153</td><td> 5.489600</td><td> 5.380810</td><td> 5.680120</td><td> 5.517040</td><td> 5.584820</td><td> 5.356170</td><td>⋯</td><td> 5.325340</td><td> 5.715690</td><td> 5.265600</td><td> 5.503580</td><td> 5.343190</td><td> 5.357820</td><td> 5.126700</td><td> 5.644050</td><td> 5.284600</td><td> 5.54677</td></tr>
	<tr><th scope=row>AAR2</th><td> 7.46769</td><td> 7.46391</td><td> 7.18075</td><td> 7.25863</td><td> 7.629600</td><td> 7.250830</td><td> 7.409180</td><td> 7.372140</td><td> 7.532570</td><td> 7.334860</td><td>⋯</td><td> 7.240740</td><td> 7.304730</td><td> 7.674550</td><td> 7.463160</td><td> 7.572510</td><td> 7.367340</td><td> 7.417170</td><td> 7.451330</td><td> 7.471430</td><td> 7.30658</td></tr>
	<tr><th scope=row>AARD</th><td> 5.08894</td><td> 5.17774</td><td> 5.13358</td><td> 5.13923</td><td> 5.215750</td><td> 4.884590</td><td> 5.094010</td><td> 4.856330</td><td> 4.774160</td><td> 4.850970</td><td>⋯</td><td> 4.852540</td><td> 4.850680</td><td> 4.903870</td><td> 5.351960</td><td> 5.208560</td><td> 4.954350</td><td> 5.153520</td><td> 4.940130</td><td> 5.100770</td><td> 5.20342</td></tr>
	<tr><th scope=row>AARS</th><td> 9.19336</td><td> 9.39812</td><td> 8.85656</td><td> 8.95554</td><td> 9.076860</td><td> 9.148440</td><td> 9.099850</td><td> 9.183370</td><td> 9.190180</td><td> 9.296830</td><td>⋯</td><td> 9.057630</td><td> 9.122210</td><td> 9.228770</td><td> 9.127030</td><td> 9.001740</td><td> 8.955210</td><td> 8.941960</td><td> 9.033640</td><td> 8.994620</td><td> 8.97487</td></tr>
	<tr><th scope=row>AARS2</th><td> 6.23486</td><td> 6.39646</td><td> 5.88791</td><td> 5.67745</td><td> 6.065220</td><td> 6.385780</td><td> 6.105200</td><td> 6.082310</td><td> 6.120040</td><td> 6.041290</td><td>⋯</td><td> 5.824280</td><td> 5.862870</td><td> 5.701840</td><td> 5.591840</td><td> 5.890080</td><td> 5.823160</td><td> 5.595190</td><td> 6.080590</td><td> 6.213680</td><td> 6.11416</td></tr>
	<tr><th scope=row>AASDH</th><td> 6.49165</td><td> 6.90503</td><td> 7.15174</td><td> 7.24475</td><td> 6.874940</td><td> 7.125700</td><td> 6.985320</td><td> 7.014320</td><td> 6.732850</td><td> 7.013900</td><td>⋯</td><td> 7.198320</td><td> 7.197440</td><td> 7.666840</td><td> 7.589570</td><td> 7.365940</td><td> 7.256180</td><td> 7.423760</td><td> 6.511200</td><td> 7.169750</td><td> 6.91563</td></tr>
	<tr><th scope=row>AASDHPPT</th><td> 5.49130</td><td> 5.50187</td><td> 5.44725</td><td> 5.23175</td><td> 5.489490</td><td> 5.362130</td><td> 5.431920</td><td> 5.311595</td><td> 5.343005</td><td> 5.442175</td><td>⋯</td><td> 5.463835</td><td> 5.524950</td><td> 5.574115</td><td> 5.531265</td><td> 5.399075</td><td> 5.289875</td><td> 5.614125</td><td> 5.227605</td><td> 5.403960</td><td> 5.43015</td></tr>
	<tr><th scope=row>AASS</th><td> 8.08306</td><td> 7.61825</td><td> 9.08592</td><td> 9.14799</td><td> 8.994880</td><td> 8.893730</td><td> 9.194020</td><td> 8.624360</td><td> 8.232910</td><td> 8.663470</td><td>⋯</td><td> 8.520820</td><td> 7.958940</td><td> 9.166090</td><td> 9.234190</td><td> 8.991620</td><td> 8.765830</td><td> 8.637630</td><td> 7.642380</td><td> 9.034170</td><td> 8.59184</td></tr>
	<tr><th scope=row>AATF</th><td> 8.20921</td><td> 8.33213</td><td> 8.10343</td><td> 7.93821</td><td> 8.259060</td><td> 8.302820</td><td> 8.241540</td><td> 8.213570</td><td> 8.316550</td><td> 8.342720</td><td>⋯</td><td> 8.073900</td><td> 8.125390</td><td> 8.173360</td><td> 8.115940</td><td> 8.202090</td><td> 8.036140</td><td> 8.214020</td><td> 8.237470</td><td> 8.103410</td><td> 8.09945</td></tr>
	<tr><th scope=row>AATK</th><td> 6.11253</td><td> 6.35266</td><td> 5.79944</td><td> 5.87782</td><td> 6.227190</td><td> 6.051150</td><td> 6.081640</td><td> 5.927340</td><td> 6.070090</td><td> 5.818400</td><td>⋯</td><td> 5.761910</td><td> 5.922940</td><td> 5.669620</td><td> 5.897280</td><td> 5.988400</td><td> 5.985300</td><td> 5.703370</td><td> 6.300550</td><td> 6.284870</td><td> 6.26960</td></tr>
	<tr><th scope=row>AATK-AS1</th><td> 6.35996</td><td> 6.45766</td><td> 6.13680</td><td> 5.98517</td><td> 6.162350</td><td> 6.092930</td><td> 6.216260</td><td> 6.006140</td><td> 6.158240</td><td> 6.033200</td><td>⋯</td><td> 5.921370</td><td> 6.126150</td><td> 5.867100</td><td> 6.149820</td><td> 6.224110</td><td> 6.028670</td><td> 6.050850</td><td> 6.381870</td><td> 6.358470</td><td> 6.29537</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>ZSCAN23</th><td>4.739290</td><td>4.44092</td><td>5.322820</td><td>5.360020</td><td>4.829890</td><td>4.831710</td><td>4.994780</td><td>4.952940</td><td>4.926850</td><td>5.052570</td><td>⋯</td><td>4.519370</td><td>4.361260</td><td>4.795260</td><td>4.951550</td><td>5.256740</td><td>4.99614</td><td>4.488520</td><td>4.642920</td><td>4.940440</td><td>4.441750</td></tr>
	<tr><th scope=row>ZSCAN25</th><td>7.255930</td><td>6.98109</td><td>7.557610</td><td>7.452420</td><td>7.342360</td><td>7.340570</td><td>7.306950</td><td>7.422280</td><td>7.279090</td><td>7.353910</td><td>⋯</td><td>6.920140</td><td>6.724410</td><td>6.944710</td><td>7.144050</td><td>7.202520</td><td>7.12622</td><td>6.845080</td><td>7.080030</td><td>6.941390</td><td>7.280240</td></tr>
	<tr><th scope=row>ZSCAN26</th><td>6.541880</td><td>6.70649</td><td>7.292780</td><td>7.093080</td><td>7.202300</td><td>7.112510</td><td>7.353800</td><td>6.982480</td><td>7.042130</td><td>7.019270</td><td>⋯</td><td>7.211700</td><td>7.222990</td><td>7.532700</td><td>7.452160</td><td>7.654810</td><td>7.23261</td><td>7.065450</td><td>6.937430</td><td>7.157190</td><td>7.225920</td></tr>
	<tr><th scope=row>ZSCAN29</th><td>7.341850</td><td>7.30534</td><td>7.587300</td><td>7.445370</td><td>7.597690</td><td>7.352240</td><td>7.439310</td><td>7.267930</td><td>7.208070</td><td>7.319900</td><td>⋯</td><td>7.166810</td><td>7.235580</td><td>7.549300</td><td>7.264300</td><td>7.268140</td><td>7.11719</td><td>7.236160</td><td>7.386040</td><td>7.251770</td><td>7.457940</td></tr>
	<tr><th scope=row>ZSCAN30</th><td>5.442550</td><td>5.32821</td><td>5.629910</td><td>5.853860</td><td>5.608910</td><td>5.523470</td><td>5.584390</td><td>5.652740</td><td>5.576710</td><td>5.622540</td><td>⋯</td><td>5.460460</td><td>5.318840</td><td>5.777450</td><td>5.556350</td><td>5.556090</td><td>5.72232</td><td>5.319860</td><td>5.590800</td><td>5.708400</td><td>5.304390</td></tr>
	<tr><th scope=row>ZSCAN31</th><td>5.055560</td><td>4.93817</td><td>5.596350</td><td>5.517290</td><td>5.702040</td><td>5.313130</td><td>5.919250</td><td>5.272140</td><td>5.461130</td><td>6.033240</td><td>⋯</td><td>5.273260</td><td>5.061410</td><td>5.989300</td><td>5.607330</td><td>5.964670</td><td>5.82656</td><td>5.777770</td><td>5.628260</td><td>5.478270</td><td>5.147750</td></tr>
	<tr><th scope=row>ZSCAN32</th><td>6.275010</td><td>6.41561</td><td>6.262170</td><td>6.061880</td><td>6.232320</td><td>6.119590</td><td>6.337180</td><td>6.236720</td><td>6.304530</td><td>6.374180</td><td>⋯</td><td>5.945780</td><td>6.179630</td><td>6.296460</td><td>6.145740</td><td>6.235130</td><td>6.05691</td><td>6.229960</td><td>6.093360</td><td>6.300080</td><td>6.140160</td></tr>
	<tr><th scope=row>ZSCAN4</th><td>3.075580</td><td>3.03255</td><td>2.733240</td><td>2.902610</td><td>2.848370</td><td>2.658400</td><td>2.869320</td><td>2.916180</td><td>3.083690</td><td>2.915580</td><td>⋯</td><td>2.856440</td><td>3.177350</td><td>2.756660</td><td>3.108830</td><td>3.005950</td><td>2.87814</td><td>2.897430</td><td>2.999380</td><td>3.024120</td><td>3.216140</td></tr>
	<tr><th scope=row>ZSCAN5A</th><td>6.010695</td><td>6.25011</td><td>5.819710</td><td>5.772685</td><td>5.862080</td><td>5.970995</td><td>5.817635</td><td>5.933545</td><td>6.100115</td><td>5.942610</td><td>⋯</td><td>5.716315</td><td>6.054235</td><td>5.790415</td><td>5.823925</td><td>5.736070</td><td>5.78424</td><td>5.702880</td><td>5.950375</td><td>5.993690</td><td>5.935370</td></tr>
	<tr><th scope=row>ZSCAN5B</th><td>3.847550</td><td>4.58524</td><td>3.580100</td><td>3.909380</td><td>3.924590</td><td>3.874920</td><td>4.009860</td><td>3.881090</td><td>3.686350</td><td>3.907280</td><td>⋯</td><td>4.030190</td><td>4.051830</td><td>3.692560</td><td>3.575180</td><td>3.669300</td><td>3.59561</td><td>3.950280</td><td>3.922240</td><td>3.656240</td><td>3.721820</td></tr>
	<tr><th scope=row>ZSCAN9</th><td>6.046260</td><td>5.96850</td><td>6.302460</td><td>6.477780</td><td>6.464210</td><td>6.432070</td><td>6.477620</td><td>6.403120</td><td>6.227690</td><td>6.463900</td><td>⋯</td><td>6.660720</td><td>6.572350</td><td>6.975170</td><td>6.543870</td><td>6.579070</td><td>6.42234</td><td>6.334720</td><td>6.553440</td><td>6.600680</td><td>6.593470</td></tr>
	<tr><th scope=row>ZSWIM1</th><td>6.114610</td><td>6.21384</td><td>6.390560</td><td>6.134040</td><td>6.421410</td><td>6.146060</td><td>6.295100</td><td>6.366980</td><td>6.387090</td><td>6.506010</td><td>⋯</td><td>6.578690</td><td>6.463700</td><td>6.612010</td><td>6.526830</td><td>6.671320</td><td>6.77247</td><td>6.365320</td><td>6.407550</td><td>6.626760</td><td>6.563900</td></tr>
	<tr><th scope=row>ZSWIM2</th><td>2.699320</td><td>2.79870</td><td>2.661550</td><td>2.695030</td><td>2.798440</td><td>2.755660</td><td>2.598900</td><td>2.535180</td><td>2.747370</td><td>2.645670</td><td>⋯</td><td>2.691240</td><td>2.592970</td><td>2.495090</td><td>2.660870</td><td>2.677300</td><td>2.67362</td><td>2.700920</td><td>2.817000</td><td>2.717310</td><td>2.668020</td></tr>
	<tr><th scope=row>ZSWIM3</th><td>4.804540</td><td>4.82555</td><td>4.729410</td><td>4.498370</td><td>5.003850</td><td>4.844300</td><td>4.690370</td><td>4.826390</td><td>4.814820</td><td>4.969180</td><td>⋯</td><td>4.779870</td><td>4.757850</td><td>4.757170</td><td>4.862460</td><td>4.706120</td><td>4.85834</td><td>4.926400</td><td>4.985530</td><td>4.819000</td><td>4.754920</td></tr>
	<tr><th scope=row>ZSWIM4</th><td>6.274790</td><td>6.67088</td><td>6.135360</td><td>6.159800</td><td>6.149640</td><td>6.312880</td><td>6.313830</td><td>6.336560</td><td>6.431320</td><td>6.282780</td><td>⋯</td><td>5.795520</td><td>6.276530</td><td>5.821120</td><td>5.958720</td><td>5.943340</td><td>5.88275</td><td>5.930330</td><td>6.672070</td><td>6.145630</td><td>6.184300</td></tr>
	<tr><th scope=row>ZSWIM5</th><td>5.813450</td><td>5.42529</td><td>6.104120</td><td>6.035930</td><td>6.081840</td><td>5.851010</td><td>6.348580</td><td>6.054770</td><td>6.160560</td><td>5.933700</td><td>⋯</td><td>5.879680</td><td>5.915860</td><td>6.120320</td><td>6.211970</td><td>6.348330</td><td>6.14024</td><td>5.826440</td><td>5.954130</td><td>6.021740</td><td>5.970800</td></tr>
	<tr><th scope=row>ZSWIM6</th><td>9.017940</td><td>9.25979</td><td>9.153690</td><td>9.113050</td><td>8.846030</td><td>9.013890</td><td>8.494060</td><td>9.051230</td><td>8.789680</td><td>9.071290</td><td>⋯</td><td>9.152630</td><td>9.486710</td><td>8.881170</td><td>8.956680</td><td>8.675820</td><td>8.86422</td><td>8.938520</td><td>9.495070</td><td>8.936340</td><td>9.556640</td></tr>
	<tr><th scope=row>ZSWIM8</th><td>6.630595</td><td>6.48948</td><td>6.233210</td><td>6.300660</td><td>6.756585</td><td>6.452325</td><td>6.538550</td><td>6.518880</td><td>6.711720</td><td>6.640125</td><td>⋯</td><td>6.336675</td><td>6.353280</td><td>6.231415</td><td>6.348865</td><td>6.464705</td><td>6.46188</td><td>6.332295</td><td>6.551470</td><td>6.507205</td><td>6.569705</td></tr>
	<tr><th scope=row>ZUFSP</th><td>5.507730</td><td>5.65530</td><td>5.917350</td><td>5.797510</td><td>5.712160</td><td>5.607790</td><td>5.680950</td><td>5.699750</td><td>5.444790</td><td>5.672050</td><td>⋯</td><td>6.261170</td><td>6.048710</td><td>5.804440</td><td>6.011210</td><td>5.742370</td><td>5.57482</td><td>6.453930</td><td>5.775680</td><td>5.550680</td><td>6.104750</td></tr>
	<tr><th scope=row>ZW10</th><td>7.345310</td><td>7.47667</td><td>7.603950</td><td>7.599710</td><td>7.459620</td><td>7.410450</td><td>7.327460</td><td>7.317760</td><td>7.283100</td><td>7.351820</td><td>⋯</td><td>7.590200</td><td>7.658260</td><td>7.761480</td><td>7.550790</td><td>7.646170</td><td>7.47439</td><td>7.567280</td><td>7.405930</td><td>7.355990</td><td>7.322510</td></tr>
	<tr><th scope=row>ZWILCH</th><td>7.659200</td><td>7.76540</td><td>7.853750</td><td>7.896000</td><td>7.453270</td><td>7.634980</td><td>6.888660</td><td>7.419830</td><td>7.188390</td><td>7.487350</td><td>⋯</td><td>7.673230</td><td>7.754120</td><td>7.688300</td><td>7.705060</td><td>7.395470</td><td>7.20876</td><td>7.513730</td><td>7.323560</td><td>7.370120</td><td>7.576060</td></tr>
	<tr><th scope=row>ZWINT</th><td>6.204810</td><td>6.36891</td><td>6.200190</td><td>5.954580</td><td>5.790340</td><td>5.776370</td><td>6.184190</td><td>6.218960</td><td>6.216740</td><td>5.872860</td><td>⋯</td><td>6.493810</td><td>6.536290</td><td>5.738670</td><td>5.769890</td><td>5.988710</td><td>5.91986</td><td>6.364110</td><td>6.582960</td><td>6.120220</td><td>6.164980</td></tr>
	<tr><th scope=row>ZXDA</th><td>8.602590</td><td>8.62517</td><td>7.879890</td><td>8.124170</td><td>8.388830</td><td>8.302900</td><td>8.617190</td><td>8.468910</td><td>8.568600</td><td>8.414160</td><td>⋯</td><td>8.411400</td><td>8.387040</td><td>8.475430</td><td>8.546340</td><td>8.659400</td><td>8.60623</td><td>8.833910</td><td>8.671970</td><td>8.733610</td><td>8.906800</td></tr>
	<tr><th scope=row>ZXDB</th><td>7.477360</td><td>7.50947</td><td>7.987030</td><td>7.660490</td><td>7.792710</td><td>7.717490</td><td>7.928070</td><td>7.664610</td><td>7.638030</td><td>7.754150</td><td>⋯</td><td>7.601770</td><td>7.809770</td><td>7.502010</td><td>7.454340</td><td>8.245160</td><td>7.88590</td><td>7.794930</td><td>7.228620</td><td>7.990840</td><td>7.681480</td></tr>
	<tr><th scope=row>ZXDC</th><td>7.062435</td><td>6.86584</td><td>7.203055</td><td>7.284215</td><td>7.172705</td><td>7.339605</td><td>7.401120</td><td>7.318995</td><td>7.023970</td><td>7.153655</td><td>⋯</td><td>7.001105</td><td>6.981355</td><td>6.932020</td><td>6.933645</td><td>7.000900</td><td>7.08762</td><td>6.783065</td><td>7.167445</td><td>6.746320</td><td>7.024060</td></tr>
	<tr><th scope=row>ZYG11A</th><td>4.085680</td><td>4.06466</td><td>3.843450</td><td>3.706710</td><td>3.686710</td><td>3.900160</td><td>3.718290</td><td>3.901310</td><td>3.799180</td><td>3.757950</td><td>⋯</td><td>3.901450</td><td>3.974580</td><td>3.942110</td><td>3.854310</td><td>3.713220</td><td>4.09394</td><td>3.802650</td><td>3.904440</td><td>3.983960</td><td>3.981740</td></tr>
	<tr><th scope=row>ZYG11B</th><td>8.712840</td><td>8.61784</td><td>8.925850</td><td>9.013900</td><td>8.939690</td><td>8.626650</td><td>8.820020</td><td>8.613150</td><td>8.771180</td><td>8.652050</td><td>⋯</td><td>8.644580</td><td>8.760430</td><td>8.815220</td><td>8.876860</td><td>9.003580</td><td>8.71795</td><td>8.657630</td><td>8.624880</td><td>8.759080</td><td>8.521820</td></tr>
	<tr><th scope=row>ZYX</th><td>9.323210</td><td>9.48536</td><td>9.450090</td><td>9.582230</td><td>9.879900</td><td>9.477050</td><td>9.555810</td><td>9.518880</td><td>9.555710</td><td>9.497450</td><td>⋯</td><td>9.165310</td><td>9.261220</td><td>9.055100</td><td>9.066690</td><td>8.810810</td><td>8.95743</td><td>8.990700</td><td>9.322780</td><td>9.136240</td><td>9.185200</td></tr>
	<tr><th scope=row>ZZEF1</th><td>8.094900</td><td>8.17186</td><td>8.237850</td><td>8.164240</td><td>8.061520</td><td>8.082000</td><td>7.976300</td><td>8.141990</td><td>7.928100</td><td>8.083440</td><td>⋯</td><td>8.078860</td><td>8.198660</td><td>7.627300</td><td>7.766150</td><td>7.815190</td><td>7.82520</td><td>7.559610</td><td>8.104940</td><td>7.651580</td><td>7.863040</td></tr>
	<tr><th scope=row>ZZZ3</th><td>8.551840</td><td>8.59790</td><td>9.379490</td><td>9.296950</td><td>8.795220</td><td>8.837990</td><td>8.830620</td><td>8.747940</td><td>8.625360</td><td>8.773120</td><td>⋯</td><td>8.761410</td><td>8.738990</td><td>9.029910</td><td>8.993710</td><td>9.159190</td><td>8.91158</td><td>8.821660</td><td>8.680830</td><td>8.911710</td><td>8.868150</td></tr>
</tbody>
</table>




```R
modules_ME_yellow  <- data.table::fread('~/AS_HG/08.module_result/yellow-module-gene.txt')
modules_ME_brown  <- data.table::fread('~/AS_HG/08.module_result/brown-module-gene.txt')
```


```R
yellow_top100  <- modules_ME_yellow  %>% 
  arrange(-connectivity)  %>% 
   .[1:100,] %>% 
    pull(gene) %>%
     toupper()
brown_top100  <- modules_ME_brown  %>% 
  arrange(-connectivity)  %>% 
   .[1:100,] %>% 
   pull(gene) %>% 
   toupper()%>%
     toupper()
```


```R
library(GSVA)
```


```R
module_list  <- list(yellow_top100,brown_top100)
names(module_list)  <- c('Yellow_Module', 'Brown_Module')
```


```R
normalized_gset_mean
```


<table class="dataframe">
<caption>A data.frame: 23307 × 64</caption>
<thead>
	<tr><th></th><th scope=col>GSM1060117</th><th scope=col>GSM1060118</th><th scope=col>GSM1060119</th><th scope=col>GSM1060120</th><th scope=col>GSM1060121</th><th scope=col>GSM1060122</th><th scope=col>GSM1060123</th><th scope=col>GSM1060124</th><th scope=col>GSM1060125</th><th scope=col>GSM1060126</th><th scope=col>⋯</th><th scope=col>GSM1060171</th><th scope=col>GSM1060172</th><th scope=col>GSM1060173</th><th scope=col>GSM1060174</th><th scope=col>GSM1060175</th><th scope=col>GSM1060176</th><th scope=col>GSM1060177</th><th scope=col>GSM1060178</th><th scope=col>GSM1060179</th><th scope=col>GSM1060180</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>A1BG</th><td> 6.33413</td><td> 6.23136</td><td> 6.17699</td><td> 6.28800</td><td> 6.205660</td><td> 5.995310</td><td> 6.162800</td><td> 5.914780</td><td> 6.172340</td><td> 6.187750</td><td>⋯</td><td> 5.884050</td><td> 6.097230</td><td> 5.809070</td><td> 5.962630</td><td> 5.929280</td><td> 6.054380</td><td> 5.831460</td><td> 6.352940</td><td> 6.144100</td><td> 6.11485</td></tr>
	<tr><th scope=row>A1CF</th><td> 4.10072</td><td> 4.01097</td><td> 3.96745</td><td> 4.02464</td><td> 4.006060</td><td> 4.014880</td><td> 4.276880</td><td> 3.950280</td><td> 4.148390</td><td> 3.829580</td><td>⋯</td><td> 3.849670</td><td> 3.874750</td><td> 3.768630</td><td> 3.833250</td><td> 4.138480</td><td> 3.845660</td><td> 3.928350</td><td> 3.939520</td><td> 4.118220</td><td> 3.96717</td></tr>
	<tr><th scope=row>A2M</th><td>12.49600</td><td>12.26350</td><td>12.92450</td><td>12.97690</td><td>12.933100</td><td>12.958000</td><td>12.692200</td><td>12.813700</td><td>12.734900</td><td>12.929600</td><td>⋯</td><td>12.636900</td><td>12.602800</td><td>12.873400</td><td>12.867000</td><td>12.930700</td><td>12.865700</td><td>12.335300</td><td>12.131500</td><td>12.977600</td><td>12.78230</td></tr>
	<tr><th scope=row>A2ML1</th><td> 4.03661</td><td> 4.10083</td><td> 3.81421</td><td> 4.16249</td><td> 4.030510</td><td> 4.130190</td><td> 4.105710</td><td> 3.739370</td><td> 3.740730</td><td> 3.956200</td><td>⋯</td><td> 4.140920</td><td> 3.956720</td><td> 3.950060</td><td> 4.059610</td><td> 4.022620</td><td> 3.899440</td><td> 4.085170</td><td> 4.255340</td><td> 4.189320</td><td> 3.95617</td></tr>
	<tr><th scope=row>A3GALT2</th><td> 4.83119</td><td> 4.82828</td><td> 4.95262</td><td> 4.93172</td><td> 4.727280</td><td> 4.836350</td><td> 4.982820</td><td> 4.730410</td><td> 4.607110</td><td> 4.670480</td><td>⋯</td><td> 4.481430</td><td> 4.880340</td><td> 4.651580</td><td> 4.920190</td><td> 5.253800</td><td> 5.302240</td><td> 4.900220</td><td> 5.061330</td><td> 5.091100</td><td> 4.88142</td></tr>
	<tr><th scope=row>A4GALT</th><td> 7.33187</td><td> 7.26400</td><td> 6.81011</td><td> 6.81327</td><td> 6.873140</td><td> 6.717460</td><td> 6.654660</td><td> 7.011190</td><td> 7.038280</td><td> 6.839140</td><td>⋯</td><td> 6.578780</td><td> 6.307620</td><td> 6.593280</td><td> 6.963220</td><td> 6.698600</td><td> 6.508300</td><td> 7.005950</td><td> 6.676420</td><td> 6.846110</td><td> 6.55254</td></tr>
	<tr><th scope=row>A4GNT</th><td> 3.57305</td><td> 3.74132</td><td> 3.44454</td><td> 3.56300</td><td> 3.536210</td><td> 3.670310</td><td> 3.520020</td><td> 3.640490</td><td> 3.651710</td><td> 3.493550</td><td>⋯</td><td> 3.366450</td><td> 3.507460</td><td> 3.458870</td><td> 3.171950</td><td> 3.447070</td><td> 3.573210</td><td> 3.406360</td><td> 3.744170</td><td> 3.683410</td><td> 3.58425</td></tr>
	<tr><th scope=row>AAAS</th><td> 7.25019</td><td> 7.73544</td><td> 7.84654</td><td> 7.83944</td><td> 8.057220</td><td> 7.962680</td><td> 7.852890</td><td> 8.029380</td><td> 7.854730</td><td> 8.011330</td><td>⋯</td><td> 7.899490</td><td> 7.900170</td><td> 7.866410</td><td> 7.802050</td><td> 7.939140</td><td> 7.964970</td><td> 7.880880</td><td> 7.863170</td><td> 8.069680</td><td> 7.78399</td></tr>
	<tr><th scope=row>AACS</th><td> 6.52156</td><td> 6.56949</td><td> 6.50518</td><td> 6.16380</td><td> 6.794600</td><td> 6.630790</td><td> 6.713580</td><td> 6.689270</td><td> 6.841120</td><td> 6.664670</td><td>⋯</td><td> 6.558760</td><td> 6.312700</td><td> 6.389090</td><td> 6.908250</td><td> 6.357360</td><td> 6.694810</td><td> 6.294770</td><td> 6.590790</td><td> 6.427850</td><td> 6.32270</td></tr>
	<tr><th scope=row>AADAC</th><td> 2.84985</td><td> 2.82268</td><td> 2.79593</td><td> 2.78707</td><td> 2.768170</td><td> 2.703740</td><td> 3.043500</td><td> 3.035650</td><td> 2.780560</td><td> 2.614740</td><td>⋯</td><td> 2.779440</td><td> 2.899550</td><td> 2.668780</td><td> 2.890830</td><td> 2.848480</td><td> 2.701900</td><td> 2.746710</td><td> 2.709050</td><td> 2.676910</td><td> 2.95428</td></tr>
	<tr><th scope=row>AADACL2</th><td> 3.30588</td><td> 3.25592</td><td> 3.30354</td><td> 3.27986</td><td> 3.350960</td><td> 3.078130</td><td> 3.152720</td><td> 3.158270</td><td> 3.101850</td><td> 3.129590</td><td>⋯</td><td> 3.142570</td><td> 3.212370</td><td> 3.140080</td><td> 3.102040</td><td> 2.942280</td><td> 3.071950</td><td> 3.365670</td><td> 3.292630</td><td> 3.110640</td><td> 3.19635</td></tr>
	<tr><th scope=row>AADACL3</th><td> 3.90194</td><td> 4.28039</td><td> 4.01338</td><td> 4.37966</td><td> 4.308860</td><td> 4.302990</td><td> 4.719950</td><td> 4.420110</td><td> 4.291320</td><td> 4.082200</td><td>⋯</td><td> 3.989950</td><td> 4.300350</td><td> 4.447350</td><td> 4.065000</td><td> 4.244920</td><td> 4.321740</td><td> 4.029460</td><td> 4.611170</td><td> 4.362510</td><td> 4.16190</td></tr>
	<tr><th scope=row>AADACL4</th><td> 4.05240</td><td> 3.87569</td><td> 3.76920</td><td> 3.91195</td><td> 3.835640</td><td> 4.040580</td><td> 4.054060</td><td> 3.703250</td><td> 3.915150</td><td> 3.958500</td><td>⋯</td><td> 3.938130</td><td> 4.060640</td><td> 3.776560</td><td> 4.054650</td><td> 3.997450</td><td> 3.833180</td><td> 3.870620</td><td> 3.994450</td><td> 4.283060</td><td> 3.98729</td></tr>
	<tr><th scope=row>AADAT</th><td> 5.34443</td><td> 5.18716</td><td> 6.42704</td><td> 6.36462</td><td> 6.586680</td><td> 5.647800</td><td> 7.056810</td><td> 5.847700</td><td> 5.805940</td><td> 6.033930</td><td>⋯</td><td> 5.761140</td><td> 5.258240</td><td> 6.527180</td><td> 6.496410</td><td> 6.367350</td><td> 5.929410</td><td> 6.246180</td><td> 5.152180</td><td> 5.931100</td><td> 5.23681</td></tr>
	<tr><th scope=row>AAED1</th><td> 7.77722</td><td> 7.87651</td><td> 8.05750</td><td> 7.77185</td><td> 7.627100</td><td> 7.727870</td><td> 7.479700</td><td> 7.742790</td><td> 7.911490</td><td> 7.935770</td><td>⋯</td><td> 7.826170</td><td> 7.691350</td><td> 7.861950</td><td> 7.514610</td><td> 7.477840</td><td> 7.398790</td><td> 7.851200</td><td> 7.622080</td><td> 7.656300</td><td> 7.34300</td></tr>
	<tr><th scope=row>AAGAB</th><td> 7.85477</td><td> 7.70956</td><td> 7.19151</td><td> 7.22967</td><td> 7.488280</td><td> 7.314570</td><td> 7.575450</td><td> 7.494540</td><td> 7.642140</td><td> 7.432660</td><td>⋯</td><td> 7.758110</td><td> 7.923780</td><td> 7.718050</td><td> 7.748940</td><td> 7.570000</td><td> 7.654970</td><td> 8.127740</td><td> 7.788350</td><td> 7.501440</td><td> 7.84698</td></tr>
	<tr><th scope=row>AAK1</th><td> 8.22717</td><td> 8.31815</td><td> 8.74060</td><td> 8.66952</td><td> 8.095265</td><td> 8.313355</td><td> 8.128905</td><td> 8.447590</td><td> 8.398370</td><td> 8.440045</td><td>⋯</td><td> 8.375805</td><td> 8.456655</td><td> 8.012455</td><td> 8.164350</td><td> 7.946180</td><td> 8.132335</td><td> 7.754315</td><td> 8.327455</td><td> 7.955945</td><td> 8.44443</td></tr>
	<tr><th scope=row>AAMDC</th><td> 5.65657</td><td> 5.57825</td><td> 6.05968</td><td> 6.35897</td><td> 6.521880</td><td> 6.125580</td><td> 6.337650</td><td> 6.076220</td><td> 6.213170</td><td> 6.116430</td><td>⋯</td><td> 6.518710</td><td> 6.249870</td><td> 6.586000</td><td> 6.512410</td><td> 6.691830</td><td> 6.426220</td><td> 6.431450</td><td> 5.910120</td><td> 6.358170</td><td> 6.15294</td></tr>
	<tr><th scope=row>AAMP</th><td> 8.18316</td><td> 8.33303</td><td> 7.66461</td><td> 7.71843</td><td> 8.297810</td><td> 7.988820</td><td> 8.147270</td><td> 8.171650</td><td> 8.147440</td><td> 8.082130</td><td>⋯</td><td> 7.850910</td><td> 7.861290</td><td> 8.165260</td><td> 7.989590</td><td> 7.910050</td><td> 7.907020</td><td> 8.024960</td><td> 7.811020</td><td> 8.089790</td><td> 7.70452</td></tr>
	<tr><th scope=row>AANAT</th><td> 5.63924</td><td> 5.65591</td><td> 5.41210</td><td> 5.43153</td><td> 5.489600</td><td> 5.380810</td><td> 5.680120</td><td> 5.517040</td><td> 5.584820</td><td> 5.356170</td><td>⋯</td><td> 5.325340</td><td> 5.715690</td><td> 5.265600</td><td> 5.503580</td><td> 5.343190</td><td> 5.357820</td><td> 5.126700</td><td> 5.644050</td><td> 5.284600</td><td> 5.54677</td></tr>
	<tr><th scope=row>AAR2</th><td> 7.46769</td><td> 7.46391</td><td> 7.18075</td><td> 7.25863</td><td> 7.629600</td><td> 7.250830</td><td> 7.409180</td><td> 7.372140</td><td> 7.532570</td><td> 7.334860</td><td>⋯</td><td> 7.240740</td><td> 7.304730</td><td> 7.674550</td><td> 7.463160</td><td> 7.572510</td><td> 7.367340</td><td> 7.417170</td><td> 7.451330</td><td> 7.471430</td><td> 7.30658</td></tr>
	<tr><th scope=row>AARD</th><td> 5.08894</td><td> 5.17774</td><td> 5.13358</td><td> 5.13923</td><td> 5.215750</td><td> 4.884590</td><td> 5.094010</td><td> 4.856330</td><td> 4.774160</td><td> 4.850970</td><td>⋯</td><td> 4.852540</td><td> 4.850680</td><td> 4.903870</td><td> 5.351960</td><td> 5.208560</td><td> 4.954350</td><td> 5.153520</td><td> 4.940130</td><td> 5.100770</td><td> 5.20342</td></tr>
	<tr><th scope=row>AARS</th><td> 9.19336</td><td> 9.39812</td><td> 8.85656</td><td> 8.95554</td><td> 9.076860</td><td> 9.148440</td><td> 9.099850</td><td> 9.183370</td><td> 9.190180</td><td> 9.296830</td><td>⋯</td><td> 9.057630</td><td> 9.122210</td><td> 9.228770</td><td> 9.127030</td><td> 9.001740</td><td> 8.955210</td><td> 8.941960</td><td> 9.033640</td><td> 8.994620</td><td> 8.97487</td></tr>
	<tr><th scope=row>AARS2</th><td> 6.23486</td><td> 6.39646</td><td> 5.88791</td><td> 5.67745</td><td> 6.065220</td><td> 6.385780</td><td> 6.105200</td><td> 6.082310</td><td> 6.120040</td><td> 6.041290</td><td>⋯</td><td> 5.824280</td><td> 5.862870</td><td> 5.701840</td><td> 5.591840</td><td> 5.890080</td><td> 5.823160</td><td> 5.595190</td><td> 6.080590</td><td> 6.213680</td><td> 6.11416</td></tr>
	<tr><th scope=row>AASDH</th><td> 6.49165</td><td> 6.90503</td><td> 7.15174</td><td> 7.24475</td><td> 6.874940</td><td> 7.125700</td><td> 6.985320</td><td> 7.014320</td><td> 6.732850</td><td> 7.013900</td><td>⋯</td><td> 7.198320</td><td> 7.197440</td><td> 7.666840</td><td> 7.589570</td><td> 7.365940</td><td> 7.256180</td><td> 7.423760</td><td> 6.511200</td><td> 7.169750</td><td> 6.91563</td></tr>
	<tr><th scope=row>AASDHPPT</th><td> 5.49130</td><td> 5.50187</td><td> 5.44725</td><td> 5.23175</td><td> 5.489490</td><td> 5.362130</td><td> 5.431920</td><td> 5.311595</td><td> 5.343005</td><td> 5.442175</td><td>⋯</td><td> 5.463835</td><td> 5.524950</td><td> 5.574115</td><td> 5.531265</td><td> 5.399075</td><td> 5.289875</td><td> 5.614125</td><td> 5.227605</td><td> 5.403960</td><td> 5.43015</td></tr>
	<tr><th scope=row>AASS</th><td> 8.08306</td><td> 7.61825</td><td> 9.08592</td><td> 9.14799</td><td> 8.994880</td><td> 8.893730</td><td> 9.194020</td><td> 8.624360</td><td> 8.232910</td><td> 8.663470</td><td>⋯</td><td> 8.520820</td><td> 7.958940</td><td> 9.166090</td><td> 9.234190</td><td> 8.991620</td><td> 8.765830</td><td> 8.637630</td><td> 7.642380</td><td> 9.034170</td><td> 8.59184</td></tr>
	<tr><th scope=row>AATF</th><td> 8.20921</td><td> 8.33213</td><td> 8.10343</td><td> 7.93821</td><td> 8.259060</td><td> 8.302820</td><td> 8.241540</td><td> 8.213570</td><td> 8.316550</td><td> 8.342720</td><td>⋯</td><td> 8.073900</td><td> 8.125390</td><td> 8.173360</td><td> 8.115940</td><td> 8.202090</td><td> 8.036140</td><td> 8.214020</td><td> 8.237470</td><td> 8.103410</td><td> 8.09945</td></tr>
	<tr><th scope=row>AATK</th><td> 6.11253</td><td> 6.35266</td><td> 5.79944</td><td> 5.87782</td><td> 6.227190</td><td> 6.051150</td><td> 6.081640</td><td> 5.927340</td><td> 6.070090</td><td> 5.818400</td><td>⋯</td><td> 5.761910</td><td> 5.922940</td><td> 5.669620</td><td> 5.897280</td><td> 5.988400</td><td> 5.985300</td><td> 5.703370</td><td> 6.300550</td><td> 6.284870</td><td> 6.26960</td></tr>
	<tr><th scope=row>AATK-AS1</th><td> 6.35996</td><td> 6.45766</td><td> 6.13680</td><td> 5.98517</td><td> 6.162350</td><td> 6.092930</td><td> 6.216260</td><td> 6.006140</td><td> 6.158240</td><td> 6.033200</td><td>⋯</td><td> 5.921370</td><td> 6.126150</td><td> 5.867100</td><td> 6.149820</td><td> 6.224110</td><td> 6.028670</td><td> 6.050850</td><td> 6.381870</td><td> 6.358470</td><td> 6.29537</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>ZSCAN23</th><td>4.739290</td><td>4.44092</td><td>5.322820</td><td>5.360020</td><td>4.829890</td><td>4.831710</td><td>4.994780</td><td>4.952940</td><td>4.926850</td><td>5.052570</td><td>⋯</td><td>4.519370</td><td>4.361260</td><td>4.795260</td><td>4.951550</td><td>5.256740</td><td>4.99614</td><td>4.488520</td><td>4.642920</td><td>4.940440</td><td>4.441750</td></tr>
	<tr><th scope=row>ZSCAN25</th><td>7.255930</td><td>6.98109</td><td>7.557610</td><td>7.452420</td><td>7.342360</td><td>7.340570</td><td>7.306950</td><td>7.422280</td><td>7.279090</td><td>7.353910</td><td>⋯</td><td>6.920140</td><td>6.724410</td><td>6.944710</td><td>7.144050</td><td>7.202520</td><td>7.12622</td><td>6.845080</td><td>7.080030</td><td>6.941390</td><td>7.280240</td></tr>
	<tr><th scope=row>ZSCAN26</th><td>6.541880</td><td>6.70649</td><td>7.292780</td><td>7.093080</td><td>7.202300</td><td>7.112510</td><td>7.353800</td><td>6.982480</td><td>7.042130</td><td>7.019270</td><td>⋯</td><td>7.211700</td><td>7.222990</td><td>7.532700</td><td>7.452160</td><td>7.654810</td><td>7.23261</td><td>7.065450</td><td>6.937430</td><td>7.157190</td><td>7.225920</td></tr>
	<tr><th scope=row>ZSCAN29</th><td>7.341850</td><td>7.30534</td><td>7.587300</td><td>7.445370</td><td>7.597690</td><td>7.352240</td><td>7.439310</td><td>7.267930</td><td>7.208070</td><td>7.319900</td><td>⋯</td><td>7.166810</td><td>7.235580</td><td>7.549300</td><td>7.264300</td><td>7.268140</td><td>7.11719</td><td>7.236160</td><td>7.386040</td><td>7.251770</td><td>7.457940</td></tr>
	<tr><th scope=row>ZSCAN30</th><td>5.442550</td><td>5.32821</td><td>5.629910</td><td>5.853860</td><td>5.608910</td><td>5.523470</td><td>5.584390</td><td>5.652740</td><td>5.576710</td><td>5.622540</td><td>⋯</td><td>5.460460</td><td>5.318840</td><td>5.777450</td><td>5.556350</td><td>5.556090</td><td>5.72232</td><td>5.319860</td><td>5.590800</td><td>5.708400</td><td>5.304390</td></tr>
	<tr><th scope=row>ZSCAN31</th><td>5.055560</td><td>4.93817</td><td>5.596350</td><td>5.517290</td><td>5.702040</td><td>5.313130</td><td>5.919250</td><td>5.272140</td><td>5.461130</td><td>6.033240</td><td>⋯</td><td>5.273260</td><td>5.061410</td><td>5.989300</td><td>5.607330</td><td>5.964670</td><td>5.82656</td><td>5.777770</td><td>5.628260</td><td>5.478270</td><td>5.147750</td></tr>
	<tr><th scope=row>ZSCAN32</th><td>6.275010</td><td>6.41561</td><td>6.262170</td><td>6.061880</td><td>6.232320</td><td>6.119590</td><td>6.337180</td><td>6.236720</td><td>6.304530</td><td>6.374180</td><td>⋯</td><td>5.945780</td><td>6.179630</td><td>6.296460</td><td>6.145740</td><td>6.235130</td><td>6.05691</td><td>6.229960</td><td>6.093360</td><td>6.300080</td><td>6.140160</td></tr>
	<tr><th scope=row>ZSCAN4</th><td>3.075580</td><td>3.03255</td><td>2.733240</td><td>2.902610</td><td>2.848370</td><td>2.658400</td><td>2.869320</td><td>2.916180</td><td>3.083690</td><td>2.915580</td><td>⋯</td><td>2.856440</td><td>3.177350</td><td>2.756660</td><td>3.108830</td><td>3.005950</td><td>2.87814</td><td>2.897430</td><td>2.999380</td><td>3.024120</td><td>3.216140</td></tr>
	<tr><th scope=row>ZSCAN5A</th><td>6.010695</td><td>6.25011</td><td>5.819710</td><td>5.772685</td><td>5.862080</td><td>5.970995</td><td>5.817635</td><td>5.933545</td><td>6.100115</td><td>5.942610</td><td>⋯</td><td>5.716315</td><td>6.054235</td><td>5.790415</td><td>5.823925</td><td>5.736070</td><td>5.78424</td><td>5.702880</td><td>5.950375</td><td>5.993690</td><td>5.935370</td></tr>
	<tr><th scope=row>ZSCAN5B</th><td>3.847550</td><td>4.58524</td><td>3.580100</td><td>3.909380</td><td>3.924590</td><td>3.874920</td><td>4.009860</td><td>3.881090</td><td>3.686350</td><td>3.907280</td><td>⋯</td><td>4.030190</td><td>4.051830</td><td>3.692560</td><td>3.575180</td><td>3.669300</td><td>3.59561</td><td>3.950280</td><td>3.922240</td><td>3.656240</td><td>3.721820</td></tr>
	<tr><th scope=row>ZSCAN9</th><td>6.046260</td><td>5.96850</td><td>6.302460</td><td>6.477780</td><td>6.464210</td><td>6.432070</td><td>6.477620</td><td>6.403120</td><td>6.227690</td><td>6.463900</td><td>⋯</td><td>6.660720</td><td>6.572350</td><td>6.975170</td><td>6.543870</td><td>6.579070</td><td>6.42234</td><td>6.334720</td><td>6.553440</td><td>6.600680</td><td>6.593470</td></tr>
	<tr><th scope=row>ZSWIM1</th><td>6.114610</td><td>6.21384</td><td>6.390560</td><td>6.134040</td><td>6.421410</td><td>6.146060</td><td>6.295100</td><td>6.366980</td><td>6.387090</td><td>6.506010</td><td>⋯</td><td>6.578690</td><td>6.463700</td><td>6.612010</td><td>6.526830</td><td>6.671320</td><td>6.77247</td><td>6.365320</td><td>6.407550</td><td>6.626760</td><td>6.563900</td></tr>
	<tr><th scope=row>ZSWIM2</th><td>2.699320</td><td>2.79870</td><td>2.661550</td><td>2.695030</td><td>2.798440</td><td>2.755660</td><td>2.598900</td><td>2.535180</td><td>2.747370</td><td>2.645670</td><td>⋯</td><td>2.691240</td><td>2.592970</td><td>2.495090</td><td>2.660870</td><td>2.677300</td><td>2.67362</td><td>2.700920</td><td>2.817000</td><td>2.717310</td><td>2.668020</td></tr>
	<tr><th scope=row>ZSWIM3</th><td>4.804540</td><td>4.82555</td><td>4.729410</td><td>4.498370</td><td>5.003850</td><td>4.844300</td><td>4.690370</td><td>4.826390</td><td>4.814820</td><td>4.969180</td><td>⋯</td><td>4.779870</td><td>4.757850</td><td>4.757170</td><td>4.862460</td><td>4.706120</td><td>4.85834</td><td>4.926400</td><td>4.985530</td><td>4.819000</td><td>4.754920</td></tr>
	<tr><th scope=row>ZSWIM4</th><td>6.274790</td><td>6.67088</td><td>6.135360</td><td>6.159800</td><td>6.149640</td><td>6.312880</td><td>6.313830</td><td>6.336560</td><td>6.431320</td><td>6.282780</td><td>⋯</td><td>5.795520</td><td>6.276530</td><td>5.821120</td><td>5.958720</td><td>5.943340</td><td>5.88275</td><td>5.930330</td><td>6.672070</td><td>6.145630</td><td>6.184300</td></tr>
	<tr><th scope=row>ZSWIM5</th><td>5.813450</td><td>5.42529</td><td>6.104120</td><td>6.035930</td><td>6.081840</td><td>5.851010</td><td>6.348580</td><td>6.054770</td><td>6.160560</td><td>5.933700</td><td>⋯</td><td>5.879680</td><td>5.915860</td><td>6.120320</td><td>6.211970</td><td>6.348330</td><td>6.14024</td><td>5.826440</td><td>5.954130</td><td>6.021740</td><td>5.970800</td></tr>
	<tr><th scope=row>ZSWIM6</th><td>9.017940</td><td>9.25979</td><td>9.153690</td><td>9.113050</td><td>8.846030</td><td>9.013890</td><td>8.494060</td><td>9.051230</td><td>8.789680</td><td>9.071290</td><td>⋯</td><td>9.152630</td><td>9.486710</td><td>8.881170</td><td>8.956680</td><td>8.675820</td><td>8.86422</td><td>8.938520</td><td>9.495070</td><td>8.936340</td><td>9.556640</td></tr>
	<tr><th scope=row>ZSWIM8</th><td>6.630595</td><td>6.48948</td><td>6.233210</td><td>6.300660</td><td>6.756585</td><td>6.452325</td><td>6.538550</td><td>6.518880</td><td>6.711720</td><td>6.640125</td><td>⋯</td><td>6.336675</td><td>6.353280</td><td>6.231415</td><td>6.348865</td><td>6.464705</td><td>6.46188</td><td>6.332295</td><td>6.551470</td><td>6.507205</td><td>6.569705</td></tr>
	<tr><th scope=row>ZUFSP</th><td>5.507730</td><td>5.65530</td><td>5.917350</td><td>5.797510</td><td>5.712160</td><td>5.607790</td><td>5.680950</td><td>5.699750</td><td>5.444790</td><td>5.672050</td><td>⋯</td><td>6.261170</td><td>6.048710</td><td>5.804440</td><td>6.011210</td><td>5.742370</td><td>5.57482</td><td>6.453930</td><td>5.775680</td><td>5.550680</td><td>6.104750</td></tr>
	<tr><th scope=row>ZW10</th><td>7.345310</td><td>7.47667</td><td>7.603950</td><td>7.599710</td><td>7.459620</td><td>7.410450</td><td>7.327460</td><td>7.317760</td><td>7.283100</td><td>7.351820</td><td>⋯</td><td>7.590200</td><td>7.658260</td><td>7.761480</td><td>7.550790</td><td>7.646170</td><td>7.47439</td><td>7.567280</td><td>7.405930</td><td>7.355990</td><td>7.322510</td></tr>
	<tr><th scope=row>ZWILCH</th><td>7.659200</td><td>7.76540</td><td>7.853750</td><td>7.896000</td><td>7.453270</td><td>7.634980</td><td>6.888660</td><td>7.419830</td><td>7.188390</td><td>7.487350</td><td>⋯</td><td>7.673230</td><td>7.754120</td><td>7.688300</td><td>7.705060</td><td>7.395470</td><td>7.20876</td><td>7.513730</td><td>7.323560</td><td>7.370120</td><td>7.576060</td></tr>
	<tr><th scope=row>ZWINT</th><td>6.204810</td><td>6.36891</td><td>6.200190</td><td>5.954580</td><td>5.790340</td><td>5.776370</td><td>6.184190</td><td>6.218960</td><td>6.216740</td><td>5.872860</td><td>⋯</td><td>6.493810</td><td>6.536290</td><td>5.738670</td><td>5.769890</td><td>5.988710</td><td>5.91986</td><td>6.364110</td><td>6.582960</td><td>6.120220</td><td>6.164980</td></tr>
	<tr><th scope=row>ZXDA</th><td>8.602590</td><td>8.62517</td><td>7.879890</td><td>8.124170</td><td>8.388830</td><td>8.302900</td><td>8.617190</td><td>8.468910</td><td>8.568600</td><td>8.414160</td><td>⋯</td><td>8.411400</td><td>8.387040</td><td>8.475430</td><td>8.546340</td><td>8.659400</td><td>8.60623</td><td>8.833910</td><td>8.671970</td><td>8.733610</td><td>8.906800</td></tr>
	<tr><th scope=row>ZXDB</th><td>7.477360</td><td>7.50947</td><td>7.987030</td><td>7.660490</td><td>7.792710</td><td>7.717490</td><td>7.928070</td><td>7.664610</td><td>7.638030</td><td>7.754150</td><td>⋯</td><td>7.601770</td><td>7.809770</td><td>7.502010</td><td>7.454340</td><td>8.245160</td><td>7.88590</td><td>7.794930</td><td>7.228620</td><td>7.990840</td><td>7.681480</td></tr>
	<tr><th scope=row>ZXDC</th><td>7.062435</td><td>6.86584</td><td>7.203055</td><td>7.284215</td><td>7.172705</td><td>7.339605</td><td>7.401120</td><td>7.318995</td><td>7.023970</td><td>7.153655</td><td>⋯</td><td>7.001105</td><td>6.981355</td><td>6.932020</td><td>6.933645</td><td>7.000900</td><td>7.08762</td><td>6.783065</td><td>7.167445</td><td>6.746320</td><td>7.024060</td></tr>
	<tr><th scope=row>ZYG11A</th><td>4.085680</td><td>4.06466</td><td>3.843450</td><td>3.706710</td><td>3.686710</td><td>3.900160</td><td>3.718290</td><td>3.901310</td><td>3.799180</td><td>3.757950</td><td>⋯</td><td>3.901450</td><td>3.974580</td><td>3.942110</td><td>3.854310</td><td>3.713220</td><td>4.09394</td><td>3.802650</td><td>3.904440</td><td>3.983960</td><td>3.981740</td></tr>
	<tr><th scope=row>ZYG11B</th><td>8.712840</td><td>8.61784</td><td>8.925850</td><td>9.013900</td><td>8.939690</td><td>8.626650</td><td>8.820020</td><td>8.613150</td><td>8.771180</td><td>8.652050</td><td>⋯</td><td>8.644580</td><td>8.760430</td><td>8.815220</td><td>8.876860</td><td>9.003580</td><td>8.71795</td><td>8.657630</td><td>8.624880</td><td>8.759080</td><td>8.521820</td></tr>
	<tr><th scope=row>ZYX</th><td>9.323210</td><td>9.48536</td><td>9.450090</td><td>9.582230</td><td>9.879900</td><td>9.477050</td><td>9.555810</td><td>9.518880</td><td>9.555710</td><td>9.497450</td><td>⋯</td><td>9.165310</td><td>9.261220</td><td>9.055100</td><td>9.066690</td><td>8.810810</td><td>8.95743</td><td>8.990700</td><td>9.322780</td><td>9.136240</td><td>9.185200</td></tr>
	<tr><th scope=row>ZZEF1</th><td>8.094900</td><td>8.17186</td><td>8.237850</td><td>8.164240</td><td>8.061520</td><td>8.082000</td><td>7.976300</td><td>8.141990</td><td>7.928100</td><td>8.083440</td><td>⋯</td><td>8.078860</td><td>8.198660</td><td>7.627300</td><td>7.766150</td><td>7.815190</td><td>7.82520</td><td>7.559610</td><td>8.104940</td><td>7.651580</td><td>7.863040</td></tr>
	<tr><th scope=row>ZZZ3</th><td>8.551840</td><td>8.59790</td><td>9.379490</td><td>9.296950</td><td>8.795220</td><td>8.837990</td><td>8.830620</td><td>8.747940</td><td>8.625360</td><td>8.773120</td><td>⋯</td><td>8.761410</td><td>8.738990</td><td>9.029910</td><td>8.993710</td><td>9.159190</td><td>8.91158</td><td>8.821660</td><td>8.680830</td><td>8.911710</td><td>8.868150</td></tr>
</tbody>
</table>




```R
gsva_matrix<- gsva(as.matrix(normalized_gset_mean), module_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
```

    Estimating ssGSEA scores for 2 gene sets.
      |======================================================================| 100%
    



```R
library(ggpubr)
```


```R
gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list) 
```


<table class="dataframe">
<caption>A data.frame: 64 × 4</caption>
<thead>
	<tr><th scope=col>sample</th><th scope=col>Yellow_Module</th><th scope=col>Brown_Module</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM1060117</td><td>2.413568</td><td>2.930204</td><td>MIT</td></tr>
	<tr><td>GSM1060118</td><td>2.473144</td><td>3.061218</td><td>ATH</td></tr>
	<tr><td>GSM1060119</td><td>2.379032</td><td>2.577177</td><td>MIT</td></tr>
	<tr><td>GSM1060120</td><td>2.419105</td><td>2.593116</td><td>ATH</td></tr>
	<tr><td>GSM1060121</td><td>2.344249</td><td>2.362440</td><td>MIT</td></tr>
	<tr><td>GSM1060122</td><td>2.448625</td><td>2.691610</td><td>ATH</td></tr>
	<tr><td>GSM1060123</td><td>2.212691</td><td>2.061218</td><td>MIT</td></tr>
	<tr><td>GSM1060124</td><td>2.476706</td><td>2.779070</td><td>ATH</td></tr>
	<tr><td>GSM1060125</td><td>2.421892</td><td>2.709830</td><td>MIT</td></tr>
	<tr><td>GSM1060126</td><td>2.404670</td><td>2.700017</td><td>ATH</td></tr>
	<tr><td>GSM1060127</td><td>2.456928</td><td>2.505985</td><td>MIT</td></tr>
	<tr><td>GSM1060128</td><td>2.637750</td><td>3.041844</td><td>ATH</td></tr>
	<tr><td>GSM1060129</td><td>2.499054</td><td>2.625484</td><td>MIT</td></tr>
	<tr><td>GSM1060130</td><td>2.516723</td><td>2.968342</td><td>ATH</td></tr>
	<tr><td>GSM1060131</td><td>2.521294</td><td>2.716727</td><td>MIT</td></tr>
	<tr><td>GSM1060132</td><td>2.677688</td><td>3.031473</td><td>ATH</td></tr>
	<tr><td>GSM1060133</td><td>2.291868</td><td>2.419839</td><td>MIT</td></tr>
	<tr><td>GSM1060134</td><td>2.447359</td><td>2.743517</td><td>ATH</td></tr>
	<tr><td>GSM1060135</td><td>2.340889</td><td>2.513890</td><td>MIT</td></tr>
	<tr><td>GSM1060136</td><td>2.463149</td><td>2.843104</td><td>ATH</td></tr>
	<tr><td>GSM1060137</td><td>2.322210</td><td>2.357709</td><td>MIT</td></tr>
	<tr><td>GSM1060138</td><td>2.426669</td><td>2.756912</td><td>ATH</td></tr>
	<tr><td>GSM1060139</td><td>2.472265</td><td>2.987640</td><td>MIT</td></tr>
	<tr><td>GSM1060140</td><td>2.462679</td><td>2.904020</td><td>ATH</td></tr>
	<tr><td>GSM1060141</td><td>2.434314</td><td>2.826986</td><td>MIT</td></tr>
	<tr><td>GSM1060142</td><td>2.566418</td><td>2.966012</td><td>ATH</td></tr>
	<tr><td>GSM1060143</td><td>2.292374</td><td>2.580535</td><td>MIT</td></tr>
	<tr><td>GSM1060144</td><td>2.494033</td><td>3.001745</td><td>ATH</td></tr>
	<tr><td>GSM1060145</td><td>2.409568</td><td>2.669450</td><td>MIT</td></tr>
	<tr><td>GSM1060146</td><td>2.628845</td><td>3.027742</td><td>ATH</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>GSM1060151</td><td>2.394780</td><td>2.767477</td><td>MIT</td></tr>
	<tr><td>GSM1060152</td><td>2.508729</td><td>2.873901</td><td>ATH</td></tr>
	<tr><td>GSM1060153</td><td>2.361131</td><td>2.388048</td><td>MIT</td></tr>
	<tr><td>GSM1060154</td><td>2.499367</td><td>2.791676</td><td>ATH</td></tr>
	<tr><td>GSM1060155</td><td>2.298204</td><td>2.478005</td><td>MIT</td></tr>
	<tr><td>GSM1060156</td><td>2.456052</td><td>2.875328</td><td>ATH</td></tr>
	<tr><td>GSM1060157</td><td>2.396072</td><td>2.699308</td><td>MIT</td></tr>
	<tr><td>GSM1060158</td><td>2.401507</td><td>2.728450</td><td>ATH</td></tr>
	<tr><td>GSM1060159</td><td>2.384107</td><td>2.699930</td><td>MIT</td></tr>
	<tr><td>GSM1060160</td><td>2.482489</td><td>2.797697</td><td>ATH</td></tr>
	<tr><td>GSM1060161</td><td>2.219465</td><td>2.362664</td><td>MIT</td></tr>
	<tr><td>GSM1060162</td><td>2.207843</td><td>2.427069</td><td>ATH</td></tr>
	<tr><td>GSM1060163</td><td>2.232270</td><td>2.473000</td><td>MIT</td></tr>
	<tr><td>GSM1060164</td><td>2.410908</td><td>2.740741</td><td>ATH</td></tr>
	<tr><td>GSM1060165</td><td>2.270437</td><td>2.436001</td><td>MIT</td></tr>
	<tr><td>GSM1060166</td><td>2.541955</td><td>2.915907</td><td>ATH</td></tr>
	<tr><td>GSM1060167</td><td>2.319819</td><td>2.558937</td><td>MIT</td></tr>
	<tr><td>GSM1060168</td><td>2.367182</td><td>2.577066</td><td>ATH</td></tr>
	<tr><td>GSM1060169</td><td>2.299078</td><td>2.198046</td><td>MIT</td></tr>
	<tr><td>GSM1060170</td><td>2.331730</td><td>2.674969</td><td>ATH</td></tr>
	<tr><td>GSM1060171</td><td>2.384890</td><td>2.696406</td><td>MIT</td></tr>
	<tr><td>GSM1060172</td><td>2.416842</td><td>2.923808</td><td>ATH</td></tr>
	<tr><td>GSM1060173</td><td>2.175337</td><td>2.410616</td><td>MIT</td></tr>
	<tr><td>GSM1060174</td><td>2.266244</td><td>2.478458</td><td>ATH</td></tr>
	<tr><td>GSM1060175</td><td>2.260628</td><td>2.277613</td><td>MIT</td></tr>
	<tr><td>GSM1060176</td><td>2.331283</td><td>2.570472</td><td>ATH</td></tr>
	<tr><td>GSM1060177</td><td>2.269510</td><td>2.602500</td><td>MIT</td></tr>
	<tr><td>GSM1060178</td><td>2.530203</td><td>2.930627</td><td>ATH</td></tr>
	<tr><td>GSM1060179</td><td>2.318810</td><td>2.558938</td><td>MIT</td></tr>
	<tr><td>GSM1060180</td><td>2.407185</td><td>2.729593</td><td>ATH</td></tr>
</tbody>
</table>




```R
p_yellow  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       ggboxplot(x = 'group', y = 'Yellow_Module', ylab = 'GSVA Score',title = 'Yellow Module',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('MIT','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_yellow
```


    
![png](Step5_Two_Group_Compare_GSE43292_files/Step5_Two_Group_Compare_GSE43292_35_0.png)
    



```R
gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  
```


<table class="dataframe">
<caption>A data.frame: 64 × 4</caption>
<thead>
	<tr><th scope=col>sample</th><th scope=col>Yellow_Module</th><th scope=col>Brown_Module</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM1060117</td><td>2.413568</td><td>2.930204</td><td>MIT</td></tr>
	<tr><td>GSM1060118</td><td>2.473144</td><td>3.061218</td><td>ATH</td></tr>
	<tr><td>GSM1060119</td><td>2.379032</td><td>2.577177</td><td>MIT</td></tr>
	<tr><td>GSM1060120</td><td>2.419105</td><td>2.593116</td><td>ATH</td></tr>
	<tr><td>GSM1060121</td><td>2.344249</td><td>2.362440</td><td>MIT</td></tr>
	<tr><td>GSM1060122</td><td>2.448625</td><td>2.691610</td><td>ATH</td></tr>
	<tr><td>GSM1060123</td><td>2.212691</td><td>2.061218</td><td>MIT</td></tr>
	<tr><td>GSM1060124</td><td>2.476706</td><td>2.779070</td><td>ATH</td></tr>
	<tr><td>GSM1060125</td><td>2.421892</td><td>2.709830</td><td>MIT</td></tr>
	<tr><td>GSM1060126</td><td>2.404670</td><td>2.700017</td><td>ATH</td></tr>
	<tr><td>GSM1060127</td><td>2.456928</td><td>2.505985</td><td>MIT</td></tr>
	<tr><td>GSM1060128</td><td>2.637750</td><td>3.041844</td><td>ATH</td></tr>
	<tr><td>GSM1060129</td><td>2.499054</td><td>2.625484</td><td>MIT</td></tr>
	<tr><td>GSM1060130</td><td>2.516723</td><td>2.968342</td><td>ATH</td></tr>
	<tr><td>GSM1060131</td><td>2.521294</td><td>2.716727</td><td>MIT</td></tr>
	<tr><td>GSM1060132</td><td>2.677688</td><td>3.031473</td><td>ATH</td></tr>
	<tr><td>GSM1060133</td><td>2.291868</td><td>2.419839</td><td>MIT</td></tr>
	<tr><td>GSM1060134</td><td>2.447359</td><td>2.743517</td><td>ATH</td></tr>
	<tr><td>GSM1060135</td><td>2.340889</td><td>2.513890</td><td>MIT</td></tr>
	<tr><td>GSM1060136</td><td>2.463149</td><td>2.843104</td><td>ATH</td></tr>
	<tr><td>GSM1060137</td><td>2.322210</td><td>2.357709</td><td>MIT</td></tr>
	<tr><td>GSM1060138</td><td>2.426669</td><td>2.756912</td><td>ATH</td></tr>
	<tr><td>GSM1060139</td><td>2.472265</td><td>2.987640</td><td>MIT</td></tr>
	<tr><td>GSM1060140</td><td>2.462679</td><td>2.904020</td><td>ATH</td></tr>
	<tr><td>GSM1060141</td><td>2.434314</td><td>2.826986</td><td>MIT</td></tr>
	<tr><td>GSM1060142</td><td>2.566418</td><td>2.966012</td><td>ATH</td></tr>
	<tr><td>GSM1060143</td><td>2.292374</td><td>2.580535</td><td>MIT</td></tr>
	<tr><td>GSM1060144</td><td>2.494033</td><td>3.001745</td><td>ATH</td></tr>
	<tr><td>GSM1060145</td><td>2.409568</td><td>2.669450</td><td>MIT</td></tr>
	<tr><td>GSM1060146</td><td>2.628845</td><td>3.027742</td><td>ATH</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>GSM1060151</td><td>2.394780</td><td>2.767477</td><td>MIT</td></tr>
	<tr><td>GSM1060152</td><td>2.508729</td><td>2.873901</td><td>ATH</td></tr>
	<tr><td>GSM1060153</td><td>2.361131</td><td>2.388048</td><td>MIT</td></tr>
	<tr><td>GSM1060154</td><td>2.499367</td><td>2.791676</td><td>ATH</td></tr>
	<tr><td>GSM1060155</td><td>2.298204</td><td>2.478005</td><td>MIT</td></tr>
	<tr><td>GSM1060156</td><td>2.456052</td><td>2.875328</td><td>ATH</td></tr>
	<tr><td>GSM1060157</td><td>2.396072</td><td>2.699308</td><td>MIT</td></tr>
	<tr><td>GSM1060158</td><td>2.401507</td><td>2.728450</td><td>ATH</td></tr>
	<tr><td>GSM1060159</td><td>2.384107</td><td>2.699930</td><td>MIT</td></tr>
	<tr><td>GSM1060160</td><td>2.482489</td><td>2.797697</td><td>ATH</td></tr>
	<tr><td>GSM1060161</td><td>2.219465</td><td>2.362664</td><td>MIT</td></tr>
	<tr><td>GSM1060162</td><td>2.207843</td><td>2.427069</td><td>ATH</td></tr>
	<tr><td>GSM1060163</td><td>2.232270</td><td>2.473000</td><td>MIT</td></tr>
	<tr><td>GSM1060164</td><td>2.410908</td><td>2.740741</td><td>ATH</td></tr>
	<tr><td>GSM1060165</td><td>2.270437</td><td>2.436001</td><td>MIT</td></tr>
	<tr><td>GSM1060166</td><td>2.541955</td><td>2.915907</td><td>ATH</td></tr>
	<tr><td>GSM1060167</td><td>2.319819</td><td>2.558937</td><td>MIT</td></tr>
	<tr><td>GSM1060168</td><td>2.367182</td><td>2.577066</td><td>ATH</td></tr>
	<tr><td>GSM1060169</td><td>2.299078</td><td>2.198046</td><td>MIT</td></tr>
	<tr><td>GSM1060170</td><td>2.331730</td><td>2.674969</td><td>ATH</td></tr>
	<tr><td>GSM1060171</td><td>2.384890</td><td>2.696406</td><td>MIT</td></tr>
	<tr><td>GSM1060172</td><td>2.416842</td><td>2.923808</td><td>ATH</td></tr>
	<tr><td>GSM1060173</td><td>2.175337</td><td>2.410616</td><td>MIT</td></tr>
	<tr><td>GSM1060174</td><td>2.266244</td><td>2.478458</td><td>ATH</td></tr>
	<tr><td>GSM1060175</td><td>2.260628</td><td>2.277613</td><td>MIT</td></tr>
	<tr><td>GSM1060176</td><td>2.331283</td><td>2.570472</td><td>ATH</td></tr>
	<tr><td>GSM1060177</td><td>2.269510</td><td>2.602500</td><td>MIT</td></tr>
	<tr><td>GSM1060178</td><td>2.530203</td><td>2.930627</td><td>ATH</td></tr>
	<tr><td>GSM1060179</td><td>2.318810</td><td>2.558938</td><td>MIT</td></tr>
	<tr><td>GSM1060180</td><td>2.407185</td><td>2.729593</td><td>ATH</td></tr>
</tbody>
</table>




```R
p_brown  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       ggboxplot(x = 'group', y = 'Brown_Module', ylab = 'GSVA Score',title = 'Brown Module',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('MIT','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_brown
```


    
![png](Step5_Two_Group_Compare_GSE43292_files/Step5_Two_Group_Compare_GSE43292_37_0.png)
    



```R
library(patchwork)
```


```R
p0 <- p_yellow + p_brown + plot_layout(guides='collect') +  plot_annotation(title = 'GSE43292',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE43292_Two_Module.pdf', height = 6, width = 9, units = 'cm')
```


    
![png](Step5_Two_Group_Compare_GSE43292_files/Step5_Two_Group_Compare_GSE43292_39_0.png)
    

