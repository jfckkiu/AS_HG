```R
library(tidyverse)
library(GSVA)
library(GEOquery) 
library(limma)
library(ggpubr)
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
gset  <-  GEOquery::getGEO('GSE163154',getGPL = F)
```

    Found 1 file(s)
    
    GSE163154_series_matrix.txt.gz
    
    [1mRows: [22m[34m22184[39m [1mColumns: [22m[34m44[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m "\t"
    [31mchr[39m  (1): ID_REF
    [32mdbl[39m (43): GSM4973394, GSM4973395, GSM4973397, GSM4973398, GSM4973400, GSM497...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
normalized_gset  <- gset$GSE163154_series_matrix.txt.gz@assayData$exprs
```


```R
range(normalized_gset)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>7.097614903</li><li>15.45184144</li></ol>




```R
normalized_gset
```


<table class="dataframe">
<caption>A matrix: 22184 × 43 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>GSM4973394</th><th scope=col>GSM4973395</th><th scope=col>GSM4973397</th><th scope=col>GSM4973398</th><th scope=col>GSM4973400</th><th scope=col>GSM4973401</th><th scope=col>GSM4973403</th><th scope=col>GSM4973404</th><th scope=col>GSM4973406</th><th scope=col>GSM4973408</th><th scope=col>⋯</th><th scope=col>GSM4973443</th><th scope=col>GSM4973444</th><th scope=col>GSM4973445</th><th scope=col>GSM4973446</th><th scope=col>GSM4973448</th><th scope=col>GSM4973449</th><th scope=col>GSM4973451</th><th scope=col>GSM4973452</th><th scope=col>GSM4973454</th><th scope=col>GSM4973456</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ILMN_1343291</th><td>15.397685</td><td>15.397685</td><td>15.363611</td><td>15.397685</td><td>15.326456</td><td>15.363611</td><td>15.326456</td><td>15.397685</td><td>15.363611</td><td>15.253761</td><td>⋯</td><td>15.363611</td><td>15.451841</td><td>15.222846</td><td>15.451841</td><td>15.397685</td><td>15.451841</td><td>15.363611</td><td>15.363611</td><td>15.363611</td><td>15.451841</td></tr>
	<tr><th scope=row>ILMN_1343292</th><td>10.574274</td><td>11.283223</td><td>10.932375</td><td>10.952329</td><td>10.620145</td><td>10.958478</td><td>10.998223</td><td>11.137739</td><td>11.108668</td><td>10.572747</td><td>⋯</td><td>10.771249</td><td>10.451724</td><td>10.661556</td><td>10.443294</td><td>10.458157</td><td>10.990131</td><td>10.650094</td><td>10.908000</td><td>10.898592</td><td>11.075648</td></tr>
	<tr><th scope=row>ILMN_1343293</th><td>12.170338</td><td>12.394059</td><td>12.454637</td><td>11.881277</td><td>12.152035</td><td>12.514034</td><td>12.235775</td><td>12.034007</td><td>12.124991</td><td>11.795001</td><td>⋯</td><td>12.745573</td><td>12.599769</td><td>12.230211</td><td>12.443148</td><td>12.536527</td><td>12.916573</td><td>12.199100</td><td>12.471040</td><td>12.786596</td><td>11.730697</td></tr>
	<tr><th scope=row>ILMN_1343294</th><td>14.837210</td><td>14.979876</td><td>15.090183</td><td>15.029515</td><td>15.050317</td><td>15.192185</td><td>15.004360</td><td>15.070705</td><td>14.908203</td><td>15.029515</td><td>⋯</td><td>14.733469</td><td>15.253761</td><td>15.070705</td><td>15.029515</td><td>15.222846</td><td>15.070705</td><td>15.050317</td><td>14.908203</td><td>15.029515</td><td>15.137000</td></tr>
	<tr><th scope=row>ILMN_1651209</th><td> 7.327692</td><td> 7.297238</td><td> 7.368454</td><td> 7.440558</td><td> 7.374302</td><td> 7.407157</td><td> 7.238493</td><td> 7.511092</td><td> 7.356138</td><td> 7.465387</td><td>⋯</td><td> 7.300760</td><td> 7.425922</td><td> 7.321693</td><td> 7.347743</td><td> 7.411170</td><td> 7.320686</td><td> 7.290813</td><td> 7.320089</td><td> 7.324418</td><td> 7.396416</td></tr>
	<tr><th scope=row>ILMN_1651217</th><td> 7.339957</td><td> 7.451542</td><td> 7.425922</td><td> 7.520726</td><td> 7.608858</td><td> 7.471193</td><td> 7.506215</td><td> 7.344494</td><td> 7.435864</td><td> 7.448614</td><td>⋯</td><td> 7.954441</td><td> 8.146389</td><td> 7.818563</td><td> 7.636750</td><td> 7.790764</td><td> 7.688711</td><td> 7.683776</td><td> 7.784164</td><td> 7.926848</td><td> 7.472395</td></tr>
	<tr><th scope=row>ILMN_1651228</th><td>12.504059</td><td>12.815069</td><td>13.028573</td><td>12.203615</td><td>12.524005</td><td>13.020192</td><td>12.651050</td><td>12.124991</td><td>12.592494</td><td>12.855117</td><td>⋯</td><td>13.238425</td><td>12.588441</td><td>12.172591</td><td>12.067803</td><td>11.586562</td><td>12.062480</td><td>12.745573</td><td>12.309862</td><td>11.813200</td><td>12.865577</td></tr>
	<tr><th scope=row>ILMN_1651229</th><td> 8.705993</td><td> 8.932186</td><td> 8.936920</td><td> 8.612675</td><td> 8.825292</td><td> 8.742217</td><td> 8.725928</td><td> 8.816435</td><td> 8.298536</td><td> 8.550943</td><td>⋯</td><td> 8.913726</td><td> 8.502492</td><td> 8.827003</td><td> 8.690535</td><td> 8.877260</td><td> 8.878386</td><td> 8.392278</td><td> 8.888091</td><td> 8.941801</td><td> 9.224220</td></tr>
	<tr><th scope=row>ILMN_1651234</th><td> 7.205112</td><td> 7.207519</td><td> 7.208132</td><td> 7.231922</td><td> 7.196864</td><td> 7.211683</td><td> 7.209276</td><td> 7.167135</td><td> 7.218240</td><td> 7.203934</td><td>⋯</td><td> 7.204305</td><td> 7.225860</td><td> 7.202406</td><td> 7.187971</td><td> 7.212011</td><td> 7.197252</td><td> 7.183154</td><td> 7.197405</td><td> 7.229501</td><td> 7.184438</td></tr>
	<tr><th scope=row>ILMN_1651235</th><td> 7.322992</td><td> 7.302700</td><td> 7.278236</td><td> 7.293895</td><td> 7.268788</td><td> 7.304510</td><td> 7.233158</td><td> 7.294420</td><td> 7.264209</td><td> 7.321807</td><td>⋯</td><td> 7.298919</td><td> 7.271995</td><td> 7.299993</td><td> 7.270564</td><td> 7.274080</td><td> 7.263196</td><td> 7.271906</td><td> 7.302987</td><td> 7.300631</td><td> 7.263695</td></tr>
	<tr><th scope=row>ILMN_1651236</th><td> 7.281923</td><td> 7.291364</td><td> 7.291776</td><td> 7.248889</td><td> 7.266008</td><td> 7.279390</td><td> 7.232793</td><td> 7.303736</td><td> 7.245969</td><td> 7.299639</td><td>⋯</td><td> 7.293446</td><td> 7.260502</td><td> 7.260496</td><td> 7.230649</td><td> 7.242725</td><td> 7.292930</td><td> 7.258358</td><td> 7.247163</td><td> 7.244898</td><td> 7.239107</td></tr>
	<tr><th scope=row>ILMN_1651237</th><td> 7.379488</td><td> 7.415809</td><td> 7.196129</td><td> 7.314899</td><td> 7.330289</td><td> 7.352524</td><td> 7.275882</td><td> 7.227536</td><td> 7.317302</td><td> 7.208287</td><td>⋯</td><td> 7.600373</td><td> 7.433631</td><td> 7.450794</td><td> 7.525683</td><td> 7.544508</td><td> 7.509081</td><td> 7.495575</td><td> 7.504392</td><td> 7.477984</td><td> 7.389196</td></tr>
	<tr><th scope=row>ILMN_1651238</th><td> 7.226368</td><td> 7.270779</td><td> 7.271324</td><td> 7.265822</td><td> 7.279299</td><td> 7.244067</td><td> 7.315254</td><td> 7.278951</td><td> 7.241511</td><td> 7.293166</td><td>⋯</td><td> 7.265188</td><td> 7.247178</td><td> 7.259644</td><td> 7.289411</td><td> 7.295263</td><td> 7.249965</td><td> 7.219809</td><td> 7.252342</td><td> 7.269912</td><td> 7.272036</td></tr>
	<tr><th scope=row>ILMN_1651254</th><td>13.806967</td><td>13.539681</td><td>14.146502</td><td>13.637754</td><td>13.918210</td><td>13.521846</td><td>13.610869</td><td>14.481605</td><td>14.307499</td><td>13.918210</td><td>⋯</td><td>11.944286</td><td>12.678690</td><td>12.779663</td><td>12.172591</td><td>13.033257</td><td>11.715358</td><td>12.418204</td><td>12.783150</td><td>12.484264</td><td>13.688502</td></tr>
	<tr><th scope=row>ILMN_1651259</th><td> 7.423039</td><td> 7.579675</td><td> 7.657879</td><td> 7.532986</td><td> 8.030292</td><td> 7.551360</td><td> 8.311267</td><td> 7.772092</td><td> 7.826081</td><td> 8.020862</td><td>⋯</td><td> 7.448659</td><td> 7.894967</td><td> 7.657879</td><td> 7.604679</td><td> 7.702206</td><td> 7.485217</td><td> 7.746034</td><td> 7.736513</td><td> 7.756585</td><td> 7.410618</td></tr>
	<tr><th scope=row>ILMN_1651260</th><td> 7.281846</td><td> 7.195189</td><td> 7.249309</td><td> 7.215029</td><td> 7.248469</td><td> 7.279541</td><td> 7.308241</td><td> 7.182207</td><td> 7.261013</td><td> 7.220732</td><td>⋯</td><td> 7.228335</td><td> 7.236523</td><td> 7.244799</td><td> 7.256900</td><td> 7.248220</td><td> 7.268602</td><td> 7.255618</td><td> 7.234740</td><td> 7.240352</td><td> 7.246395</td></tr>
	<tr><th scope=row>ILMN_1651261</th><td> 8.518367</td><td> 8.294897</td><td> 8.586178</td><td> 8.653123</td><td> 9.279381</td><td> 8.484754</td><td> 9.339030</td><td> 8.959649</td><td> 8.568260</td><td> 9.326065</td><td>⋯</td><td> 8.356062</td><td> 9.077118</td><td> 8.724043</td><td> 8.892479</td><td> 8.450167</td><td> 8.755510</td><td> 8.493736</td><td> 8.726614</td><td> 8.444546</td><td> 7.786847</td></tr>
	<tr><th scope=row>ILMN_1651262</th><td>12.085026</td><td>12.206546</td><td>12.118636</td><td>12.054038</td><td>12.282975</td><td>11.912645</td><td>12.233056</td><td>12.031678</td><td>11.621361</td><td>12.090703</td><td>⋯</td><td>12.183902</td><td>12.192782</td><td>12.065365</td><td>12.105938</td><td>12.088591</td><td>12.133091</td><td>12.090703</td><td>12.422287</td><td>12.378062</td><td>12.113812</td></tr>
	<tr><th scope=row>ILMN_1651268</th><td> 7.414458</td><td> 7.443333</td><td> 7.589320</td><td> 7.374094</td><td> 7.333371</td><td> 7.383858</td><td> 7.413800</td><td> 7.543010</td><td> 7.479206</td><td> 7.392411</td><td>⋯</td><td> 7.267249</td><td> 7.287595</td><td> 7.338392</td><td> 7.321902</td><td> 7.399667</td><td> 7.373791</td><td> 7.321732</td><td> 7.343050</td><td> 7.305435</td><td> 7.319823</td></tr>
	<tr><th scope=row>ILMN_1651278</th><td> 9.165145</td><td> 9.052038</td><td> 9.275565</td><td> 8.994327</td><td> 8.987752</td><td> 9.078937</td><td> 9.119647</td><td> 9.197400</td><td> 9.178676</td><td> 9.020888</td><td>⋯</td><td> 9.015577</td><td> 8.750802</td><td> 8.797608</td><td> 8.925028</td><td> 9.008742</td><td> 8.925424</td><td> 8.939920</td><td> 8.765578</td><td> 8.933972</td><td> 8.984666</td></tr>
	<tr><th scope=row>ILMN_1651282</th><td> 7.345924</td><td> 7.591585</td><td> 7.311619</td><td> 7.339182</td><td> 7.377170</td><td> 7.266417</td><td> 7.295035</td><td> 7.299882</td><td> 7.336421</td><td> 7.534868</td><td>⋯</td><td> 7.486561</td><td> 7.580888</td><td> 7.403186</td><td> 7.884225</td><td> 7.573875</td><td> 7.523448</td><td> 7.807115</td><td> 7.531030</td><td> 7.699306</td><td> 7.511566</td></tr>
	<tr><th scope=row>ILMN_1651286</th><td> 7.319776</td><td> 7.319418</td><td> 7.299010</td><td> 7.316436</td><td> 7.349118</td><td> 7.365936</td><td> 7.272071</td><td> 7.389325</td><td> 7.369159</td><td> 7.376327</td><td>⋯</td><td> 7.354768</td><td> 7.314297</td><td> 7.372413</td><td> 7.376327</td><td> 7.284528</td><td> 7.348340</td><td> 7.301692</td><td> 7.360388</td><td> 7.335202</td><td> 7.303337</td></tr>
	<tr><th scope=row>ILMN_1651296</th><td> 7.725394</td><td> 7.799781</td><td> 7.715952</td><td> 7.901950</td><td> 7.691666</td><td> 7.793375</td><td> 7.676954</td><td> 7.798431</td><td> 7.653237</td><td> 7.854445</td><td>⋯</td><td> 7.816222</td><td> 7.763374</td><td> 7.956398</td><td> 7.820769</td><td> 7.883409</td><td> 7.998838</td><td> 7.966458</td><td> 7.715893</td><td> 7.916947</td><td> 7.830545</td></tr>
	<tr><th scope=row>ILMN_1651298</th><td> 7.253394</td><td> 7.246765</td><td> 7.278693</td><td> 7.264011</td><td> 7.238399</td><td> 7.256813</td><td> 7.210356</td><td> 7.287905</td><td> 7.244855</td><td> 7.279576</td><td>⋯</td><td> 7.235998</td><td> 7.270805</td><td> 7.273290</td><td> 7.242019</td><td> 7.272134</td><td> 7.281915</td><td> 7.269409</td><td> 7.235215</td><td> 7.290236</td><td> 7.259300</td></tr>
	<tr><th scope=row>ILMN_1651303</th><td> 7.222736</td><td> 7.247647</td><td> 7.232780</td><td> 7.297565</td><td> 7.254969</td><td> 7.236302</td><td> 7.195875</td><td> 7.268146</td><td> 7.296475</td><td> 7.254227</td><td>⋯</td><td> 7.229281</td><td> 7.252944</td><td> 7.236959</td><td> 7.238020</td><td> 7.237577</td><td> 7.268429</td><td> 7.219242</td><td> 7.235599</td><td> 7.223959</td><td> 7.220732</td></tr>
	<tr><th scope=row>ILMN_1651316</th><td> 7.717240</td><td> 7.575777</td><td> 7.351520</td><td> 7.757911</td><td> 7.424758</td><td> 7.413289</td><td> 7.325682</td><td> 7.336730</td><td> 7.603029</td><td> 7.594832</td><td>⋯</td><td> 7.877162</td><td> 8.065558</td><td> 7.878759</td><td> 7.845327</td><td> 7.502832</td><td> 7.444193</td><td> 7.576983</td><td> 7.725902</td><td> 8.354902</td><td> 7.665653</td></tr>
	<tr><th scope=row>ILMN_1651330</th><td> 7.328595</td><td> 7.272354</td><td> 7.334431</td><td> 7.290530</td><td> 7.298802</td><td> 7.353108</td><td> 7.308056</td><td> 7.301047</td><td> 7.324117</td><td> 7.246123</td><td>⋯</td><td> 7.252330</td><td> 7.310028</td><td> 7.305109</td><td> 7.287330</td><td> 7.247985</td><td> 7.284730</td><td> 7.311427</td><td> 7.270429</td><td> 7.310882</td><td> 7.295372</td></tr>
	<tr><th scope=row>ILMN_1651336</th><td> 7.797771</td><td> 7.661095</td><td> 8.029717</td><td> 7.694936</td><td> 7.826081</td><td> 7.761593</td><td> 7.757811</td><td> 7.846034</td><td> 8.326006</td><td> 7.943099</td><td>⋯</td><td> 7.836122</td><td> 7.669512</td><td> 7.634384</td><td> 7.763549</td><td> 7.736182</td><td> 7.793086</td><td> 7.726998</td><td> 7.775012</td><td> 7.690243</td><td> 7.763138</td></tr>
	<tr><th scope=row>ILMN_1651339</th><td> 7.324689</td><td> 7.318143</td><td> 7.352173</td><td> 7.345603</td><td> 7.352943</td><td> 7.300459</td><td> 7.416398</td><td> 7.292177</td><td> 7.318286</td><td> 7.319646</td><td>⋯</td><td> 7.398716</td><td> 7.327692</td><td> 7.325016</td><td> 7.303049</td><td> 7.334726</td><td> 7.340936</td><td> 7.249891</td><td> 7.339242</td><td> 7.335891</td><td> 7.293805</td></tr>
	<tr><th scope=row>ILMN_1651343</th><td>11.046277</td><td>10.596099</td><td>11.763076</td><td>10.866014</td><td>10.504696</td><td>10.532986</td><td>10.978685</td><td>10.949018</td><td>10.661556</td><td>10.312945</td><td>⋯</td><td> 9.338134</td><td> 8.696550</td><td> 9.684871</td><td> 9.175033</td><td> 9.533054</td><td> 9.130688</td><td> 9.937437</td><td>10.008520</td><td> 9.170616</td><td>10.219657</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>ILMN_1815719</th><td> 9.179718</td><td> 9.587973</td><td> 7.750780</td><td> 9.546750</td><td> 8.997113</td><td> 9.325212</td><td> 8.590985</td><td> 8.110609</td><td> 8.772737</td><td> 9.657057</td><td>⋯</td><td>10.120668</td><td>10.245688</td><td> 9.796820</td><td> 9.874710</td><td> 9.903536</td><td>10.358298</td><td>10.256998</td><td>10.033637</td><td>10.388640</td><td>10.324245</td></tr>
	<tr><th scope=row>ILMN_1815723</th><td> 7.784783</td><td> 7.728034</td><td> 8.009661</td><td> 7.861535</td><td> 7.993029</td><td> 7.736773</td><td> 8.340097</td><td> 7.973944</td><td> 8.024763</td><td> 8.129636</td><td>⋯</td><td> 7.769880</td><td> 8.002647</td><td> 8.033048</td><td> 7.884821</td><td> 7.940472</td><td> 7.899521</td><td> 7.937791</td><td> 8.100283</td><td> 7.937791</td><td> 7.637534</td></tr>
	<tr><th scope=row>ILMN_1815733</th><td> 9.800157</td><td> 9.857971</td><td> 9.484091</td><td> 9.775808</td><td> 9.554961</td><td>10.061231</td><td> 9.620250</td><td> 9.483040</td><td> 9.488859</td><td> 9.593567</td><td>⋯</td><td> 9.754671</td><td> 9.856830</td><td> 9.083321</td><td> 9.065944</td><td> 9.789384</td><td> 9.850954</td><td> 9.716034</td><td> 9.601086</td><td> 9.643289</td><td> 9.600608</td></tr>
	<tr><th scope=row>ILMN_1815734</th><td>10.573525</td><td>10.239281</td><td>10.616719</td><td>10.658819</td><td>10.372935</td><td>10.520772</td><td>10.282220</td><td>10.683957</td><td>10.813227</td><td>10.429463</td><td>⋯</td><td>10.055726</td><td>10.098576</td><td>10.226348</td><td>10.162095</td><td>10.264295</td><td> 9.942608</td><td>10.310882</td><td>10.321336</td><td> 9.964624</td><td> 9.938990</td></tr>
	<tr><th scope=row>ILMN_1815745</th><td> 8.935056</td><td> 8.858813</td><td> 9.254408</td><td> 8.488450</td><td> 8.774254</td><td> 8.824948</td><td> 8.681662</td><td> 8.998565</td><td> 8.924661</td><td> 8.906758</td><td>⋯</td><td> 8.639559</td><td> 8.787256</td><td> 8.517182</td><td> 8.810712</td><td> 8.820930</td><td> 8.748685</td><td> 8.952298</td><td> 8.816727</td><td> 8.602641</td><td> 8.405424</td></tr>
	<tr><th scope=row>ILMN_1815759</th><td> 7.896605</td><td> 7.987875</td><td> 7.962446</td><td> 7.915164</td><td> 7.914356</td><td> 7.895223</td><td> 8.062391</td><td> 7.952463</td><td> 7.807262</td><td> 7.821944</td><td>⋯</td><td> 8.130087</td><td> 8.104837</td><td> 8.297802</td><td> 8.158375</td><td> 8.149317</td><td> 8.132279</td><td> 8.000541</td><td> 8.135865</td><td> 8.013749</td><td> 7.969219</td></tr>
	<tr><th scope=row>ILMN_1815773</th><td> 7.262818</td><td> 7.256086</td><td> 7.271000</td><td> 7.247518</td><td> 7.249999</td><td> 7.243351</td><td> 7.304826</td><td> 7.229646</td><td> 7.269641</td><td> 7.244386</td><td>⋯</td><td> 7.249389</td><td> 7.254027</td><td> 7.247485</td><td> 7.251096</td><td> 7.227695</td><td> 7.221567</td><td> 7.269739</td><td> 7.272506</td><td> 7.223959</td><td> 7.254033</td></tr>
	<tr><th scope=row>ILMN_1815780</th><td> 7.204965</td><td> 7.221138</td><td> 7.193708</td><td> 7.217951</td><td> 7.400085</td><td> 7.214288</td><td> 7.440597</td><td> 7.244921</td><td> 7.287981</td><td> 7.271233</td><td>⋯</td><td> 7.260283</td><td> 7.325452</td><td> 7.260323</td><td> 7.204522</td><td> 7.243569</td><td> 7.203947</td><td> 7.216379</td><td> 7.286438</td><td> 7.246184</td><td> 7.184053</td></tr>
	<tr><th scope=row>ILMN_1815795</th><td> 8.713454</td><td> 9.229663</td><td> 8.590365</td><td> 8.660969</td><td> 8.850444</td><td> 8.695236</td><td> 8.770857</td><td> 8.659100</td><td> 8.557428</td><td> 8.470516</td><td>⋯</td><td> 9.013579</td><td> 8.980781</td><td> 9.246998</td><td> 9.188676</td><td> 9.001238</td><td> 9.310173</td><td> 8.996418</td><td> 9.151311</td><td> 9.010747</td><td> 9.021468</td></tr>
	<tr><th scope=row>ILMN_1815796</th><td> 7.261643</td><td> 7.228978</td><td> 7.260609</td><td> 7.249508</td><td> 7.242441</td><td> 7.247732</td><td> 7.169445</td><td> 7.269726</td><td> 7.246200</td><td> 7.251775</td><td>⋯</td><td> 7.258633</td><td> 7.304888</td><td> 7.264362</td><td> 7.234619</td><td> 7.290411</td><td> 7.272024</td><td> 7.227777</td><td> 7.257393</td><td> 7.236692</td><td> 7.285235</td></tr>
	<tr><th scope=row>ILMN_1815812</th><td> 8.092522</td><td> 7.951864</td><td> 8.576777</td><td> 8.017761</td><td> 8.077110</td><td> 7.816763</td><td> 8.297483</td><td> 8.539556</td><td> 8.264162</td><td> 8.055589</td><td>⋯</td><td> 7.653366</td><td> 7.726238</td><td> 7.799419</td><td> 7.764777</td><td> 7.810924</td><td> 7.735018</td><td> 7.722065</td><td> 7.754737</td><td> 7.820164</td><td> 7.894886</td></tr>
	<tr><th scope=row>ILMN_1815840</th><td> 7.340552</td><td> 7.334196</td><td> 7.405242</td><td> 7.353365</td><td> 7.328280</td><td> 7.359995</td><td> 7.333042</td><td> 7.365936</td><td> 7.334075</td><td> 7.320950</td><td>⋯</td><td> 7.385634</td><td> 7.299850</td><td> 7.302125</td><td> 7.327393</td><td> 7.323107</td><td> 7.414249</td><td> 7.314594</td><td> 7.332471</td><td> 7.369037</td><td> 7.463679</td></tr>
	<tr><th scope=row>ILMN_1815849</th><td> 7.242466</td><td> 7.177864</td><td> 7.192287</td><td> 7.177105</td><td> 7.154340</td><td> 7.172473</td><td> 7.220077</td><td> 7.180943</td><td> 7.178536</td><td> 7.176244</td><td>⋯</td><td> 7.213631</td><td> 7.166032</td><td> 7.164261</td><td> 7.201924</td><td> 7.165135</td><td> 7.209092</td><td> 7.177843</td><td> 7.171619</td><td> 7.184355</td><td> 7.183927</td></tr>
	<tr><th scope=row>ILMN_1815859</th><td> 7.961346</td><td> 8.055121</td><td> 7.959272</td><td> 8.099019</td><td> 7.954359</td><td> 8.042481</td><td> 8.085968</td><td> 8.201399</td><td> 8.159219</td><td> 7.794292</td><td>⋯</td><td> 7.934172</td><td> 8.104294</td><td> 8.053927</td><td> 8.177804</td><td> 8.154593</td><td> 8.180488</td><td> 7.990242</td><td> 8.015188</td><td> 7.886453</td><td> 8.367435</td></tr>
	<tr><th scope=row>ILMN_1815874</th><td> 8.997746</td><td> 9.171955</td><td> 8.675459</td><td> 8.833961</td><td> 9.099671</td><td> 8.988140</td><td> 8.959330</td><td> 8.706638</td><td> 8.940344</td><td> 8.673141</td><td>⋯</td><td> 9.476413</td><td> 9.136816</td><td> 9.382961</td><td> 9.232052</td><td> 9.201560</td><td> 9.234646</td><td> 9.266216</td><td> 9.291357</td><td> 9.166262</td><td> 9.272135</td></tr>
	<tr><th scope=row>ILMN_1815878</th><td>10.589751</td><td>10.736022</td><td>10.583114</td><td>10.369826</td><td>10.556491</td><td>10.951441</td><td>10.559109</td><td>10.240076</td><td>10.038857</td><td>10.395699</td><td>⋯</td><td>10.783207</td><td>11.160426</td><td>10.829010</td><td>11.142211</td><td>10.993595</td><td>11.284070</td><td>11.075648</td><td>10.808434</td><td>10.636503</td><td>10.896888</td></tr>
	<tr><th scope=row>ILMN_1815882</th><td> 8.077857</td><td> 7.917237</td><td> 8.374968</td><td> 8.135050</td><td> 8.446158</td><td> 8.062890</td><td> 8.616791</td><td> 8.418736</td><td> 8.396691</td><td> 8.292994</td><td>⋯</td><td> 7.756433</td><td> 8.332349</td><td> 8.053686</td><td> 8.094642</td><td> 8.094066</td><td> 7.826081</td><td> 8.230748</td><td> 8.045888</td><td> 7.915406</td><td> 7.711444</td></tr>
	<tr><th scope=row>ILMN_1815885</th><td> 8.099842</td><td> 8.084597</td><td> 8.281198</td><td> 8.093233</td><td> 8.040160</td><td> 8.099408</td><td> 8.124230</td><td> 8.349547</td><td> 8.488450</td><td> 8.133877</td><td>⋯</td><td> 7.835938</td><td> 7.947613</td><td> 8.007087</td><td> 7.964730</td><td> 8.046987</td><td> 8.038710</td><td> 8.193979</td><td> 8.025866</td><td> 7.885480</td><td> 8.142848</td></tr>
	<tr><th scope=row>ILMN_1815890</th><td> 7.476659</td><td> 7.580629</td><td> 7.214657</td><td> 7.563308</td><td> 7.463240</td><td> 7.428503</td><td> 7.313171</td><td> 7.311707</td><td> 7.460932</td><td> 7.344041</td><td>⋯</td><td> 7.781759</td><td> 7.789116</td><td> 7.641788</td><td> 7.707770</td><td> 7.527691</td><td> 7.596968</td><td> 7.615278</td><td> 7.591839</td><td> 7.538505</td><td> 7.653301</td></tr>
	<tr><th scope=row>ILMN_1815908</th><td> 8.210155</td><td> 7.832124</td><td> 8.396115</td><td> 8.134584</td><td> 8.022541</td><td> 7.902125</td><td> 8.321463</td><td> 8.222483</td><td> 7.936109</td><td> 8.024493</td><td>⋯</td><td> 7.703721</td><td> 7.824096</td><td> 7.751923</td><td> 7.582565</td><td> 7.599136</td><td> 7.571474</td><td> 8.009066</td><td> 7.879705</td><td> 7.596559</td><td> 7.880327</td></tr>
	<tr><th scope=row>ILMN_1815923</th><td> 7.170334</td><td> 7.215901</td><td> 7.202828</td><td> 7.232299</td><td> 7.237880</td><td> 7.244166</td><td> 7.246504</td><td> 7.296397</td><td> 7.233887</td><td> 7.186446</td><td>⋯</td><td> 7.239220</td><td> 7.278047</td><td> 7.168953</td><td> 7.241946</td><td> 7.231832</td><td> 7.182938</td><td> 7.206737</td><td> 7.215314</td><td> 7.214067</td><td> 7.258633</td></tr>
	<tr><th scope=row>ILMN_1815924</th><td> 9.365133</td><td> 9.353383</td><td> 9.375704</td><td> 9.162812</td><td> 9.223520</td><td> 9.206855</td><td> 9.251731</td><td> 9.358738</td><td> 9.790366</td><td> 9.265230</td><td>⋯</td><td> 9.118502</td><td> 9.203148</td><td> 9.157930</td><td> 9.289189</td><td> 9.343734</td><td> 9.388720</td><td> 9.297218</td><td> 9.462036</td><td> 9.535842</td><td> 9.318426</td></tr>
	<tr><th scope=row>ILMN_1815933</th><td> 8.087673</td><td> 7.976249</td><td> 7.943099</td><td> 8.108333</td><td> 7.924111</td><td> 8.065171</td><td> 7.832018</td><td> 8.005984</td><td> 8.163325</td><td> 7.796148</td><td>⋯</td><td> 8.233175</td><td> 8.094740</td><td> 8.142293</td><td> 8.158009</td><td> 8.279694</td><td> 8.286033</td><td> 8.022717</td><td> 8.134822</td><td> 8.025610</td><td> 8.474256</td></tr>
	<tr><th scope=row>ILMN_1815937</th><td> 7.252676</td><td> 7.293665</td><td> 7.371083</td><td> 7.311351</td><td> 7.282919</td><td> 7.295841</td><td> 7.259993</td><td> 7.285671</td><td> 7.315877</td><td> 7.275810</td><td>⋯</td><td> 7.284502</td><td> 7.257732</td><td> 7.296783</td><td> 7.303240</td><td> 7.310290</td><td> 7.265408</td><td> 7.283729</td><td> 7.295294</td><td> 7.315117</td><td> 7.275810</td></tr>
	<tr><th scope=row>ILMN_1815938</th><td> 7.230973</td><td> 7.251768</td><td> 7.250046</td><td> 7.253949</td><td> 7.210804</td><td> 7.240074</td><td> 7.236095</td><td> 7.252055</td><td> 7.213958</td><td> 7.278176</td><td>⋯</td><td> 7.223423</td><td> 7.238224</td><td> 7.244267</td><td> 7.272532</td><td> 7.244075</td><td> 7.215729</td><td> 7.212128</td><td> 7.207027</td><td> 7.223263</td><td> 7.258715</td></tr>
	<tr><th scope=row>ILMN_1815941</th><td> 9.798893</td><td> 9.694271</td><td> 9.907072</td><td> 9.768732</td><td> 9.895686</td><td> 9.764151</td><td>10.004447</td><td>10.544573</td><td> 9.749108</td><td> 9.595786</td><td>⋯</td><td> 9.256205</td><td> 9.098114</td><td> 9.498280</td><td> 9.132293</td><td> 9.248486</td><td> 8.758983</td><td> 9.800157</td><td> 9.500341</td><td> 9.520802</td><td> 9.415060</td></tr>
	<tr><th scope=row>ILMN_1815951</th><td> 7.966225</td><td> 7.953338</td><td> 7.804860</td><td> 7.807950</td><td> 7.807399</td><td> 7.812572</td><td> 7.878314</td><td> 7.683947</td><td> 7.699060</td><td> 7.849360</td><td>⋯</td><td> 7.907033</td><td> 7.762312</td><td> 7.762570</td><td> 7.698209</td><td> 7.607156</td><td> 7.593200</td><td> 7.693124</td><td> 7.747622</td><td> 7.677472</td><td> 7.663077</td></tr>
	<tr><th scope=row>ILMN_2038774</th><td>14.928846</td><td>15.050317</td><td>15.165482</td><td>15.137000</td><td>15.363611</td><td>14.908203</td><td>15.397685</td><td>15.192185</td><td>15.117043</td><td>15.290683</td><td>⋯</td><td>15.326456</td><td>15.290683</td><td>15.363611</td><td>15.363611</td><td>15.137000</td><td>15.253761</td><td>15.326456</td><td>15.397685</td><td>15.326456</td><td>14.885728</td></tr>
	<tr><th scope=row>ILMN_2038777</th><td>15.004360</td><td>15.290683</td><td>14.758077</td><td>15.165482</td><td>15.192185</td><td>14.733469</td><td>15.137000</td><td>15.029515</td><td>14.664041</td><td>15.070705</td><td>⋯</td><td>14.908203</td><td>15.029515</td><td>15.253761</td><td>14.837210</td><td>15.029515</td><td>14.634348</td><td>14.908203</td><td>15.290683</td><td>15.290683</td><td>14.712066</td></tr>
	<tr><th scope=row>ILMN_2038778</th><td>14.077380</td><td>13.952513</td><td>13.661227</td><td>13.999062</td><td>13.964369</td><td>14.229570</td><td>14.098730</td><td>13.888182</td><td>13.194778</td><td>13.821948</td><td>⋯</td><td>14.134699</td><td>13.850222</td><td>14.684127</td><td>14.216344</td><td>14.206824</td><td>14.335820</td><td>14.066954</td><td>14.265875</td><td>14.010573</td><td>13.734694</td></tr>
</tbody>
</table>




```R
ids <- AnnoProbe::idmap('GPL6104',type = 'soft')
```

    file downloaded in /home/jfckkiu/AS_HG/Final_Results/Figure2
    



```R
ids   <- ids  %>% 
           mutate(ID = as.character(ID))
```


```R
normalized_gset %>% 
 as.data.frame() %>% 
  rownames_to_column(var = 'Probe') %>%  
   left_join(ids, by = c('Probe'= 'ID'))
```


<table class="dataframe">
<caption>A data.frame: 22184 × 45</caption>
<thead>
	<tr><th scope=col>Probe</th><th scope=col>GSM4973394</th><th scope=col>GSM4973395</th><th scope=col>GSM4973397</th><th scope=col>GSM4973398</th><th scope=col>GSM4973400</th><th scope=col>GSM4973401</th><th scope=col>GSM4973403</th><th scope=col>GSM4973404</th><th scope=col>GSM4973406</th><th scope=col>⋯</th><th scope=col>GSM4973444</th><th scope=col>GSM4973445</th><th scope=col>GSM4973446</th><th scope=col>GSM4973448</th><th scope=col>GSM4973449</th><th scope=col>GSM4973451</th><th scope=col>GSM4973452</th><th scope=col>GSM4973454</th><th scope=col>GSM4973456</th><th scope=col>symbol</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>ILMN_1343291</td><td>15.397685</td><td>15.397685</td><td>15.363611</td><td>15.397685</td><td>15.326456</td><td>15.363611</td><td>15.326456</td><td>15.397685</td><td>15.363611</td><td>⋯</td><td>15.451841</td><td>15.222846</td><td>15.451841</td><td>15.397685</td><td>15.451841</td><td>15.363611</td><td>15.363611</td><td>15.363611</td><td>15.451841</td><td>EEF1A1   </td></tr>
	<tr><td>ILMN_1343292</td><td>10.574274</td><td>11.283223</td><td>10.932375</td><td>10.952329</td><td>10.620145</td><td>10.958478</td><td>10.998223</td><td>11.137739</td><td>11.108668</td><td>⋯</td><td>10.451724</td><td>10.661556</td><td>10.443294</td><td>10.458157</td><td>10.990131</td><td>10.650094</td><td>10.908000</td><td>10.898592</td><td>11.075648</td><td>TUBB     </td></tr>
	<tr><td>ILMN_1343293</td><td>12.170338</td><td>12.394059</td><td>12.454637</td><td>11.881277</td><td>12.152035</td><td>12.514034</td><td>12.235775</td><td>12.034007</td><td>12.124991</td><td>⋯</td><td>12.599769</td><td>12.230211</td><td>12.443148</td><td>12.536527</td><td>12.916573</td><td>12.199100</td><td>12.471040</td><td>12.786596</td><td>11.730697</td><td>TXN      </td></tr>
	<tr><td>ILMN_1343294</td><td>14.837210</td><td>14.979876</td><td>15.090183</td><td>15.029515</td><td>15.050317</td><td>15.192185</td><td>15.004360</td><td>15.070705</td><td>14.908203</td><td>⋯</td><td>15.253761</td><td>15.070705</td><td>15.029515</td><td>15.222846</td><td>15.070705</td><td>15.050317</td><td>14.908203</td><td>15.029515</td><td>15.137000</td><td>ACTB     </td></tr>
	<tr><td>ILMN_1651209</td><td> 7.327692</td><td> 7.297238</td><td> 7.368454</td><td> 7.440558</td><td> 7.374302</td><td> 7.407157</td><td> 7.238493</td><td> 7.511092</td><td> 7.356138</td><td>⋯</td><td> 7.425922</td><td> 7.321693</td><td> 7.347743</td><td> 7.411170</td><td> 7.320686</td><td> 7.290813</td><td> 7.320089</td><td> 7.324418</td><td> 7.396416</td><td>SLC35E2  </td></tr>
	<tr><td>ILMN_1651217</td><td> 7.339957</td><td> 7.451542</td><td> 7.425922</td><td> 7.520726</td><td> 7.608858</td><td> 7.471193</td><td> 7.506215</td><td> 7.344494</td><td> 7.435864</td><td>⋯</td><td> 8.146389</td><td> 7.818563</td><td> 7.636750</td><td> 7.790764</td><td> 7.688711</td><td> 7.683776</td><td> 7.784164</td><td> 7.926848</td><td> 7.472395</td><td>PDCD1LG2 </td></tr>
	<tr><td>ILMN_1651228</td><td>12.504059</td><td>12.815069</td><td>13.028573</td><td>12.203615</td><td>12.524005</td><td>13.020192</td><td>12.651050</td><td>12.124991</td><td>12.592494</td><td>⋯</td><td>12.588441</td><td>12.172591</td><td>12.067803</td><td>11.586562</td><td>12.062480</td><td>12.745573</td><td>12.309862</td><td>11.813200</td><td>12.865577</td><td>RPS28    </td></tr>
	<tr><td>ILMN_1651229</td><td> 8.705993</td><td> 8.932186</td><td> 8.936920</td><td> 8.612675</td><td> 8.825292</td><td> 8.742217</td><td> 8.725928</td><td> 8.816435</td><td> 8.298536</td><td>⋯</td><td> 8.502492</td><td> 8.827003</td><td> 8.690535</td><td> 8.877260</td><td> 8.878386</td><td> 8.392278</td><td> 8.888091</td><td> 8.941801</td><td> 9.224220</td><td>IPO13    </td></tr>
	<tr><td>ILMN_1651234</td><td> 7.205112</td><td> 7.207519</td><td> 7.208132</td><td> 7.231922</td><td> 7.196864</td><td> 7.211683</td><td> 7.209276</td><td> 7.167135</td><td> 7.218240</td><td>⋯</td><td> 7.225860</td><td> 7.202406</td><td> 7.187971</td><td> 7.212011</td><td> 7.197252</td><td> 7.183154</td><td> 7.197405</td><td> 7.229501</td><td> 7.184438</td><td>SYT14    </td></tr>
	<tr><td>ILMN_1651235</td><td> 7.322992</td><td> 7.302700</td><td> 7.278236</td><td> 7.293895</td><td> 7.268788</td><td> 7.304510</td><td> 7.233158</td><td> 7.294420</td><td> 7.264209</td><td>⋯</td><td> 7.271995</td><td> 7.299993</td><td> 7.270564</td><td> 7.274080</td><td> 7.263196</td><td> 7.271906</td><td> 7.302987</td><td> 7.300631</td><td> 7.263695</td><td>AFAP1    </td></tr>
	<tr><td>ILMN_1651236</td><td> 7.281923</td><td> 7.291364</td><td> 7.291776</td><td> 7.248889</td><td> 7.266008</td><td> 7.279390</td><td> 7.232793</td><td> 7.303736</td><td> 7.245969</td><td>⋯</td><td> 7.260502</td><td> 7.260496</td><td> 7.230649</td><td> 7.242725</td><td> 7.292930</td><td> 7.258358</td><td> 7.247163</td><td> 7.244898</td><td> 7.239107</td><td>GGTLA4   </td></tr>
	<tr><td>ILMN_1651237</td><td> 7.379488</td><td> 7.415809</td><td> 7.196129</td><td> 7.314899</td><td> 7.330289</td><td> 7.352524</td><td> 7.275882</td><td> 7.227536</td><td> 7.317302</td><td>⋯</td><td> 7.433631</td><td> 7.450794</td><td> 7.525683</td><td> 7.544508</td><td> 7.509081</td><td> 7.495575</td><td> 7.504392</td><td> 7.477984</td><td> 7.389196</td><td>CDT1     </td></tr>
	<tr><td>ILMN_1651238</td><td> 7.226368</td><td> 7.270779</td><td> 7.271324</td><td> 7.265822</td><td> 7.279299</td><td> 7.244067</td><td> 7.315254</td><td> 7.278951</td><td> 7.241511</td><td>⋯</td><td> 7.247178</td><td> 7.259644</td><td> 7.289411</td><td> 7.295263</td><td> 7.249965</td><td> 7.219809</td><td> 7.252342</td><td> 7.269912</td><td> 7.272036</td><td>TRPV1    </td></tr>
	<tr><td>ILMN_1651254</td><td>13.806967</td><td>13.539681</td><td>14.146502</td><td>13.637754</td><td>13.918210</td><td>13.521846</td><td>13.610869</td><td>14.481605</td><td>14.307499</td><td>⋯</td><td>12.678690</td><td>12.779663</td><td>12.172591</td><td>13.033257</td><td>11.715358</td><td>12.418204</td><td>12.783150</td><td>12.484264</td><td>13.688502</td><td>LPP      </td></tr>
	<tr><td>ILMN_1651259</td><td> 7.423039</td><td> 7.579675</td><td> 7.657879</td><td> 7.532986</td><td> 8.030292</td><td> 7.551360</td><td> 8.311267</td><td> 7.772092</td><td> 7.826081</td><td>⋯</td><td> 7.894967</td><td> 7.657879</td><td> 7.604679</td><td> 7.702206</td><td> 7.485217</td><td> 7.746034</td><td> 7.736513</td><td> 7.756585</td><td> 7.410618</td><td>UGP2     </td></tr>
	<tr><td>ILMN_1651260</td><td> 7.281846</td><td> 7.195189</td><td> 7.249309</td><td> 7.215029</td><td> 7.248469</td><td> 7.279541</td><td> 7.308241</td><td> 7.182207</td><td> 7.261013</td><td>⋯</td><td> 7.236523</td><td> 7.244799</td><td> 7.256900</td><td> 7.248220</td><td> 7.268602</td><td> 7.255618</td><td> 7.234740</td><td> 7.240352</td><td> 7.246395</td><td>CCNE2    </td></tr>
	<tr><td>ILMN_1651261</td><td> 8.518367</td><td> 8.294897</td><td> 8.586178</td><td> 8.653123</td><td> 9.279381</td><td> 8.484754</td><td> 9.339030</td><td> 8.959649</td><td> 8.568260</td><td>⋯</td><td> 9.077118</td><td> 8.724043</td><td> 8.892479</td><td> 8.450167</td><td> 8.755510</td><td> 8.493736</td><td> 8.726614</td><td> 8.444546</td><td> 7.786847</td><td>RSU1     </td></tr>
	<tr><td>ILMN_1651262</td><td>12.085026</td><td>12.206546</td><td>12.118636</td><td>12.054038</td><td>12.282975</td><td>11.912645</td><td>12.233056</td><td>12.031678</td><td>11.621361</td><td>⋯</td><td>12.192782</td><td>12.065365</td><td>12.105938</td><td>12.088591</td><td>12.133091</td><td>12.090703</td><td>12.422287</td><td>12.378062</td><td>12.113812</td><td>HNRPAB   </td></tr>
	<tr><td>ILMN_1651268</td><td> 7.414458</td><td> 7.443333</td><td> 7.589320</td><td> 7.374094</td><td> 7.333371</td><td> 7.383858</td><td> 7.413800</td><td> 7.543010</td><td> 7.479206</td><td>⋯</td><td> 7.287595</td><td> 7.338392</td><td> 7.321902</td><td> 7.399667</td><td> 7.373791</td><td> 7.321732</td><td> 7.343050</td><td> 7.305435</td><td> 7.319823</td><td>LOH12CR1 </td></tr>
	<tr><td>ILMN_1651278</td><td> 9.165145</td><td> 9.052038</td><td> 9.275565</td><td> 8.994327</td><td> 8.987752</td><td> 9.078937</td><td> 9.119647</td><td> 9.197400</td><td> 9.178676</td><td>⋯</td><td> 8.750802</td><td> 8.797608</td><td> 8.925028</td><td> 9.008742</td><td> 8.925424</td><td> 8.939920</td><td> 8.765578</td><td> 8.933972</td><td> 8.984666</td><td>SNIP1    </td></tr>
	<tr><td>ILMN_1651282</td><td> 7.345924</td><td> 7.591585</td><td> 7.311619</td><td> 7.339182</td><td> 7.377170</td><td> 7.266417</td><td> 7.295035</td><td> 7.299882</td><td> 7.336421</td><td>⋯</td><td> 7.580888</td><td> 7.403186</td><td> 7.884225</td><td> 7.573875</td><td> 7.523448</td><td> 7.807115</td><td> 7.531030</td><td> 7.699306</td><td> 7.511566</td><td>COL17A1  </td></tr>
	<tr><td>ILMN_1651286</td><td> 7.319776</td><td> 7.319418</td><td> 7.299010</td><td> 7.316436</td><td> 7.349118</td><td> 7.365936</td><td> 7.272071</td><td> 7.389325</td><td> 7.369159</td><td>⋯</td><td> 7.314297</td><td> 7.372413</td><td> 7.376327</td><td> 7.284528</td><td> 7.348340</td><td> 7.301692</td><td> 7.360388</td><td> 7.335202</td><td> 7.303337</td><td>GRHL1    </td></tr>
	<tr><td>ILMN_1651296</td><td> 7.725394</td><td> 7.799781</td><td> 7.715952</td><td> 7.901950</td><td> 7.691666</td><td> 7.793375</td><td> 7.676954</td><td> 7.798431</td><td> 7.653237</td><td>⋯</td><td> 7.763374</td><td> 7.956398</td><td> 7.820769</td><td> 7.883409</td><td> 7.998838</td><td> 7.966458</td><td> 7.715893</td><td> 7.916947</td><td> 7.830545</td><td>LOC143666</td></tr>
	<tr><td>ILMN_1651298</td><td> 7.253394</td><td> 7.246765</td><td> 7.278693</td><td> 7.264011</td><td> 7.238399</td><td> 7.256813</td><td> 7.210356</td><td> 7.287905</td><td> 7.244855</td><td>⋯</td><td> 7.270805</td><td> 7.273290</td><td> 7.242019</td><td> 7.272134</td><td> 7.281915</td><td> 7.269409</td><td> 7.235215</td><td> 7.290236</td><td> 7.259300</td><td>RAD17    </td></tr>
	<tr><td>ILMN_1651303</td><td> 7.222736</td><td> 7.247647</td><td> 7.232780</td><td> 7.297565</td><td> 7.254969</td><td> 7.236302</td><td> 7.195875</td><td> 7.268146</td><td> 7.296475</td><td>⋯</td><td> 7.252944</td><td> 7.236959</td><td> 7.238020</td><td> 7.237577</td><td> 7.268429</td><td> 7.219242</td><td> 7.235599</td><td> 7.223959</td><td> 7.220732</td><td>ATP13A4  </td></tr>
	<tr><td>ILMN_1651316</td><td> 7.717240</td><td> 7.575777</td><td> 7.351520</td><td> 7.757911</td><td> 7.424758</td><td> 7.413289</td><td> 7.325682</td><td> 7.336730</td><td> 7.603029</td><td>⋯</td><td> 8.065558</td><td> 7.878759</td><td> 7.845327</td><td> 7.502832</td><td> 7.444193</td><td> 7.576983</td><td> 7.725902</td><td> 8.354902</td><td> 7.665653</td><td>CD69     </td></tr>
	<tr><td>ILMN_1651330</td><td> 7.328595</td><td> 7.272354</td><td> 7.334431</td><td> 7.290530</td><td> 7.298802</td><td> 7.353108</td><td> 7.308056</td><td> 7.301047</td><td> 7.324117</td><td>⋯</td><td> 7.310028</td><td> 7.305109</td><td> 7.287330</td><td> 7.247985</td><td> 7.284730</td><td> 7.311427</td><td> 7.270429</td><td> 7.310882</td><td> 7.295372</td><td>KCNG4    </td></tr>
	<tr><td>ILMN_1651336</td><td> 7.797771</td><td> 7.661095</td><td> 8.029717</td><td> 7.694936</td><td> 7.826081</td><td> 7.761593</td><td> 7.757811</td><td> 7.846034</td><td> 8.326006</td><td>⋯</td><td> 7.669512</td><td> 7.634384</td><td> 7.763549</td><td> 7.736182</td><td> 7.793086</td><td> 7.726998</td><td> 7.775012</td><td> 7.690243</td><td> 7.763138</td><td>MLYCD    </td></tr>
	<tr><td>ILMN_1651339</td><td> 7.324689</td><td> 7.318143</td><td> 7.352173</td><td> 7.345603</td><td> 7.352943</td><td> 7.300459</td><td> 7.416398</td><td> 7.292177</td><td> 7.318286</td><td>⋯</td><td> 7.327692</td><td> 7.325016</td><td> 7.303049</td><td> 7.334726</td><td> 7.340936</td><td> 7.249891</td><td> 7.339242</td><td> 7.335891</td><td> 7.293805</td><td>KIAA0701 </td></tr>
	<tr><td>ILMN_1651343</td><td>11.046277</td><td>10.596099</td><td>11.763076</td><td>10.866014</td><td>10.504696</td><td>10.532986</td><td>10.978685</td><td>10.949018</td><td>10.661556</td><td>⋯</td><td> 8.696550</td><td> 9.684871</td><td> 9.175033</td><td> 9.533054</td><td> 9.130688</td><td> 9.937437</td><td>10.008520</td><td> 9.170616</td><td>10.219657</td><td>ITGA11   </td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>ILMN_1815719</td><td> 9.179718</td><td> 9.587973</td><td> 7.750780</td><td> 9.546750</td><td> 8.997113</td><td> 9.325212</td><td> 8.590985</td><td> 8.110609</td><td> 8.772737</td><td>⋯</td><td>10.245688</td><td> 9.796820</td><td> 9.874710</td><td> 9.903536</td><td>10.358298</td><td>10.256998</td><td>10.033637</td><td>10.388640</td><td>10.324245</td><td>PLCG2   </td></tr>
	<tr><td>ILMN_1815723</td><td> 7.784783</td><td> 7.728034</td><td> 8.009661</td><td> 7.861535</td><td> 7.993029</td><td> 7.736773</td><td> 8.340097</td><td> 7.973944</td><td> 8.024763</td><td>⋯</td><td> 8.002647</td><td> 8.033048</td><td> 7.884821</td><td> 7.940472</td><td> 7.899521</td><td> 7.937791</td><td> 8.100283</td><td> 7.937791</td><td> 7.637534</td><td>NUP35   </td></tr>
	<tr><td>ILMN_1815733</td><td> 9.800157</td><td> 9.857971</td><td> 9.484091</td><td> 9.775808</td><td> 9.554961</td><td>10.061231</td><td> 9.620250</td><td> 9.483040</td><td> 9.488859</td><td>⋯</td><td> 9.856830</td><td> 9.083321</td><td> 9.065944</td><td> 9.789384</td><td> 9.850954</td><td> 9.716034</td><td> 9.601086</td><td> 9.643289</td><td> 9.600608</td><td>EIF5    </td></tr>
	<tr><td>ILMN_1815734</td><td>10.573525</td><td>10.239281</td><td>10.616719</td><td>10.658819</td><td>10.372935</td><td>10.520772</td><td>10.282220</td><td>10.683957</td><td>10.813227</td><td>⋯</td><td>10.098576</td><td>10.226348</td><td>10.162095</td><td>10.264295</td><td> 9.942608</td><td>10.310882</td><td>10.321336</td><td> 9.964624</td><td> 9.938990</td><td>FCHSD2  </td></tr>
	<tr><td>ILMN_1815745</td><td> 8.935056</td><td> 8.858813</td><td> 9.254408</td><td> 8.488450</td><td> 8.774254</td><td> 8.824948</td><td> 8.681662</td><td> 8.998565</td><td> 8.924661</td><td>⋯</td><td> 8.787256</td><td> 8.517182</td><td> 8.810712</td><td> 8.820930</td><td> 8.748685</td><td> 8.952298</td><td> 8.816727</td><td> 8.602641</td><td> 8.405424</td><td>SOX4    </td></tr>
	<tr><td>ILMN_1815759</td><td> 7.896605</td><td> 7.987875</td><td> 7.962446</td><td> 7.915164</td><td> 7.914356</td><td> 7.895223</td><td> 8.062391</td><td> 7.952463</td><td> 7.807262</td><td>⋯</td><td> 8.104837</td><td> 8.297802</td><td> 8.158375</td><td> 8.149317</td><td> 8.132279</td><td> 8.000541</td><td> 8.135865</td><td> 8.013749</td><td> 7.969219</td><td>CTDP1   </td></tr>
	<tr><td>ILMN_1815773</td><td> 7.262818</td><td> 7.256086</td><td> 7.271000</td><td> 7.247518</td><td> 7.249999</td><td> 7.243351</td><td> 7.304826</td><td> 7.229646</td><td> 7.269641</td><td>⋯</td><td> 7.254027</td><td> 7.247485</td><td> 7.251096</td><td> 7.227695</td><td> 7.221567</td><td> 7.269739</td><td> 7.272506</td><td> 7.223959</td><td> 7.254033</td><td>MDAC1   </td></tr>
	<tr><td>ILMN_1815780</td><td> 7.204965</td><td> 7.221138</td><td> 7.193708</td><td> 7.217951</td><td> 7.400085</td><td> 7.214288</td><td> 7.440597</td><td> 7.244921</td><td> 7.287981</td><td>⋯</td><td> 7.325452</td><td> 7.260323</td><td> 7.204522</td><td> 7.243569</td><td> 7.203947</td><td> 7.216379</td><td> 7.286438</td><td> 7.246184</td><td> 7.184053</td><td>MAPK10  </td></tr>
	<tr><td>ILMN_1815795</td><td> 8.713454</td><td> 9.229663</td><td> 8.590365</td><td> 8.660969</td><td> 8.850444</td><td> 8.695236</td><td> 8.770857</td><td> 8.659100</td><td> 8.557428</td><td>⋯</td><td> 8.980781</td><td> 9.246998</td><td> 9.188676</td><td> 9.001238</td><td> 9.310173</td><td> 8.996418</td><td> 9.151311</td><td> 9.010747</td><td> 9.021468</td><td>TSSC1   </td></tr>
	<tr><td>ILMN_1815796</td><td> 7.261643</td><td> 7.228978</td><td> 7.260609</td><td> 7.249508</td><td> 7.242441</td><td> 7.247732</td><td> 7.169445</td><td> 7.269726</td><td> 7.246200</td><td>⋯</td><td> 7.304888</td><td> 7.264362</td><td> 7.234619</td><td> 7.290411</td><td> 7.272024</td><td> 7.227777</td><td> 7.257393</td><td> 7.236692</td><td> 7.285235</td><td>NT5C1A  </td></tr>
	<tr><td>ILMN_1815812</td><td> 8.092522</td><td> 7.951864</td><td> 8.576777</td><td> 8.017761</td><td> 8.077110</td><td> 7.816763</td><td> 8.297483</td><td> 8.539556</td><td> 8.264162</td><td>⋯</td><td> 7.726238</td><td> 7.799419</td><td> 7.764777</td><td> 7.810924</td><td> 7.735018</td><td> 7.722065</td><td> 7.754737</td><td> 7.820164</td><td> 7.894886</td><td>C16orf5 </td></tr>
	<tr><td>ILMN_1815840</td><td> 7.340552</td><td> 7.334196</td><td> 7.405242</td><td> 7.353365</td><td> 7.328280</td><td> 7.359995</td><td> 7.333042</td><td> 7.365936</td><td> 7.334075</td><td>⋯</td><td> 7.299850</td><td> 7.302125</td><td> 7.327393</td><td> 7.323107</td><td> 7.414249</td><td> 7.314594</td><td> 7.332471</td><td> 7.369037</td><td> 7.463679</td><td>IGLL1   </td></tr>
	<tr><td>ILMN_1815849</td><td> 7.242466</td><td> 7.177864</td><td> 7.192287</td><td> 7.177105</td><td> 7.154340</td><td> 7.172473</td><td> 7.220077</td><td> 7.180943</td><td> 7.178536</td><td>⋯</td><td> 7.166032</td><td> 7.164261</td><td> 7.201924</td><td> 7.165135</td><td> 7.209092</td><td> 7.177843</td><td> 7.171619</td><td> 7.184355</td><td> 7.183927</td><td>OR8B2   </td></tr>
	<tr><td>ILMN_1815859</td><td> 7.961346</td><td> 8.055121</td><td> 7.959272</td><td> 8.099019</td><td> 7.954359</td><td> 8.042481</td><td> 8.085968</td><td> 8.201399</td><td> 8.159219</td><td>⋯</td><td> 8.104294</td><td> 8.053927</td><td> 8.177804</td><td> 8.154593</td><td> 8.180488</td><td> 7.990242</td><td> 8.015188</td><td> 7.886453</td><td> 8.367435</td><td>ERCC2   </td></tr>
	<tr><td>ILMN_1815874</td><td> 8.997746</td><td> 9.171955</td><td> 8.675459</td><td> 8.833961</td><td> 9.099671</td><td> 8.988140</td><td> 8.959330</td><td> 8.706638</td><td> 8.940344</td><td>⋯</td><td> 9.136816</td><td> 9.382961</td><td> 9.232052</td><td> 9.201560</td><td> 9.234646</td><td> 9.266216</td><td> 9.291357</td><td> 9.166262</td><td> 9.272135</td><td>NANS    </td></tr>
	<tr><td>ILMN_1815878</td><td>10.589751</td><td>10.736022</td><td>10.583114</td><td>10.369826</td><td>10.556491</td><td>10.951441</td><td>10.559109</td><td>10.240076</td><td>10.038857</td><td>⋯</td><td>11.160426</td><td>10.829010</td><td>11.142211</td><td>10.993595</td><td>11.284070</td><td>11.075648</td><td>10.808434</td><td>10.636503</td><td>10.896888</td><td>C11orf59</td></tr>
	<tr><td>ILMN_1815882</td><td> 8.077857</td><td> 7.917237</td><td> 8.374968</td><td> 8.135050</td><td> 8.446158</td><td> 8.062890</td><td> 8.616791</td><td> 8.418736</td><td> 8.396691</td><td>⋯</td><td> 8.332349</td><td> 8.053686</td><td> 8.094642</td><td> 8.094066</td><td> 7.826081</td><td> 8.230748</td><td> 8.045888</td><td> 7.915406</td><td> 7.711444</td><td>HNRPA1  </td></tr>
	<tr><td>ILMN_1815885</td><td> 8.099842</td><td> 8.084597</td><td> 8.281198</td><td> 8.093233</td><td> 8.040160</td><td> 8.099408</td><td> 8.124230</td><td> 8.349547</td><td> 8.488450</td><td>⋯</td><td> 7.947613</td><td> 8.007087</td><td> 7.964730</td><td> 8.046987</td><td> 8.038710</td><td> 8.193979</td><td> 8.025866</td><td> 7.885480</td><td> 8.142848</td><td>ZNF228  </td></tr>
	<tr><td>ILMN_1815890</td><td> 7.476659</td><td> 7.580629</td><td> 7.214657</td><td> 7.563308</td><td> 7.463240</td><td> 7.428503</td><td> 7.313171</td><td> 7.311707</td><td> 7.460932</td><td>⋯</td><td> 7.789116</td><td> 7.641788</td><td> 7.707770</td><td> 7.527691</td><td> 7.596968</td><td> 7.615278</td><td> 7.591839</td><td> 7.538505</td><td> 7.653301</td><td>IL12RB1 </td></tr>
	<tr><td>ILMN_1815908</td><td> 8.210155</td><td> 7.832124</td><td> 8.396115</td><td> 8.134584</td><td> 8.022541</td><td> 7.902125</td><td> 8.321463</td><td> 8.222483</td><td> 7.936109</td><td>⋯</td><td> 7.824096</td><td> 7.751923</td><td> 7.582565</td><td> 7.599136</td><td> 7.571474</td><td> 8.009066</td><td> 7.879705</td><td> 7.596559</td><td> 7.880327</td><td>GYPE    </td></tr>
	<tr><td>ILMN_1815923</td><td> 7.170334</td><td> 7.215901</td><td> 7.202828</td><td> 7.232299</td><td> 7.237880</td><td> 7.244166</td><td> 7.246504</td><td> 7.296397</td><td> 7.233887</td><td>⋯</td><td> 7.278047</td><td> 7.168953</td><td> 7.241946</td><td> 7.231832</td><td> 7.182938</td><td> 7.206737</td><td> 7.215314</td><td> 7.214067</td><td> 7.258633</td><td>SMCR7   </td></tr>
	<tr><td>ILMN_1815924</td><td> 9.365133</td><td> 9.353383</td><td> 9.375704</td><td> 9.162812</td><td> 9.223520</td><td> 9.206855</td><td> 9.251731</td><td> 9.358738</td><td> 9.790366</td><td>⋯</td><td> 9.203148</td><td> 9.157930</td><td> 9.289189</td><td> 9.343734</td><td> 9.388720</td><td> 9.297218</td><td> 9.462036</td><td> 9.535842</td><td> 9.318426</td><td>NUP107  </td></tr>
	<tr><td>ILMN_1815933</td><td> 8.087673</td><td> 7.976249</td><td> 7.943099</td><td> 8.108333</td><td> 7.924111</td><td> 8.065171</td><td> 7.832018</td><td> 8.005984</td><td> 8.163325</td><td>⋯</td><td> 8.094740</td><td> 8.142293</td><td> 8.158009</td><td> 8.279694</td><td> 8.286033</td><td> 8.022717</td><td> 8.134822</td><td> 8.025610</td><td> 8.474256</td><td>FTSJ2   </td></tr>
	<tr><td>ILMN_1815937</td><td> 7.252676</td><td> 7.293665</td><td> 7.371083</td><td> 7.311351</td><td> 7.282919</td><td> 7.295841</td><td> 7.259993</td><td> 7.285671</td><td> 7.315877</td><td>⋯</td><td> 7.257732</td><td> 7.296783</td><td> 7.303240</td><td> 7.310290</td><td> 7.265408</td><td> 7.283729</td><td> 7.295294</td><td> 7.315117</td><td> 7.275810</td><td>MGC9712 </td></tr>
	<tr><td>ILMN_1815938</td><td> 7.230973</td><td> 7.251768</td><td> 7.250046</td><td> 7.253949</td><td> 7.210804</td><td> 7.240074</td><td> 7.236095</td><td> 7.252055</td><td> 7.213958</td><td>⋯</td><td> 7.238224</td><td> 7.244267</td><td> 7.272532</td><td> 7.244075</td><td> 7.215729</td><td> 7.212128</td><td> 7.207027</td><td> 7.223263</td><td> 7.258715</td><td>TRPM3   </td></tr>
	<tr><td>ILMN_1815941</td><td> 9.798893</td><td> 9.694271</td><td> 9.907072</td><td> 9.768732</td><td> 9.895686</td><td> 9.764151</td><td>10.004447</td><td>10.544573</td><td> 9.749108</td><td>⋯</td><td> 9.098114</td><td> 9.498280</td><td> 9.132293</td><td> 9.248486</td><td> 8.758983</td><td> 9.800157</td><td> 9.500341</td><td> 9.520802</td><td> 9.415060</td><td>SMAD7   </td></tr>
	<tr><td>ILMN_1815951</td><td> 7.966225</td><td> 7.953338</td><td> 7.804860</td><td> 7.807950</td><td> 7.807399</td><td> 7.812572</td><td> 7.878314</td><td> 7.683947</td><td> 7.699060</td><td>⋯</td><td> 7.762312</td><td> 7.762570</td><td> 7.698209</td><td> 7.607156</td><td> 7.593200</td><td> 7.693124</td><td> 7.747622</td><td> 7.677472</td><td> 7.663077</td><td>PCYOX1L </td></tr>
	<tr><td>ILMN_2038774</td><td>14.928846</td><td>15.050317</td><td>15.165482</td><td>15.137000</td><td>15.363611</td><td>14.908203</td><td>15.397685</td><td>15.192185</td><td>15.117043</td><td>⋯</td><td>15.290683</td><td>15.363611</td><td>15.363611</td><td>15.137000</td><td>15.253761</td><td>15.326456</td><td>15.397685</td><td>15.326456</td><td>14.885728</td><td>EEF1A1  </td></tr>
	<tr><td>ILMN_2038777</td><td>15.004360</td><td>15.290683</td><td>14.758077</td><td>15.165482</td><td>15.192185</td><td>14.733469</td><td>15.137000</td><td>15.029515</td><td>14.664041</td><td>⋯</td><td>15.029515</td><td>15.253761</td><td>14.837210</td><td>15.029515</td><td>14.634348</td><td>14.908203</td><td>15.290683</td><td>15.290683</td><td>14.712066</td><td>ACTB    </td></tr>
	<tr><td>ILMN_2038778</td><td>14.077380</td><td>13.952513</td><td>13.661227</td><td>13.999062</td><td>13.964369</td><td>14.229570</td><td>14.098730</td><td>13.888182</td><td>13.194778</td><td>⋯</td><td>13.850222</td><td>14.684127</td><td>14.216344</td><td>14.206824</td><td>14.335820</td><td>14.066954</td><td>14.265875</td><td>14.010573</td><td>13.734694</td><td>GAPDH   </td></tr>
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
<ol class=list-inline><li>'GSM4973394'</li><li>'GSM4973395'</li><li>'GSM4973397'</li><li>'GSM4973398'</li><li>'GSM4973400'</li><li>'GSM4973401'</li><li>'GSM4973403'</li><li>'GSM4973404'</li><li>'GSM4973406'</li><li>'GSM4973408'</li><li>'GSM4973409'</li><li>'GSM4973411'</li><li>'GSM4973412'</li><li>'GSM4973414'</li><li>'GSM4973415'</li><li>'GSM4973417'</li><li>'GSM4973418'</li><li>'GSM4973419'</li><li>'GSM4973421'</li><li>'GSM4973422'</li><li>'GSM4973424'</li><li>'GSM4973425'</li><li>'GSM4973427'</li><li>'GSM4973429'</li><li>'GSM4973430'</li><li>'GSM4973432'</li><li>'GSM4973433'</li><li>'GSM4973435'</li><li>'GSM4973436'</li><li>'GSM4973438'</li><li>'GSM4973439'</li><li>'GSM4973440'</li><li>'GSM4973441'</li><li>'GSM4973443'</li><li>'GSM4973444'</li><li>'GSM4973445'</li><li>'GSM4973446'</li><li>'GSM4973448'</li><li>'GSM4973449'</li><li>'GSM4973451'</li><li>'GSM4973452'</li><li>'GSM4973454'</li><li>'GSM4973456'</li></ol>




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


'GSM4973394=mean(GSM4973394),GSM4973395=mean(GSM4973395),GSM4973397=mean(GSM4973397),GSM4973398=mean(GSM4973398),GSM4973400=mean(GSM4973400),GSM4973401=mean(GSM4973401),GSM4973403=mean(GSM4973403),GSM4973404=mean(GSM4973404),GSM4973406=mean(GSM4973406),GSM4973408=mean(GSM4973408),GSM4973409=mean(GSM4973409),GSM4973411=mean(GSM4973411),GSM4973412=mean(GSM4973412),GSM4973414=mean(GSM4973414),GSM4973415=mean(GSM4973415),GSM4973417=mean(GSM4973417),GSM4973418=mean(GSM4973418),GSM4973419=mean(GSM4973419),GSM4973421=mean(GSM4973421),GSM4973422=mean(GSM4973422),GSM4973424=mean(GSM4973424),GSM4973425=mean(GSM4973425),GSM4973427=mean(GSM4973427),GSM4973429=mean(GSM4973429),GSM4973430=mean(GSM4973430),GSM4973432=mean(GSM4973432),GSM4973433=mean(GSM4973433),GSM4973435=mean(GSM4973435),GSM4973436=mean(GSM4973436),GSM4973438=mean(GSM4973438),GSM4973439=mean(GSM4973439),GSM4973440=mean(GSM4973440),GSM4973441=mean(GSM4973441),GSM4973443=mean(GSM4973443),GSM4973444=mean(GSM4973444),GSM4973445=mean(GSM4973445),GSM4973446=mean(GSM4973446),GSM4973448=mean(GSM4973448),GSM4973449=mean(GSM4973449),GSM4973451=mean(GSM4973451),GSM4973452=mean(GSM4973452),GSM4973454=mean(GSM4973454),GSM4973456=mean(GSM4973456),'



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
<caption>A data.frame: 6 × 44</caption>
<thead>
	<tr><th></th><th scope=col>GSM4973394</th><th scope=col>GSM4973395</th><th scope=col>GSM4973397</th><th scope=col>GSM4973398</th><th scope=col>GSM4973400</th><th scope=col>GSM4973401</th><th scope=col>GSM4973403</th><th scope=col>GSM4973404</th><th scope=col>GSM4973406</th><th scope=col>GSM4973408</th><th scope=col>⋯</th><th scope=col>GSM4973444</th><th scope=col>GSM4973445</th><th scope=col>GSM4973446</th><th scope=col>GSM4973448</th><th scope=col>GSM4973449</th><th scope=col>GSM4973451</th><th scope=col>GSM4973452</th><th scope=col>GSM4973454</th><th scope=col>GSM4973456</th><th scope=col>symbol</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>15.397685</td><td>15.397685</td><td>15.363611</td><td>15.397685</td><td>15.326456</td><td>15.363611</td><td>15.326456</td><td>15.397685</td><td>15.363611</td><td>15.253761</td><td>⋯</td><td>15.451841</td><td>15.222846</td><td>15.451841</td><td>15.397685</td><td>15.451841</td><td>15.363611</td><td>15.363611</td><td>15.363611</td><td>15.451841</td><td>EEF1A1  </td></tr>
	<tr><th scope=row>2</th><td>10.574274</td><td>11.283223</td><td>10.932375</td><td>10.952329</td><td>10.620145</td><td>10.958478</td><td>10.998223</td><td>11.137739</td><td>11.108668</td><td>10.572747</td><td>⋯</td><td>10.451724</td><td>10.661556</td><td>10.443294</td><td>10.458157</td><td>10.990131</td><td>10.650094</td><td>10.908000</td><td>10.898592</td><td>11.075648</td><td>TUBB    </td></tr>
	<tr><th scope=row>3</th><td>12.170338</td><td>12.394059</td><td>12.454637</td><td>11.881277</td><td>12.152035</td><td>12.514034</td><td>12.235775</td><td>12.034007</td><td>12.124991</td><td>11.795001</td><td>⋯</td><td>12.599769</td><td>12.230211</td><td>12.443148</td><td>12.536527</td><td>12.916573</td><td>12.199100</td><td>12.471040</td><td>12.786596</td><td>11.730697</td><td>TXN     </td></tr>
	<tr><th scope=row>4</th><td>14.837210</td><td>14.979876</td><td>15.090183</td><td>15.029515</td><td>15.050317</td><td>15.192185</td><td>15.004360</td><td>15.070705</td><td>14.908203</td><td>15.029515</td><td>⋯</td><td>15.253761</td><td>15.070705</td><td>15.029515</td><td>15.222846</td><td>15.070705</td><td>15.050317</td><td>14.908203</td><td>15.029515</td><td>15.137000</td><td>ACTB    </td></tr>
	<tr><th scope=row>5</th><td> 7.327692</td><td> 7.297238</td><td> 7.368454</td><td> 7.440558</td><td> 7.374302</td><td> 7.407157</td><td> 7.238493</td><td> 7.511092</td><td> 7.356138</td><td> 7.465387</td><td>⋯</td><td> 7.425922</td><td> 7.321693</td><td> 7.347743</td><td> 7.411170</td><td> 7.320686</td><td> 7.290813</td><td> 7.320089</td><td> 7.324418</td><td> 7.396416</td><td>SLC35E2 </td></tr>
	<tr><th scope=row>6</th><td> 7.339957</td><td> 7.451542</td><td> 7.425922</td><td> 7.520726</td><td> 7.608858</td><td> 7.471193</td><td> 7.506215</td><td> 7.344494</td><td> 7.435864</td><td> 7.448614</td><td>⋯</td><td> 8.146389</td><td> 7.818563</td><td> 7.636750</td><td> 7.790764</td><td> 7.688711</td><td> 7.683776</td><td> 7.784164</td><td> 7.926848</td><td> 7.472395</td><td>PDCD1LG2</td></tr>
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
              summarise(GSM4973394=mean(GSM4973394),GSM4973395=mean(GSM4973395),GSM4973397=mean(GSM4973397),GSM4973398=mean(GSM4973398),GSM4973400=mean(GSM4973400),GSM4973401=mean(GSM4973401),GSM4973403=mean(GSM4973403),GSM4973404=mean(GSM4973404),GSM4973406=mean(GSM4973406),GSM4973408=mean(GSM4973408),GSM4973409=mean(GSM4973409),GSM4973411=mean(GSM4973411),GSM4973412=mean(GSM4973412),GSM4973414=mean(GSM4973414),GSM4973415=mean(GSM4973415),GSM4973417=mean(GSM4973417),GSM4973418=mean(GSM4973418),GSM4973419=mean(GSM4973419),GSM4973421=mean(GSM4973421),GSM4973422=mean(GSM4973422),GSM4973424=mean(GSM4973424),GSM4973425=mean(GSM4973425),GSM4973427=mean(GSM4973427),GSM4973429=mean(GSM4973429),GSM4973430=mean(GSM4973430),GSM4973432=mean(GSM4973432),GSM4973433=mean(GSM4973433),GSM4973435=mean(GSM4973435),GSM4973436=mean(GSM4973436),GSM4973438=mean(GSM4973438),GSM4973439=mean(GSM4973439),GSM4973440=mean(GSM4973440),GSM4973441=mean(GSM4973441),GSM4973443=mean(GSM4973443),GSM4973444=mean(GSM4973444),GSM4973445=mean(GSM4973445),GSM4973446=mean(GSM4973446),GSM4973448=mean(GSM4973448),GSM4973449=mean(GSM4973449),GSM4973451=mean(GSM4973451),GSM4973452=mean(GSM4973452),GSM4973454=mean(GSM4973454),GSM4973456=mean(GSM4973456))  %>% 
                column_to_rownames(var ='symbol')
```


```R
gset$GSE163154_series_matrix.txt.gz@phenoData@data  %>%  head() 
```


<table class="dataframe">
<caption>A data.frame: 6 × 32</caption>
<thead>
	<tr><th></th><th scope=col>title</th><th scope=col>geo_accession</th><th scope=col>status</th><th scope=col>submission_date</th><th scope=col>last_update_date</th><th scope=col>type</th><th scope=col>channel_count</th><th scope=col>source_name_ch1</th><th scope=col>organism_ch1</th><th scope=col>characteristics_ch1</th><th scope=col>⋯</th><th scope=col>contact_name</th><th scope=col>contact_email</th><th scope=col>contact_institute</th><th scope=col>contact_address</th><th scope=col>contact_city</th><th scope=col>contact_zip/postal_code</th><th scope=col>contact_country</th><th scope=col>supplementary_file</th><th scope=col>data_row_count</th><th scope=col>classification:ch1</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM4973394</th><td>non-IPH 1</td><td>GSM4973394</td><td>Public on Dec 15 2020</td><td>Dec 14 2020</td><td>Dec 15 2020</td><td>RNA</td><td>1</td><td>carotid plaque</td><td>Homo sapiens</td><td>classification: non-IPH</td><td>⋯</td><td>Han,,Jin</td><td>han.jin@scilifelab.se</td><td>Maastricht University</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>6229 HX</td><td>Netherlands</td><td>NONE</td><td>22184</td><td>non-IPH</td></tr>
	<tr><th scope=row>GSM4973395</th><td>non-IPH 2</td><td>GSM4973395</td><td>Public on Dec 15 2020</td><td>Dec 14 2020</td><td>Dec 15 2020</td><td>RNA</td><td>1</td><td>carotid plaque</td><td>Homo sapiens</td><td>classification: non-IPH</td><td>⋯</td><td>Han,,Jin</td><td>han.jin@scilifelab.se</td><td>Maastricht University</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>6229 HX</td><td>Netherlands</td><td>NONE</td><td>22184</td><td>non-IPH</td></tr>
	<tr><th scope=row>GSM4973397</th><td>non-IPH 3</td><td>GSM4973397</td><td>Public on Dec 15 2020</td><td>Dec 14 2020</td><td>Dec 15 2020</td><td>RNA</td><td>1</td><td>carotid plaque</td><td>Homo sapiens</td><td>classification: non-IPH</td><td>⋯</td><td>Han,,Jin</td><td>han.jin@scilifelab.se</td><td>Maastricht University</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>6229 HX</td><td>Netherlands</td><td>NONE</td><td>22184</td><td>non-IPH</td></tr>
	<tr><th scope=row>GSM4973398</th><td>non-IPH 4</td><td>GSM4973398</td><td>Public on Dec 15 2020</td><td>Dec 14 2020</td><td>Dec 15 2020</td><td>RNA</td><td>1</td><td>carotid plaque</td><td>Homo sapiens</td><td>classification: non-IPH</td><td>⋯</td><td>Han,,Jin</td><td>han.jin@scilifelab.se</td><td>Maastricht University</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>6229 HX</td><td>Netherlands</td><td>NONE</td><td>22184</td><td>non-IPH</td></tr>
	<tr><th scope=row>GSM4973400</th><td>non-IPH 5</td><td>GSM4973400</td><td>Public on Dec 15 2020</td><td>Dec 14 2020</td><td>Dec 15 2020</td><td>RNA</td><td>1</td><td>carotid plaque</td><td>Homo sapiens</td><td>classification: non-IPH</td><td>⋯</td><td>Han,,Jin</td><td>han.jin@scilifelab.se</td><td>Maastricht University</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>6229 HX</td><td>Netherlands</td><td>NONE</td><td>22184</td><td>non-IPH</td></tr>
	<tr><th scope=row>GSM4973401</th><td>non-IPH 6</td><td>GSM4973401</td><td>Public on Dec 15 2020</td><td>Dec 14 2020</td><td>Dec 15 2020</td><td>RNA</td><td>1</td><td>carotid plaque</td><td>Homo sapiens</td><td>classification: non-IPH</td><td>⋯</td><td>Han,,Jin</td><td>han.jin@scilifelab.se</td><td>Maastricht University</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>6229 HX</td><td>Netherlands</td><td>NONE</td><td>22184</td><td>non-IPH</td></tr>
</tbody>
</table>




```R
group_list  <- factor(gset$GSE163154_series_matrix.txt.gz@phenoData@data[,c('characteristics_ch1')]) 
```


```R
levels(group_list)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'classification: IPH'</li><li>'classification: non-IPH'</li></ol>




```R
group_list  <- fct_recode(group_list, 'IPH' = 'classification: IPH',
                                      'non-IPH' = 'classification: non-IPH')
```


```R
group_list <- relevel(group_list, ref="non-IPH")
```


```R
group_list
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>non-IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li><li>IPH</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'non-IPH'</li><li>'IPH'</li></ol>
</details>



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
gsva_matrix<- gsva(as.matrix(normalized_gset_mean), module_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
```

    Estimating ssGSEA scores for 2 gene sets.
      |======================================================================| 100%
    



```R
library(ggpubr)
```


```R
p_yellow  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       ggboxplot(x = 'group', y = 'Yellow_Module', ylab = 'GSVA Score',title = 'Yellow Module',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('non-IPH','IPH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_yellow
```


    
![png](Step5.3_Two_Group_Compare_GSE163154_files/Step5.3_Two_Group_Compare_GSE163154_27_0.png)
    



```R
p_brown  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       ggboxplot(x = 'group', y = 'Brown_Module', ylab = 'GSVA Score',title = 'Brown Module',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('non-IPH','IPH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_brown
```


    
![png](Step5.3_Two_Group_Compare_GSE163154_files/Step5.3_Two_Group_Compare_GSE163154_28_0.png)
    



```R
library(patchwork)
```


```R
p0 <- p_yellow + p_brown + plot_layout(guides='collect') +  plot_annotation(title = 'GSE163154',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE163154_Two_Module.pdf', height = 6, width = 9, units = 'cm')
```


    
![png](Step5.3_Two_Group_Compare_GSE163154_files/Step5.3_Two_Group_Compare_GSE163154_30_0.png)
    

