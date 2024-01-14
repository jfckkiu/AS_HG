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
gset  <-  GEOquery::getGEO('GSE41571',getGPL = F)
```

    Found 1 file(s)
    
    GSE41571_series_matrix.txt.gz
    
    [1mRows: [22m[34m54675[39m [1mColumns: [22m[34m12[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m "\t"
    [31mchr[39m  (1): ID_REF
    [32mdbl[39m (11): GSM1019539, GSM1019540, GSM1019541, GSM1019542, GSM1019543, GSM101...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
normalized_gset  <- gset$GSE41571_series_matrix.txt.gz@assayData$exprs
```


```R
range(normalized_gset)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>0.01</li><li>3664</li></ol>




```R
normalized_gset
```


<table class="dataframe">
<caption>A matrix: 54675 × 11 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>GSM1019539</th><th scope=col>GSM1019540</th><th scope=col>GSM1019541</th><th scope=col>GSM1019542</th><th scope=col>GSM1019543</th><th scope=col>GSM1019544</th><th scope=col>GSM1019545</th><th scope=col>GSM1019546</th><th scope=col>GSM1019547</th><th scope=col>GSM1019548</th><th scope=col>GSM1019549</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1007_s_at</th><td>0.796</td><td>1.034</td><td>0.752</td><td>1.381</td><td>0.710</td><td>1.064</td><td>1.064</td><td>1.383</td><td>0.778</td><td>0.838</td><td>1.000</td></tr>
	<tr><th scope=row>1053_at</th><td>1.000</td><td>1.004</td><td>2.045</td><td>1.259</td><td>1.024</td><td>0.904</td><td>0.823</td><td>0.992</td><td>0.954</td><td>0.450</td><td>1.231</td></tr>
	<tr><th scope=row>117_at</th><td>5.436</td><td>1.000</td><td>0.694</td><td>1.246</td><td>1.697</td><td>5.746</td><td>0.733</td><td>1.030</td><td>0.907</td><td>0.519</td><td>0.986</td></tr>
	<tr><th scope=row>121_at</th><td>1.109</td><td>0.911</td><td>1.390</td><td>1.417</td><td>1.187</td><td>0.938</td><td>1.404</td><td>0.530</td><td>0.478</td><td>1.000</td><td>0.776</td></tr>
	<tr><th scope=row>1255_g_at</th><td>1.000</td><td>1.008</td><td>1.139</td><td>1.018</td><td>1.046</td><td>0.997</td><td>0.991</td><td>0.992</td><td>0.976</td><td>1.013</td><td>0.982</td></tr>
	<tr><th scope=row>1294_at</th><td>1.003</td><td>0.979</td><td>1.014</td><td>1.048</td><td>1.040</td><td>0.866</td><td>0.896</td><td>2.246</td><td>1.000</td><td>0.826</td><td>0.980</td></tr>
	<tr><th scope=row>1316_at</th><td>1.365</td><td>1.108</td><td>0.778</td><td>1.000</td><td>0.985</td><td>0.965</td><td>2.252</td><td>0.822</td><td>1.161</td><td>1.189</td><td>0.607</td></tr>
	<tr><th scope=row>1320_at</th><td>0.968</td><td>1.000</td><td>1.191</td><td>1.037</td><td>1.056</td><td>0.977</td><td>0.956</td><td>1.023</td><td>0.962</td><td>1.018</td><td>0.984</td></tr>
	<tr><th scope=row>1405_i_at</th><td>1.112</td><td>2.209</td><td>1.038</td><td>1.000</td><td>0.560</td><td>0.626</td><td>3.283</td><td>0.584</td><td>0.826</td><td>0.482</td><td>1.194</td></tr>
	<tr><th scope=row>1431_at</th><td>0.898</td><td>0.942</td><td>1.206</td><td>0.981</td><td>1.000</td><td>1.021</td><td>0.919</td><td>1.010</td><td>1.001</td><td>0.956</td><td>1.127</td></tr>
	<tr><th scope=row>1438_at</th><td>0.990</td><td>1.000</td><td>1.102</td><td>1.011</td><td>1.023</td><td>0.976</td><td>0.990</td><td>1.012</td><td>0.981</td><td>1.008</td><td>0.988</td></tr>
	<tr><th scope=row>1487_at</th><td>1.246</td><td>0.854</td><td>1.074</td><td>1.251</td><td>1.331</td><td>0.938</td><td>0.912</td><td>1.000</td><td>0.703</td><td>0.593</td><td>1.119</td></tr>
	<tr><th scope=row>1494_f_at</th><td>0.989</td><td>1.000</td><td>1.117</td><td>1.033</td><td>1.025</td><td>0.997</td><td>0.997</td><td>1.005</td><td>0.975</td><td>1.032</td><td>0.998</td></tr>
	<tr><th scope=row>1552256_a_at</th><td>1.100</td><td>1.100</td><td>1.089</td><td>1.339</td><td>1.000</td><td>0.496</td><td>0.206</td><td>0.362</td><td>0.163</td><td>0.172</td><td>1.026</td></tr>
	<tr><th scope=row>1552257_a_at</th><td>0.984</td><td>0.786</td><td>1.401</td><td>1.356</td><td>1.058</td><td>0.991</td><td>1.565</td><td>1.000</td><td>0.563</td><td>0.403</td><td>1.603</td></tr>
	<tr><th scope=row>1552258_at</th><td>2.001</td><td>1.000</td><td>0.670</td><td>1.786</td><td>1.016</td><td>0.854</td><td>0.691</td><td>0.551</td><td>1.000</td><td>1.000</td><td>0.550</td></tr>
	<tr><th scope=row>1552261_at</th><td>0.988</td><td>1.036</td><td>1.061</td><td>1.133</td><td>1.000</td><td>0.970</td><td>0.904</td><td>1.207</td><td>0.884</td><td>1.002</td><td>0.966</td></tr>
	<tr><th scope=row>1552263_at</th><td>0.413</td><td>1.220</td><td>1.700</td><td>1.480</td><td>1.628</td><td>0.783</td><td>1.098</td><td>0.445</td><td>1.000</td><td>1.000</td><td>0.861</td></tr>
	<tr><th scope=row>1552264_a_at</th><td>0.777</td><td>1.282</td><td>2.295</td><td>1.227</td><td>1.658</td><td>1.000</td><td>0.978</td><td>0.537</td><td>1.362</td><td>0.790</td><td>0.809</td></tr>
	<tr><th scope=row>1552266_at</th><td>0.842</td><td>0.965</td><td>1.330</td><td>1.021</td><td>1.069</td><td>1.400</td><td>0.973</td><td>0.921</td><td>1.000</td><td>1.000</td><td>0.901</td></tr>
	<tr><th scope=row>1552269_at</th><td>0.959</td><td>0.967</td><td>1.085</td><td>1.000</td><td>1.000</td><td>1.001</td><td>1.049</td><td>1.216</td><td>1.000</td><td>1.000</td><td>0.994</td></tr>
	<tr><th scope=row>1552271_at</th><td>0.982</td><td>0.990</td><td>1.020</td><td>1.002</td><td>1.155</td><td>0.990</td><td>1.277</td><td>1.319</td><td>0.979</td><td>1.000</td><td>0.983</td></tr>
	<tr><th scope=row>1552272_a_at</th><td>1.103</td><td>1.000</td><td>1.096</td><td>1.025</td><td>1.024</td><td>0.974</td><td>0.982</td><td>0.972</td><td>0.962</td><td>1.010</td><td>0.999</td></tr>
	<tr><th scope=row>1552274_at</th><td>1.000</td><td>0.861</td><td>0.781</td><td>0.608</td><td>0.896</td><td>1.286</td><td>1.797</td><td>1.332</td><td>1.443</td><td>3.070</td><td>0.652</td></tr>
	<tr><th scope=row>1552275_s_at</th><td>0.863</td><td>0.930</td><td>1.500</td><td>1.000</td><td>1.000</td><td>1.331</td><td>0.876</td><td>0.707</td><td>1.403</td><td>0.827</td><td>1.157</td></tr>
	<tr><th scope=row>1552276_a_at</th><td>0.994</td><td>1.030</td><td>1.136</td><td>1.089</td><td>1.000</td><td>0.978</td><td>0.883</td><td>1.447</td><td>0.876</td><td>1.015</td><td>0.961</td></tr>
	<tr><th scope=row>1552277_a_at</th><td>1.313</td><td>1.000</td><td>0.552</td><td>1.186</td><td>0.632</td><td>0.901</td><td>2.144</td><td>0.404</td><td>1.597</td><td>2.330</td><td>0.436</td></tr>
	<tr><th scope=row>1552278_a_at</th><td>1.000</td><td>1.745</td><td>4.224</td><td>1.008</td><td>3.187</td><td>0.888</td><td>0.729</td><td>0.606</td><td>0.673</td><td>0.520</td><td>1.111</td></tr>
	<tr><th scope=row>1552279_a_at</th><td>1.007</td><td>1.000</td><td>1.393</td><td>1.042</td><td>2.194</td><td>0.940</td><td>0.834</td><td>0.816</td><td>0.567</td><td>0.726</td><td>1.361</td></tr>
	<tr><th scope=row>1552280_at</th><td>0.938</td><td>4.269</td><td>1.116</td><td>1.326</td><td>1.017</td><td>0.780</td><td>1.000</td><td>0.807</td><td>0.872</td><td>0.850</td><td>1.352</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>AFFX-PheX-3_at</th><td>1.650</td><td>0.667</td><td>0.630</td><td>1.000</td><td>1.256</td><td>1.799</td><td>1.356</td><td>1.409</td><td>0.748</td><td>0.908</td><td>0.625</td></tr>
	<tr><th scope=row>AFFX-PheX-5_at</th><td>1.947</td><td>1.000</td><td>0.761</td><td>1.218</td><td>1.197</td><td>2.005</td><td>0.874</td><td>1.983</td><td>0.757</td><td>0.676</td><td>0.717</td></tr>
	<tr><th scope=row>AFFX-PheX-M_at</th><td>1.843</td><td>0.750</td><td>0.517</td><td>1.000</td><td>1.020</td><td>1.705</td><td>1.022</td><td>1.461</td><td>0.631</td><td>0.709</td><td>0.505</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-dap-3_at</th><td>1.403</td><td>0.902</td><td>0.809</td><td>1.000</td><td>1.000</td><td>1.720</td><td>1.654</td><td>1.399</td><td>0.983</td><td>0.845</td><td>0.591</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-dap-5_at</th><td>1.683</td><td>1.075</td><td>1.000</td><td>0.715</td><td>0.661</td><td>2.056</td><td>2.106</td><td>2.362</td><td>0.584</td><td>0.406</td><td>0.226</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-dap-M_at</th><td>2.564</td><td>0.944</td><td>0.741</td><td>1.000</td><td>0.836</td><td>3.294</td><td>1.833</td><td>1.550</td><td>1.143</td><td>0.760</td><td>0.474</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-lys-3_at</th><td>2.193</td><td>0.740</td><td>0.523</td><td>0.896</td><td>0.774</td><td>3.148</td><td>2.830</td><td>1.628</td><td>1.000</td><td>1.000</td><td>0.500</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-lys-5_at</th><td>2.093</td><td>1.000</td><td>0.911</td><td>1.056</td><td>0.785</td><td>1.751</td><td>1.928</td><td>1.544</td><td>0.420</td><td>0.391</td><td>0.374</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-lys-M_at</th><td>2.915</td><td>0.883</td><td>0.726</td><td>0.971</td><td>1.000</td><td>3.782</td><td>3.008</td><td>1.557</td><td>1.028</td><td>0.924</td><td>0.421</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-phe-3_at</th><td>1.573</td><td>0.717</td><td>0.559</td><td>1.000</td><td>1.131</td><td>1.765</td><td>1.345</td><td>1.342</td><td>0.722</td><td>0.761</td><td>0.644</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-phe-5_at</th><td>2.136</td><td>0.889</td><td>0.558</td><td>1.049</td><td>1.241</td><td>2.136</td><td>1.000</td><td>2.090</td><td>0.679</td><td>0.888</td><td>0.559</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-phe-M_at</th><td>1.976</td><td>0.894</td><td>0.596</td><td>1.107</td><td>1.151</td><td>1.644</td><td>1.000</td><td>1.500</td><td>0.667</td><td>0.696</td><td>0.547</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-thr-3_s_at</th><td>1.635</td><td>1.000</td><td>0.852</td><td>1.000</td><td>1.068</td><td>1.342</td><td>1.211</td><td>1.788</td><td>0.552</td><td>0.372</td><td>0.570</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-thr-5_s_at</th><td>1.451</td><td>1.204</td><td>1.033</td><td>1.116</td><td>0.920</td><td>1.000</td><td>0.600</td><td>2.324</td><td>0.143</td><td>0.078</td><td>0.538</td></tr>
	<tr><th scope=row>AFFX-r2-Bs-thr-M_s_at</th><td>1.770</td><td>0.980</td><td>0.900</td><td>1.085</td><td>1.000</td><td>1.125</td><td>1.027</td><td>2.199</td><td>0.287</td><td>0.126</td><td>0.631</td></tr>
	<tr><th scope=row>AFFX-r2-Ec-bioB-3_at</th><td>1.000</td><td>0.858</td><td>1.156</td><td>0.768</td><td>0.852</td><td>2.324</td><td>3.168</td><td>1.150</td><td>0.986</td><td>1.015</td><td>0.654</td></tr>
	<tr><th scope=row>AFFX-r2-Ec-bioB-5_at</th><td>1.000</td><td>0.928</td><td>1.190</td><td>0.721</td><td>0.890</td><td>4.509</td><td>3.424</td><td>1.500</td><td>1.000</td><td>1.386</td><td>0.681</td></tr>
	<tr><th scope=row>AFFX-r2-Ec-bioB-M_at</th><td>1.000</td><td>0.826</td><td>1.049</td><td>0.672</td><td>0.708</td><td>2.928</td><td>3.175</td><td>1.304</td><td>1.000</td><td>1.194</td><td>0.607</td></tr>
	<tr><th scope=row>AFFX-r2-Ec-bioC-3_at</th><td>1.000</td><td>0.867</td><td>1.064</td><td>0.701</td><td>0.756</td><td>2.118</td><td>2.457</td><td>1.353</td><td>0.967</td><td>1.067</td><td>0.593</td></tr>
	<tr><th scope=row>AFFX-r2-Ec-bioC-5_at</th><td>0.988</td><td>0.897</td><td>1.000</td><td>0.690</td><td>0.730</td><td>1.761</td><td>1.739</td><td>1.149</td><td>1.016</td><td>1.165</td><td>0.622</td></tr>
	<tr><th scope=row>AFFX-r2-Ec-bioD-3_at</th><td>1.132</td><td>0.848</td><td>0.983</td><td>0.668</td><td>0.645</td><td>1.840</td><td>1.818</td><td>1.165</td><td>1.000</td><td>1.035</td><td>0.496</td></tr>
	<tr><th scope=row>AFFX-r2-Ec-bioD-5_at</th><td>1.171</td><td>0.819</td><td>0.794</td><td>0.647</td><td>0.617</td><td>1.548</td><td>1.337</td><td>1.000</td><td>1.055</td><td>1.086</td><td>0.510</td></tr>
	<tr><th scope=row>AFFX-r2-P1-cre-3_at</th><td>0.999</td><td>1.023</td><td>1.024</td><td>0.807</td><td>0.869</td><td>2.063</td><td>1.995</td><td>1.305</td><td>0.968</td><td>1.000</td><td>0.772</td></tr>
	<tr><th scope=row>AFFX-r2-P1-cre-5_at</th><td>0.994</td><td>1.006</td><td>1.035</td><td>0.802</td><td>0.780</td><td>2.183</td><td>2.175</td><td>1.331</td><td>0.962</td><td>1.000</td><td>0.698</td></tr>
	<tr><th scope=row>AFFX-ThrX-3_at</th><td>1.648</td><td>0.892</td><td>0.792</td><td>1.008</td><td>1.000</td><td>1.186</td><td>1.207</td><td>1.709</td><td>0.490</td><td>0.262</td><td>0.550</td></tr>
	<tr><th scope=row>AFFX-ThrX-5_at</th><td>1.168</td><td>1.109</td><td>1.058</td><td>1.000</td><td>0.792</td><td>1.073</td><td>0.650</td><td>2.098</td><td>0.167</td><td>0.099</td><td>0.480</td></tr>
	<tr><th scope=row>AFFX-ThrX-M_at</th><td>1.868</td><td>0.983</td><td>0.908</td><td>1.134</td><td>1.000</td><td>1.248</td><td>1.023</td><td>2.136</td><td>0.245</td><td>0.109</td><td>0.594</td></tr>
	<tr><th scope=row>AFFX-TrpnX-3_at</th><td>0.981</td><td>1.027</td><td>1.240</td><td>1.059</td><td>1.069</td><td>1.000</td><td>0.990</td><td>0.998</td><td>0.978</td><td>1.045</td><td>0.982</td></tr>
	<tr><th scope=row>AFFX-TrpnX-5_at</th><td>0.999</td><td>1.008</td><td>1.195</td><td>1.036</td><td>1.056</td><td>1.000</td><td>0.996</td><td>0.994</td><td>0.986</td><td>1.016</td><td>0.986</td></tr>
	<tr><th scope=row>AFFX-TrpnX-M_at</th><td>1.000</td><td>1.038</td><td>1.226</td><td>1.071</td><td>1.084</td><td>0.985</td><td>0.979</td><td>0.990</td><td>0.962</td><td>1.025</td><td>0.973</td></tr>
</tbody>
</table>




```R
ids <- AnnoProbe::idmap('GPL570',type = 'soft')
```

    file downloaded in /home/jfckkiu/AS_HG/Final_Results/Figure7_RNA_Validation_ALDOA
    



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
<caption>A data.frame: 54675 × 13</caption>
<thead>
	<tr><th scope=col>Probe</th><th scope=col>GSM1019539</th><th scope=col>GSM1019540</th><th scope=col>GSM1019541</th><th scope=col>GSM1019542</th><th scope=col>GSM1019543</th><th scope=col>GSM1019544</th><th scope=col>GSM1019545</th><th scope=col>GSM1019546</th><th scope=col>GSM1019547</th><th scope=col>GSM1019548</th><th scope=col>GSM1019549</th><th scope=col>symbol</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1007_s_at   </td><td>0.796</td><td>1.034</td><td>0.752</td><td>1.381</td><td>0.710</td><td>1.064</td><td>1.064</td><td>1.383</td><td>0.778</td><td>0.838</td><td>1.000</td><td>DDR1 /// MIR4640          </td></tr>
	<tr><td>1053_at     </td><td>1.000</td><td>1.004</td><td>2.045</td><td>1.259</td><td>1.024</td><td>0.904</td><td>0.823</td><td>0.992</td><td>0.954</td><td>0.450</td><td>1.231</td><td>RFC2                      </td></tr>
	<tr><td>117_at      </td><td>5.436</td><td>1.000</td><td>0.694</td><td>1.246</td><td>1.697</td><td>5.746</td><td>0.733</td><td>1.030</td><td>0.907</td><td>0.519</td><td>0.986</td><td>HSPA6                     </td></tr>
	<tr><td>121_at      </td><td>1.109</td><td>0.911</td><td>1.390</td><td>1.417</td><td>1.187</td><td>0.938</td><td>1.404</td><td>0.530</td><td>0.478</td><td>1.000</td><td>0.776</td><td>PAX8                      </td></tr>
	<tr><td>1255_g_at   </td><td>1.000</td><td>1.008</td><td>1.139</td><td>1.018</td><td>1.046</td><td>0.997</td><td>0.991</td><td>0.992</td><td>0.976</td><td>1.013</td><td>0.982</td><td>GUCA1A                    </td></tr>
	<tr><td>1294_at     </td><td>1.003</td><td>0.979</td><td>1.014</td><td>1.048</td><td>1.040</td><td>0.866</td><td>0.896</td><td>2.246</td><td>1.000</td><td>0.826</td><td>0.980</td><td>MIR5193 /// UBA7          </td></tr>
	<tr><td>1316_at     </td><td>1.365</td><td>1.108</td><td>0.778</td><td>1.000</td><td>0.985</td><td>0.965</td><td>2.252</td><td>0.822</td><td>1.161</td><td>1.189</td><td>0.607</td><td>THRA                      </td></tr>
	<tr><td>1320_at     </td><td>0.968</td><td>1.000</td><td>1.191</td><td>1.037</td><td>1.056</td><td>0.977</td><td>0.956</td><td>1.023</td><td>0.962</td><td>1.018</td><td>0.984</td><td>PTPN21                    </td></tr>
	<tr><td>1405_i_at   </td><td>1.112</td><td>2.209</td><td>1.038</td><td>1.000</td><td>0.560</td><td>0.626</td><td>3.283</td><td>0.584</td><td>0.826</td><td>0.482</td><td>1.194</td><td>CCL5                      </td></tr>
	<tr><td>1431_at     </td><td>0.898</td><td>0.942</td><td>1.206</td><td>0.981</td><td>1.000</td><td>1.021</td><td>0.919</td><td>1.010</td><td>1.001</td><td>0.956</td><td>1.127</td><td>CYP2E1                    </td></tr>
	<tr><td>1438_at     </td><td>0.990</td><td>1.000</td><td>1.102</td><td>1.011</td><td>1.023</td><td>0.976</td><td>0.990</td><td>1.012</td><td>0.981</td><td>1.008</td><td>0.988</td><td>EPHB3                     </td></tr>
	<tr><td>1487_at     </td><td>1.246</td><td>0.854</td><td>1.074</td><td>1.251</td><td>1.331</td><td>0.938</td><td>0.912</td><td>1.000</td><td>0.703</td><td>0.593</td><td>1.119</td><td>ESRRA                     </td></tr>
	<tr><td>1494_f_at   </td><td>0.989</td><td>1.000</td><td>1.117</td><td>1.033</td><td>1.025</td><td>0.997</td><td>0.997</td><td>1.005</td><td>0.975</td><td>1.032</td><td>0.998</td><td>CYP2A6                    </td></tr>
	<tr><td>1552256_a_at</td><td>1.100</td><td>1.100</td><td>1.089</td><td>1.339</td><td>1.000</td><td>0.496</td><td>0.206</td><td>0.362</td><td>0.163</td><td>0.172</td><td>1.026</td><td>SCARB1                    </td></tr>
	<tr><td>1552257_a_at</td><td>0.984</td><td>0.786</td><td>1.401</td><td>1.356</td><td>1.058</td><td>0.991</td><td>1.565</td><td>1.000</td><td>0.563</td><td>0.403</td><td>1.603</td><td>TTLL12                    </td></tr>
	<tr><td>1552258_at  </td><td>2.001</td><td>1.000</td><td>0.670</td><td>1.786</td><td>1.016</td><td>0.854</td><td>0.691</td><td>0.551</td><td>1.000</td><td>1.000</td><td>0.550</td><td>LINC00152 /// LOC101930489</td></tr>
	<tr><td>1552261_at  </td><td>0.988</td><td>1.036</td><td>1.061</td><td>1.133</td><td>1.000</td><td>0.970</td><td>0.904</td><td>1.207</td><td>0.884</td><td>1.002</td><td>0.966</td><td>WFDC2                     </td></tr>
	<tr><td>1552263_at  </td><td>0.413</td><td>1.220</td><td>1.700</td><td>1.480</td><td>1.628</td><td>0.783</td><td>1.098</td><td>0.445</td><td>1.000</td><td>1.000</td><td>0.861</td><td>MAPK1                     </td></tr>
	<tr><td>1552264_a_at</td><td>0.777</td><td>1.282</td><td>2.295</td><td>1.227</td><td>1.658</td><td>1.000</td><td>0.978</td><td>0.537</td><td>1.362</td><td>0.790</td><td>0.809</td><td>MAPK1                     </td></tr>
	<tr><td>1552266_at  </td><td>0.842</td><td>0.965</td><td>1.330</td><td>1.021</td><td>1.069</td><td>1.400</td><td>0.973</td><td>0.921</td><td>1.000</td><td>1.000</td><td>0.901</td><td>ADAM32                    </td></tr>
	<tr><td>1552269_at  </td><td>0.959</td><td>0.967</td><td>1.085</td><td>1.000</td><td>1.000</td><td>1.001</td><td>1.049</td><td>1.216</td><td>1.000</td><td>1.000</td><td>0.994</td><td>SPATA17                   </td></tr>
	<tr><td>1552271_at  </td><td>0.982</td><td>0.990</td><td>1.020</td><td>1.002</td><td>1.155</td><td>0.990</td><td>1.277</td><td>1.319</td><td>0.979</td><td>1.000</td><td>0.983</td><td>PRR22                     </td></tr>
	<tr><td>1552272_a_at</td><td>1.103</td><td>1.000</td><td>1.096</td><td>1.025</td><td>1.024</td><td>0.974</td><td>0.982</td><td>0.972</td><td>0.962</td><td>1.010</td><td>0.999</td><td>PRR22                     </td></tr>
	<tr><td>1552274_at  </td><td>1.000</td><td>0.861</td><td>0.781</td><td>0.608</td><td>0.896</td><td>1.286</td><td>1.797</td><td>1.332</td><td>1.443</td><td>3.070</td><td>0.652</td><td>PXK                       </td></tr>
	<tr><td>1552275_s_at</td><td>0.863</td><td>0.930</td><td>1.500</td><td>1.000</td><td>1.000</td><td>1.331</td><td>0.876</td><td>0.707</td><td>1.403</td><td>0.827</td><td>1.157</td><td>PXK                       </td></tr>
	<tr><td>1552276_a_at</td><td>0.994</td><td>1.030</td><td>1.136</td><td>1.089</td><td>1.000</td><td>0.978</td><td>0.883</td><td>1.447</td><td>0.876</td><td>1.015</td><td>0.961</td><td>VPS18                     </td></tr>
	<tr><td>1552277_a_at</td><td>1.313</td><td>1.000</td><td>0.552</td><td>1.186</td><td>0.632</td><td>0.901</td><td>2.144</td><td>0.404</td><td>1.597</td><td>2.330</td><td>0.436</td><td>MSANTD3                   </td></tr>
	<tr><td>1552278_a_at</td><td>1.000</td><td>1.745</td><td>4.224</td><td>1.008</td><td>3.187</td><td>0.888</td><td>0.729</td><td>0.606</td><td>0.673</td><td>0.520</td><td>1.111</td><td>SLC46A1                   </td></tr>
	<tr><td>1552279_a_at</td><td>1.007</td><td>1.000</td><td>1.393</td><td>1.042</td><td>2.194</td><td>0.940</td><td>0.834</td><td>0.816</td><td>0.567</td><td>0.726</td><td>1.361</td><td>SLC46A1                   </td></tr>
	<tr><td>1552280_at  </td><td>0.938</td><td>4.269</td><td>1.116</td><td>1.326</td><td>1.017</td><td>0.780</td><td>1.000</td><td>0.807</td><td>0.872</td><td>0.850</td><td>1.352</td><td>TIMD4                     </td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>AFFX-PheX-3_at       </td><td>1.650</td><td>0.667</td><td>0.630</td><td>1.000</td><td>1.256</td><td>1.799</td><td>1.356</td><td>1.409</td><td>0.748</td><td>0.908</td><td>0.625</td><td></td></tr>
	<tr><td>AFFX-PheX-5_at       </td><td>1.947</td><td>1.000</td><td>0.761</td><td>1.218</td><td>1.197</td><td>2.005</td><td>0.874</td><td>1.983</td><td>0.757</td><td>0.676</td><td>0.717</td><td></td></tr>
	<tr><td>AFFX-PheX-M_at       </td><td>1.843</td><td>0.750</td><td>0.517</td><td>1.000</td><td>1.020</td><td>1.705</td><td>1.022</td><td>1.461</td><td>0.631</td><td>0.709</td><td>0.505</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-3_at  </td><td>1.403</td><td>0.902</td><td>0.809</td><td>1.000</td><td>1.000</td><td>1.720</td><td>1.654</td><td>1.399</td><td>0.983</td><td>0.845</td><td>0.591</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-5_at  </td><td>1.683</td><td>1.075</td><td>1.000</td><td>0.715</td><td>0.661</td><td>2.056</td><td>2.106</td><td>2.362</td><td>0.584</td><td>0.406</td><td>0.226</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-M_at  </td><td>2.564</td><td>0.944</td><td>0.741</td><td>1.000</td><td>0.836</td><td>3.294</td><td>1.833</td><td>1.550</td><td>1.143</td><td>0.760</td><td>0.474</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-3_at  </td><td>2.193</td><td>0.740</td><td>0.523</td><td>0.896</td><td>0.774</td><td>3.148</td><td>2.830</td><td>1.628</td><td>1.000</td><td>1.000</td><td>0.500</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-5_at  </td><td>2.093</td><td>1.000</td><td>0.911</td><td>1.056</td><td>0.785</td><td>1.751</td><td>1.928</td><td>1.544</td><td>0.420</td><td>0.391</td><td>0.374</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-M_at  </td><td>2.915</td><td>0.883</td><td>0.726</td><td>0.971</td><td>1.000</td><td>3.782</td><td>3.008</td><td>1.557</td><td>1.028</td><td>0.924</td><td>0.421</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-3_at  </td><td>1.573</td><td>0.717</td><td>0.559</td><td>1.000</td><td>1.131</td><td>1.765</td><td>1.345</td><td>1.342</td><td>0.722</td><td>0.761</td><td>0.644</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-5_at  </td><td>2.136</td><td>0.889</td><td>0.558</td><td>1.049</td><td>1.241</td><td>2.136</td><td>1.000</td><td>2.090</td><td>0.679</td><td>0.888</td><td>0.559</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-M_at  </td><td>1.976</td><td>0.894</td><td>0.596</td><td>1.107</td><td>1.151</td><td>1.644</td><td>1.000</td><td>1.500</td><td>0.667</td><td>0.696</td><td>0.547</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-3_s_at</td><td>1.635</td><td>1.000</td><td>0.852</td><td>1.000</td><td>1.068</td><td>1.342</td><td>1.211</td><td>1.788</td><td>0.552</td><td>0.372</td><td>0.570</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-5_s_at</td><td>1.451</td><td>1.204</td><td>1.033</td><td>1.116</td><td>0.920</td><td>1.000</td><td>0.600</td><td>2.324</td><td>0.143</td><td>0.078</td><td>0.538</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-M_s_at</td><td>1.770</td><td>0.980</td><td>0.900</td><td>1.085</td><td>1.000</td><td>1.125</td><td>1.027</td><td>2.199</td><td>0.287</td><td>0.126</td><td>0.631</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-3_at </td><td>1.000</td><td>0.858</td><td>1.156</td><td>0.768</td><td>0.852</td><td>2.324</td><td>3.168</td><td>1.150</td><td>0.986</td><td>1.015</td><td>0.654</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-5_at </td><td>1.000</td><td>0.928</td><td>1.190</td><td>0.721</td><td>0.890</td><td>4.509</td><td>3.424</td><td>1.500</td><td>1.000</td><td>1.386</td><td>0.681</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-M_at </td><td>1.000</td><td>0.826</td><td>1.049</td><td>0.672</td><td>0.708</td><td>2.928</td><td>3.175</td><td>1.304</td><td>1.000</td><td>1.194</td><td>0.607</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioC-3_at </td><td>1.000</td><td>0.867</td><td>1.064</td><td>0.701</td><td>0.756</td><td>2.118</td><td>2.457</td><td>1.353</td><td>0.967</td><td>1.067</td><td>0.593</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioC-5_at </td><td>0.988</td><td>0.897</td><td>1.000</td><td>0.690</td><td>0.730</td><td>1.761</td><td>1.739</td><td>1.149</td><td>1.016</td><td>1.165</td><td>0.622</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioD-3_at </td><td>1.132</td><td>0.848</td><td>0.983</td><td>0.668</td><td>0.645</td><td>1.840</td><td>1.818</td><td>1.165</td><td>1.000</td><td>1.035</td><td>0.496</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioD-5_at </td><td>1.171</td><td>0.819</td><td>0.794</td><td>0.647</td><td>0.617</td><td>1.548</td><td>1.337</td><td>1.000</td><td>1.055</td><td>1.086</td><td>0.510</td><td></td></tr>
	<tr><td>AFFX-r2-P1-cre-3_at  </td><td>0.999</td><td>1.023</td><td>1.024</td><td>0.807</td><td>0.869</td><td>2.063</td><td>1.995</td><td>1.305</td><td>0.968</td><td>1.000</td><td>0.772</td><td></td></tr>
	<tr><td>AFFX-r2-P1-cre-5_at  </td><td>0.994</td><td>1.006</td><td>1.035</td><td>0.802</td><td>0.780</td><td>2.183</td><td>2.175</td><td>1.331</td><td>0.962</td><td>1.000</td><td>0.698</td><td></td></tr>
	<tr><td>AFFX-ThrX-3_at       </td><td>1.648</td><td>0.892</td><td>0.792</td><td>1.008</td><td>1.000</td><td>1.186</td><td>1.207</td><td>1.709</td><td>0.490</td><td>0.262</td><td>0.550</td><td></td></tr>
	<tr><td>AFFX-ThrX-5_at       </td><td>1.168</td><td>1.109</td><td>1.058</td><td>1.000</td><td>0.792</td><td>1.073</td><td>0.650</td><td>2.098</td><td>0.167</td><td>0.099</td><td>0.480</td><td></td></tr>
	<tr><td>AFFX-ThrX-M_at       </td><td>1.868</td><td>0.983</td><td>0.908</td><td>1.134</td><td>1.000</td><td>1.248</td><td>1.023</td><td>2.136</td><td>0.245</td><td>0.109</td><td>0.594</td><td></td></tr>
	<tr><td>AFFX-TrpnX-3_at      </td><td>0.981</td><td>1.027</td><td>1.240</td><td>1.059</td><td>1.069</td><td>1.000</td><td>0.990</td><td>0.998</td><td>0.978</td><td>1.045</td><td>0.982</td><td></td></tr>
	<tr><td>AFFX-TrpnX-5_at      </td><td>0.999</td><td>1.008</td><td>1.195</td><td>1.036</td><td>1.056</td><td>1.000</td><td>0.996</td><td>0.994</td><td>0.986</td><td>1.016</td><td>0.986</td><td></td></tr>
	<tr><td>AFFX-TrpnX-M_at      </td><td>1.000</td><td>1.038</td><td>1.226</td><td>1.071</td><td>1.084</td><td>0.985</td><td>0.979</td><td>0.990</td><td>0.962</td><td>1.025</td><td>0.973</td><td></td></tr>
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
<ol class=list-inline><li>'GSM1019539'</li><li>'GSM1019540'</li><li>'GSM1019541'</li><li>'GSM1019542'</li><li>'GSM1019543'</li><li>'GSM1019544'</li><li>'GSM1019545'</li><li>'GSM1019546'</li><li>'GSM1019547'</li><li>'GSM1019548'</li><li>'GSM1019549'</li></ol>




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


'GSM1019539=mean(GSM1019539),GSM1019540=mean(GSM1019540),GSM1019541=mean(GSM1019541),GSM1019542=mean(GSM1019542),GSM1019543=mean(GSM1019543),GSM1019544=mean(GSM1019544),GSM1019545=mean(GSM1019545),GSM1019546=mean(GSM1019546),GSM1019547=mean(GSM1019547),GSM1019548=mean(GSM1019548),GSM1019549=mean(GSM1019549),'



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
<caption>A data.frame: 6 × 12</caption>
<thead>
	<tr><th></th><th scope=col>GSM1019539</th><th scope=col>GSM1019540</th><th scope=col>GSM1019541</th><th scope=col>GSM1019542</th><th scope=col>GSM1019543</th><th scope=col>GSM1019544</th><th scope=col>GSM1019545</th><th scope=col>GSM1019546</th><th scope=col>GSM1019547</th><th scope=col>GSM1019548</th><th scope=col>GSM1019549</th><th scope=col>symbol</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>0.796</td><td>1.034</td><td>0.752</td><td>1.381</td><td>0.710</td><td>1.064</td><td>1.064</td><td>1.383</td><td>0.778</td><td>0.838</td><td>1.000</td><td>DDR1 /// MIR4640</td></tr>
	<tr><th scope=row>2</th><td>1.000</td><td>1.004</td><td>2.045</td><td>1.259</td><td>1.024</td><td>0.904</td><td>0.823</td><td>0.992</td><td>0.954</td><td>0.450</td><td>1.231</td><td>RFC2            </td></tr>
	<tr><th scope=row>3</th><td>5.436</td><td>1.000</td><td>0.694</td><td>1.246</td><td>1.697</td><td>5.746</td><td>0.733</td><td>1.030</td><td>0.907</td><td>0.519</td><td>0.986</td><td>HSPA6           </td></tr>
	<tr><th scope=row>4</th><td>1.109</td><td>0.911</td><td>1.390</td><td>1.417</td><td>1.187</td><td>0.938</td><td>1.404</td><td>0.530</td><td>0.478</td><td>1.000</td><td>0.776</td><td>PAX8            </td></tr>
	<tr><th scope=row>5</th><td>1.000</td><td>1.008</td><td>1.139</td><td>1.018</td><td>1.046</td><td>0.997</td><td>0.991</td><td>0.992</td><td>0.976</td><td>1.013</td><td>0.982</td><td>GUCA1A          </td></tr>
	<tr><th scope=row>6</th><td>1.003</td><td>0.979</td><td>1.014</td><td>1.048</td><td>1.040</td><td>0.866</td><td>0.896</td><td>2.246</td><td>1.000</td><td>0.826</td><td>0.980</td><td>MIR5193 /// UBA7</td></tr>
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
              summarise(GSM1019539=mean(GSM1019539),GSM1019540=mean(GSM1019540),GSM1019541=mean(GSM1019541),GSM1019542=mean(GSM1019542),GSM1019543=mean(GSM1019543),GSM1019544=mean(GSM1019544),GSM1019545=mean(GSM1019545),GSM1019546=mean(GSM1019546),GSM1019547=mean(GSM1019547),GSM1019548=mean(GSM1019548),GSM1019549=mean(GSM1019549))  %>% 
         column_to_rownames(var ='symbol')
```


```R
gset$GSE41571_series_matrix.txt.gz@phenoData@data  %>%  head() 
```


<table class="dataframe">
<caption>A data.frame: 6 × 47</caption>
<thead>
	<tr><th></th><th scope=col>title</th><th scope=col>geo_accession</th><th scope=col>status</th><th scope=col>submission_date</th><th scope=col>last_update_date</th><th scope=col>type</th><th scope=col>channel_count</th><th scope=col>source_name_ch1</th><th scope=col>organism_ch1</th><th scope=col>characteristics_ch1</th><th scope=col>⋯</th><th scope=col>contact_country</th><th scope=col>supplementary_file</th><th scope=col>data_row_count</th><th scope=col>age:ch1</th><th scope=col>cell type:ch1</th><th scope=col>days from last symptom:ch1</th><th scope=col>gender:ch1</th><th scope=col>plaque histology:ch1</th><th scope=col>radiological stenosis:ch1</th><th scope=col>symptom event:ch1</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM1019539</th><td>Ruptured Plaque, biological rep 1</td><td>GSM1019539</td><td>Public on Dec 31 2012</td><td>Oct 14 2012</td><td>Dec 31 2012</td><td>RNA</td><td>1</td><td>Macrophages, LMDed from Ruptured Carotid Atheromatous Plaque</td><td>Homo sapiens</td><td>age: 67</td><td>⋯</td><td>United Kingdom</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1019nnn/GSM1019539/suppl/GSM1019539_CEAr17.CEL.gz </td><td>54675</td><td>67</td><td>macrophage</td><td>30 </td><td>male  </td><td>Ruptured Thin Fibrous Cap Atheroma</td><td>90.0 (%)</td><td>CVA         </td></tr>
	<tr><th scope=row>GSM1019540</th><td>Ruptured Plaque, biological rep 2</td><td>GSM1019540</td><td>Public on Dec 31 2012</td><td>Oct 14 2012</td><td>Dec 31 2012</td><td>RNA</td><td>1</td><td>Macrophages, LMDed from Ruptured Carotid Atheromatous Plaque</td><td>Homo sapiens</td><td>age: 63</td><td>⋯</td><td>United Kingdom</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1019nnn/GSM1019540/suppl/GSM1019540_CEAr33.CEL.gz </td><td>54675</td><td>63</td><td>macrophage</td><td>3  </td><td>male  </td><td>Ruptured Thin Fibrous Cap Atheroma</td><td>70.0 (%)</td><td>AFugax      </td></tr>
	<tr><th scope=row>GSM1019541</th><td>Ruptured Plaque, biological rep 3</td><td>GSM1019541</td><td>Public on Dec 31 2012</td><td>Oct 14 2012</td><td>Dec 31 2012</td><td>RNA</td><td>1</td><td>Macrophages, LMDed from Ruptured Carotid Atheromatous Plaque</td><td>Homo sapiens</td><td>age: 53</td><td>⋯</td><td>United Kingdom</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1019nnn/GSM1019541/suppl/GSM1019541_CEAr35.CEL.gz </td><td>54675</td><td>53</td><td>macrophage</td><td>28 </td><td>male  </td><td>Ruptured Thin Fibrous Cap Atheroma</td><td>70.0 (%)</td><td>TIA         </td></tr>
	<tr><th scope=row>GSM1019542</th><td>Ruptured Plaque, biological rep 4</td><td>GSM1019542</td><td>Public on Dec 31 2012</td><td>Oct 14 2012</td><td>Dec 31 2012</td><td>RNA</td><td>1</td><td>Macrophages, LMDed from Ruptured Carotid Atheromatous Plaque</td><td>Homo sapiens</td><td>age: 74</td><td>⋯</td><td>United Kingdom</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1019nnn/GSM1019542/suppl/GSM1019542_CEAr57.CEL.gz </td><td>54675</td><td>74</td><td>macrophage</td><td>180</td><td>female</td><td>Ruptured Thin Fibrous Cap Atheroma</td><td>95.0 (%)</td><td>AFugax      </td></tr>
	<tr><th scope=row>GSM1019543</th><td>Ruptured Plaque, biological rep 5</td><td>GSM1019543</td><td>Public on Dec 31 2012</td><td>Oct 14 2012</td><td>Dec 31 2012</td><td>RNA</td><td>1</td><td>Macrophages, LMDed from Ruptured Carotid Atheromatous Plaque</td><td>Homo sapiens</td><td>age: 73</td><td>⋯</td><td>United Kingdom</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1019nnn/GSM1019543/suppl/GSM1019543_CEAr100.CEL.gz</td><td>54675</td><td>73</td><td>macrophage</td><td>35 </td><td>female</td><td>Ruptured Thin Fibrous Cap Atheroma</td><td>60.0 (%)</td><td>CVA         </td></tr>
	<tr><th scope=row>GSM1019544</th><td>Stable Plaque, biological rep 1  </td><td>GSM1019544</td><td>Public on Dec 31 2012</td><td>Oct 14 2012</td><td>Dec 31 2012</td><td>RNA</td><td>1</td><td>Macrophages, LMDed from Stable Carotid Atheromatous Plaque  </td><td>Homo sapiens</td><td>age: 45</td><td>⋯</td><td>United Kingdom</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1019nnn/GSM1019544/suppl/GSM1019544_CEAs14.CEL.gz </td><td>54675</td><td>45</td><td>macrophage</td><td>365</td><td>female</td><td>Stable Thick Fibrous Cap Atheroma </td><td>85.0 (%)</td><td>CVA(lacunar)</td></tr>
</tbody>
</table>




```R
group_list <- c(rep(c('RUP','STA'),c(5,6)))
```


```R
group_list <- factor(group_list, levels = c('STA','RUP'))
```


```R
library(ggpubr)
```


```R
selected_gene <- c('ALDOA','CREG1','LGMN','PKM2')
```


```R
normalized_gset_mean[selected_gene,]
```


<table class="dataframe">
<caption>A data.frame: 4 × 11</caption>
<thead>
	<tr><th></th><th scope=col>GSM1019539</th><th scope=col>GSM1019540</th><th scope=col>GSM1019541</th><th scope=col>GSM1019542</th><th scope=col>GSM1019543</th><th scope=col>GSM1019544</th><th scope=col>GSM1019545</th><th scope=col>GSM1019546</th><th scope=col>GSM1019547</th><th scope=col>GSM1019548</th><th scope=col>GSM1019549</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ALDOA</th><td>1.231333</td><td>1.262</td><td>1.321</td><td>1.401</td><td>1.000</td><td>0.435</td><td>1.252333</td><td>0.723</td><td>0.6603333</td><td>0.770</td><td>0.7426667</td></tr>
	<tr><th scope=row>CREG1</th><td>1.000000</td><td>1.560</td><td>1.365</td><td>1.456</td><td>1.398</td><td>0.727</td><td>0.660000</td><td>0.663</td><td>0.8000000</td><td>0.742</td><td>1.3020000</td></tr>
	<tr><th scope=row>LGMN</th><td>1.000000</td><td>1.165</td><td>1.183</td><td>1.117</td><td>1.297</td><td>0.911</td><td>0.347000</td><td>0.289</td><td>0.6380000</td><td>0.465</td><td>1.0770000</td></tr>
	<tr><th scope=row>NA</th><td>      NA</td><td>   NA</td><td>   NA</td><td>   NA</td><td>   NA</td><td>   NA</td><td>      NA</td><td>   NA</td><td>       NA</td><td>   NA</td><td>       NA</td></tr>
</tbody>
</table>




```R
normalized_gset_mean[selected_gene,]  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)
```


<table class="dataframe">
<caption>A data.frame: 11 × 6</caption>
<thead>
	<tr><th scope=col>sample</th><th scope=col>ALDOA</th><th scope=col>CREG1</th><th scope=col>LGMN</th><th scope=col>PKM</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM1019539</td><td>1.2313333</td><td>1.000</td><td>1.000</td><td>1.9900</td><td>RUP</td></tr>
	<tr><td>GSM1019540</td><td>1.2620000</td><td>1.560</td><td>1.165</td><td>1.0000</td><td>RUP</td></tr>
	<tr><td>GSM1019541</td><td>1.3210000</td><td>1.365</td><td>1.183</td><td>1.1155</td><td>RUP</td></tr>
	<tr><td>GSM1019542</td><td>1.4010000</td><td>1.456</td><td>1.117</td><td>1.1365</td><td>RUP</td></tr>
	<tr><td>GSM1019543</td><td>1.0000000</td><td>1.398</td><td>1.297</td><td>0.5615</td><td>RUP</td></tr>
	<tr><td>GSM1019544</td><td>0.4350000</td><td>0.727</td><td>0.911</td><td>0.5855</td><td>STA</td></tr>
	<tr><td>GSM1019545</td><td>1.2523333</td><td>0.660</td><td>0.347</td><td>1.2705</td><td>STA</td></tr>
	<tr><td>GSM1019546</td><td>0.7230000</td><td>0.663</td><td>0.289</td><td>0.6915</td><td>STA</td></tr>
	<tr><td>GSM1019547</td><td>0.6603333</td><td>0.800</td><td>0.638</td><td>1.1580</td><td>STA</td></tr>
	<tr><td>GSM1019548</td><td>0.7700000</td><td>0.742</td><td>0.465</td><td>0.9420</td><td>STA</td></tr>
	<tr><td>GSM1019549</td><td>0.7426667</td><td>1.302</td><td>1.077</td><td>0.1895</td><td>STA</td></tr>
</tbody>
</table>




```R
normalized_gset_mean[selected_gene,]  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       pivot_longer(-c(sample,group),names_to = "gene",values_to = "value")  %>% 
        mutate(gene = factor(gene, levels = c('ALDOA','CREG1','LGMN','PKM'))) %>% 
         ggboxplot(x = 'group', y = 'value',add = 'point',fill = 'group',legend = 'right',facet.by = 'gene',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('STA','RUP')),method = 't.test', paired = F,label = 'p.signif') 
```


    
![png](Step11.3_GSE41571_files/Step11.3_GSE41571_22_0.png)
    



```R
library(patchwork)
```


```R
p1 <- normalized_gset_mean[selected_gene,]  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       pivot_longer(-c(sample,group),names_to = "gene",values_to = "value")  %>% 
        mutate(gene = factor(gene, levels = c('ALDOA','CREG1','LGMN','PKM'))) %>% 
         filter(gene %in% 'ALDOA')%>% 
         ggboxplot(x = 'group', y = 'value', xlab = '',ylab = 'expression value',title = 'ALDOA',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('STA','RUP')),method = 't.test', paired = F,label = 'p.signif')  + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6),axis.title.y = element_text(size=6), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p1
p2 <- normalized_gset_mean[selected_gene,]  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       pivot_longer(-c(sample,group),names_to = "gene",values_to = "value")  %>% 
        mutate(gene = factor(gene, levels = c('ALDOA','CREG1','LGMN','PKM'))) %>% 
         filter(gene %in% 'CREG1')%>% 
         ggboxplot(x = 'group', y = 'value', xlab = '',ylab = 'expression value',title = 'CREG1',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('STA','RUP')),method = 't.test', paired = F,label = 'p.signif')  + 
           theme(axis.text.x = element_text(size=6),axis.text.y = element_text(size=6), axis.title.y = element_text(size=0), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p2
p3 <- normalized_gset_mean[selected_gene,]  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       pivot_longer(-c(sample,group),names_to = "gene",values_to = "value")  %>% 
        mutate(gene = factor(gene, levels = c('ALDOA','CREG1','LGMN','PKM'))) %>% 
         filter(gene %in% 'LGMN')%>% 
         ggboxplot(x = 'group', y = 'value', xlab = '',ylab = 'expression value',title = 'LGMN',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('STA','RUP')),method = 't.test', paired = F,label = 'p.signif')  + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6),axis.title.y = element_text(size=0), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p3
p4 <- normalized_gset_mean[selected_gene,]  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       pivot_longer(-c(sample,group),names_to = "gene",values_to = "value")  %>% 
        mutate(gene = factor(gene, levels = c('ALDOA','CREG1','LGMN','PKM'))) %>% 
         filter(gene %in% 'PKM')%>% 
         ggboxplot(x = 'group', y = 'value', xlab = '',ylab = 'expression value',title = 'PKM',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('STA','RUP')),method = 't.test', paired = F,label = 'p.signif')  + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), axis.title.y = element_text(size=0), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p4
```


    
![png](Step11.3_GSE41571_files/Step11.3_GSE41571_24_0.png)
    



    
![png](Step11.3_GSE41571_files/Step11.3_GSE41571_24_1.png)
    



    
![png](Step11.3_GSE41571_files/Step11.3_GSE41571_24_2.png)
    



    
![png](Step11.3_GSE41571_files/Step11.3_GSE41571_24_3.png)
    



```R
p0 <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 4,guides='collect')  +  
    plot_annotation(title = 'GSE41571',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE41571_4Gene.pdf', height = 5, width = 18, units = 'cm')
```


    
![png](Step11.3_GSE41571_files/Step11.3_GSE41571_25_0.png)
    

