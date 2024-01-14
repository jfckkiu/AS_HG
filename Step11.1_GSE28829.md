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
gset  <-  GEOquery::getGEO('GSE28829',getGPL = F)
```

    Found 1 file(s)
    
    GSE28829_series_matrix.txt.gz
    
    [1mRows: [22m[34m54675[39m [1mColumns: [22m[34m30[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m "\t"
    [31mchr[39m  (1): ID_REF
    [32mdbl[39m (29): GSM714070, GSM714071, GSM714072, GSM714073, GSM714074, GSM714075, ...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
normalized_gset  <- gset$GSE28829_series_matrix.txt.gz@assayData$exprs
```


```R
range(normalized_gset)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>2.4977043494</li><li>14.1413126635</li></ol>




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
<caption>A data.frame: 54675 × 31</caption>
<thead>
	<tr><th scope=col>Probe</th><th scope=col>GSM714070</th><th scope=col>GSM714071</th><th scope=col>GSM714072</th><th scope=col>GSM714073</th><th scope=col>GSM714074</th><th scope=col>GSM714075</th><th scope=col>GSM714076</th><th scope=col>GSM714077</th><th scope=col>GSM714078</th><th scope=col>⋯</th><th scope=col>GSM714090</th><th scope=col>GSM714091</th><th scope=col>GSM714092</th><th scope=col>GSM714093</th><th scope=col>GSM714094</th><th scope=col>GSM714095</th><th scope=col>GSM714096</th><th scope=col>GSM714097</th><th scope=col>GSM714098</th><th scope=col>symbol</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1007_s_at   </td><td>7.476734</td><td>7.560470</td><td>7.545101</td><td>7.137021</td><td>7.313099</td><td>7.320076</td><td>6.757667</td><td>7.379789</td><td>7.545558</td><td>⋯</td><td>7.427269</td><td>7.046685</td><td>7.362381</td><td>7.537203</td><td>7.467394</td><td>7.281994</td><td>7.340640</td><td>7.311836</td><td>7.066598</td><td>DDR1 /// MIR4640          </td></tr>
	<tr><td>1053_at     </td><td>5.231710</td><td>5.571284</td><td>5.502338</td><td>5.406858</td><td>5.961533</td><td>5.719497</td><td>5.464369</td><td>5.687824</td><td>5.465359</td><td>⋯</td><td>5.662993</td><td>5.803273</td><td>5.406180</td><td>5.713884</td><td>5.389508</td><td>5.171387</td><td>5.768445</td><td>5.735633</td><td>5.617395</td><td>RFC2                      </td></tr>
	<tr><td>117_at      </td><td>6.407551</td><td>6.456310</td><td>5.912837</td><td>6.138684</td><td>6.431826</td><td>6.274274</td><td>5.657739</td><td>5.846036</td><td>5.710486</td><td>⋯</td><td>5.860925</td><td>6.435217</td><td>5.536453</td><td>6.086993</td><td>5.458502</td><td>5.751894</td><td>5.863148</td><td>4.882628</td><td>5.738776</td><td>HSPA6                     </td></tr>
	<tr><td>121_at      </td><td>8.086023</td><td>8.056792</td><td>7.715602</td><td>7.423653</td><td>7.747922</td><td>7.833796</td><td>7.633695</td><td>7.562925</td><td>7.985645</td><td>⋯</td><td>7.660248</td><td>7.200921</td><td>8.062172</td><td>7.596732</td><td>7.435648</td><td>7.382841</td><td>7.493919</td><td>7.492598</td><td>7.106533</td><td>PAX8                      </td></tr>
	<tr><td>1255_g_at   </td><td>3.334447</td><td>3.174585</td><td>3.156232</td><td>3.277893</td><td>3.197367</td><td>3.195213</td><td>3.268658</td><td>3.442809</td><td>3.269492</td><td>⋯</td><td>3.269480</td><td>3.368173</td><td>3.337614</td><td>3.249723</td><td>3.268982</td><td>3.363065</td><td>3.181215</td><td>3.228107</td><td>3.360925</td><td>GUCA1A                    </td></tr>
	<tr><td>1294_at     </td><td>7.727043</td><td>7.573688</td><td>7.457024</td><td>7.626756</td><td>7.616204</td><td>7.706446</td><td>7.274741</td><td>7.540043</td><td>7.495845</td><td>⋯</td><td>7.439052</td><td>7.486536</td><td>8.113033</td><td>7.569026</td><td>7.275367</td><td>7.250188</td><td>7.376961</td><td>7.086117</td><td>7.022475</td><td>MIR5193 /// UBA7          </td></tr>
	<tr><td>1316_at     </td><td>5.922763</td><td>5.531696</td><td>5.614733</td><td>5.445089</td><td>5.669549</td><td>5.673246</td><td>5.566560</td><td>5.609161</td><td>5.890963</td><td>⋯</td><td>5.718211</td><td>5.659084</td><td>5.732430</td><td>5.768923</td><td>5.723051</td><td>5.833857</td><td>5.424068</td><td>5.509707</td><td>5.550262</td><td>THRA                      </td></tr>
	<tr><td>1320_at     </td><td>4.925586</td><td>4.779562</td><td>4.759425</td><td>4.761010</td><td>4.821824</td><td>4.692128</td><td>4.890008</td><td>4.874585</td><td>4.668920</td><td>⋯</td><td>5.682762</td><td>5.110186</td><td>5.165035</td><td>5.050457</td><td>5.065661</td><td>4.937367</td><td>5.007145</td><td>4.861808</td><td>4.894809</td><td>PTPN21                    </td></tr>
	<tr><td>1405_i_at   </td><td>6.226800</td><td>6.568606</td><td>6.474123</td><td>6.184142</td><td>6.355437</td><td>6.345883</td><td>5.339767</td><td>7.343654</td><td>6.941797</td><td>⋯</td><td>5.206461</td><td>5.166857</td><td>5.770775</td><td>6.233021</td><td>4.750236</td><td>5.036030</td><td>7.693500</td><td>4.446433</td><td>4.709236</td><td>CCL5                      </td></tr>
	<tr><td>1431_at     </td><td>3.393143</td><td>3.404954</td><td>3.488778</td><td>3.379285</td><td>3.419052</td><td>3.407246</td><td>3.385841</td><td>3.399413</td><td>3.529619</td><td>⋯</td><td>3.490949</td><td>3.592954</td><td>3.491016</td><td>3.390710</td><td>3.500920</td><td>3.470048</td><td>3.291529</td><td>3.333013</td><td>3.440773</td><td>CYP2E1                    </td></tr>
	<tr><td>1438_at     </td><td>5.382629</td><td>5.210232</td><td>5.135249</td><td>5.069951</td><td>5.400649</td><td>5.397060</td><td>5.333609</td><td>5.088325</td><td>5.215748</td><td>⋯</td><td>5.089140</td><td>5.232828</td><td>5.536511</td><td>5.334206</td><td>5.170820</td><td>5.033012</td><td>5.218836</td><td>5.159302</td><td>5.148694</td><td>EPHB3                     </td></tr>
	<tr><td>1487_at     </td><td>6.253220</td><td>6.439284</td><td>6.420131</td><td>6.562199</td><td>6.470325</td><td>6.486059</td><td>6.109339</td><td>6.477231</td><td>6.508121</td><td>⋯</td><td>6.329902</td><td>6.493921</td><td>6.070508</td><td>6.521998</td><td>6.322114</td><td>6.292051</td><td>6.301568</td><td>6.337027</td><td>6.629769</td><td>ESRRA                     </td></tr>
	<tr><td>1494_f_at   </td><td>5.312670</td><td>5.343201</td><td>5.253478</td><td>5.103065</td><td>5.004329</td><td>5.261961</td><td>5.218809</td><td>5.236136</td><td>5.506457</td><td>⋯</td><td>5.075365</td><td>5.328564</td><td>6.138902</td><td>5.173727</td><td>5.300710</td><td>5.142706</td><td>5.020928</td><td>5.124512</td><td>5.254878</td><td>CYP2A6                    </td></tr>
	<tr><td>1552256_a_at</td><td>7.676968</td><td>7.704522</td><td>7.831077</td><td>7.647633</td><td>7.617398</td><td>7.348637</td><td>7.119217</td><td>7.372458</td><td>7.256456</td><td>⋯</td><td>6.805328</td><td>6.458248</td><td>6.870781</td><td>7.503047</td><td>7.066669</td><td>6.190748</td><td>7.366874</td><td>6.596461</td><td>6.849367</td><td>SCARB1                    </td></tr>
	<tr><td>1552257_a_at</td><td>6.090183</td><td>6.083343</td><td>6.601154</td><td>6.498408</td><td>6.170246</td><td>6.169182</td><td>6.311429</td><td>6.078742</td><td>6.319092</td><td>⋯</td><td>6.121865</td><td>6.077653</td><td>6.939850</td><td>5.985401</td><td>6.497193</td><td>6.138559</td><td>5.851820</td><td>6.187926</td><td>6.384868</td><td>TTLL12                    </td></tr>
	<tr><td>1552258_at  </td><td>3.658457</td><td>3.771177</td><td>3.750931</td><td>3.806265</td><td>3.852486</td><td>3.528730</td><td>3.663895</td><td>3.781000</td><td>3.624378</td><td>⋯</td><td>3.877058</td><td>3.822665</td><td>4.415067</td><td>3.750931</td><td>3.727187</td><td>3.800974</td><td>3.750931</td><td>3.849578</td><td>3.728923</td><td>LINC00152 /// LOC101930489</td></tr>
	<tr><td>1552261_at  </td><td>4.871800</td><td>4.461214</td><td>4.508945</td><td>4.754654</td><td>4.529026</td><td>4.397758</td><td>4.976930</td><td>4.859077</td><td>4.946909</td><td>⋯</td><td>4.647098</td><td>5.512455</td><td>4.857301</td><td>4.658151</td><td>4.297735</td><td>4.625253</td><td>4.243038</td><td>4.718078</td><td>4.891017</td><td>WFDC2                     </td></tr>
	<tr><td>1552263_at  </td><td>5.066076</td><td>5.839615</td><td>5.373385</td><td>5.030533</td><td>5.312583</td><td>5.323242</td><td>4.596259</td><td>4.758841</td><td>5.047684</td><td>⋯</td><td>5.290104</td><td>4.074749</td><td>4.875813</td><td>5.424100</td><td>5.085033</td><td>4.984591</td><td>5.444044</td><td>4.560666</td><td>4.478742</td><td>MAPK1                     </td></tr>
	<tr><td>1552264_a_at</td><td>6.820752</td><td>6.897460</td><td>7.059128</td><td>6.923568</td><td>6.906530</td><td>6.942615</td><td>6.249993</td><td>6.695527</td><td>6.860637</td><td>⋯</td><td>7.034056</td><td>7.036639</td><td>7.151766</td><td>7.223901</td><td>6.367623</td><td>6.449209</td><td>7.198109</td><td>6.643571</td><td>7.127962</td><td>MAPK1                     </td></tr>
	<tr><td>1552266_at  </td><td>3.808949</td><td>3.938168</td><td>3.566033</td><td>4.069212</td><td>3.907701</td><td>4.035658</td><td>3.905610</td><td>3.902500</td><td>3.980783</td><td>⋯</td><td>3.589107</td><td>4.211306</td><td>3.935823</td><td>3.894312</td><td>4.212896</td><td>3.911250</td><td>4.033726</td><td>3.765515</td><td>4.025326</td><td>ADAM32                    </td></tr>
	<tr><td>1552269_at  </td><td>3.315224</td><td>3.351903</td><td>3.343055</td><td>3.176276</td><td>3.440586</td><td>3.279698</td><td>3.264766</td><td>3.289414</td><td>3.426627</td><td>⋯</td><td>3.328380</td><td>3.393927</td><td>3.384863</td><td>3.466472</td><td>3.439239</td><td>3.483471</td><td>3.522193</td><td>3.295778</td><td>3.437778</td><td>SPATA17                   </td></tr>
	<tr><td>1552271_at  </td><td>4.998514</td><td>5.021860</td><td>4.698744</td><td>4.663108</td><td>5.040071</td><td>5.041507</td><td>4.780876</td><td>4.867960</td><td>5.185665</td><td>⋯</td><td>4.771040</td><td>5.651885</td><td>5.212445</td><td>4.998514</td><td>5.169074</td><td>5.119749</td><td>4.840627</td><td>4.675780</td><td>5.074544</td><td>PRR22                     </td></tr>
	<tr><td>1552272_a_at</td><td>4.502365</td><td>4.514138</td><td>4.396829</td><td>4.594947</td><td>4.514138</td><td>4.705061</td><td>4.471497</td><td>4.514138</td><td>4.541600</td><td>⋯</td><td>4.329922</td><td>5.107246</td><td>5.215447</td><td>4.490106</td><td>4.851125</td><td>4.276867</td><td>4.383185</td><td>4.514138</td><td>4.587427</td><td>PRR22                     </td></tr>
	<tr><td>1552274_at  </td><td>6.706973</td><td>6.962160</td><td>6.896566</td><td>7.134047</td><td>6.936095</td><td>6.945230</td><td>6.663959</td><td>6.469898</td><td>6.826577</td><td>⋯</td><td>6.805328</td><td>5.467232</td><td>5.849547</td><td>6.761226</td><td>6.948154</td><td>6.210566</td><td>7.239242</td><td>6.681639</td><td>6.246376</td><td>PXK                       </td></tr>
	<tr><td>1552275_s_at</td><td>5.525243</td><td>5.729755</td><td>5.878115</td><td>5.703896</td><td>5.672967</td><td>5.459811</td><td>5.362306</td><td>5.745094</td><td>5.880654</td><td>⋯</td><td>5.748506</td><td>4.896907</td><td>5.315977</td><td>6.469112</td><td>5.502112</td><td>5.349649</td><td>6.328668</td><td>5.629886</td><td>5.305351</td><td>PXK                       </td></tr>
	<tr><td>1552276_a_at</td><td>5.262259</td><td>4.796568</td><td>4.820012</td><td>4.988796</td><td>4.653722</td><td>4.848155</td><td>4.816378</td><td>5.056708</td><td>4.997640</td><td>⋯</td><td>4.670812</td><td>5.237672</td><td>5.551484</td><td>4.793710</td><td>4.652157</td><td>4.943135</td><td>4.661699</td><td>5.144857</td><td>5.061299</td><td>VPS18                     </td></tr>
	<tr><td>1552277_a_at</td><td>7.109956</td><td>7.096568</td><td>7.004540</td><td>6.629028</td><td>7.096687</td><td>6.939557</td><td>6.147171</td><td>7.195497</td><td>7.017300</td><td>⋯</td><td>7.049342</td><td>5.989069</td><td>6.608316</td><td>7.325639</td><td>6.577773</td><td>6.276648</td><td>6.980042</td><td>6.554436</td><td>6.496104</td><td>MSANTD3                   </td></tr>
	<tr><td>1552278_a_at</td><td>3.977090</td><td>4.325501</td><td>4.129754</td><td>4.089251</td><td>4.402472</td><td>4.299672</td><td>3.998264</td><td>4.198220</td><td>4.478939</td><td>⋯</td><td>4.073133</td><td>4.569323</td><td>4.612156</td><td>4.212502</td><td>3.985845</td><td>4.137673</td><td>4.198220</td><td>3.972445</td><td>3.861258</td><td>SLC46A1                   </td></tr>
	<tr><td>1552279_a_at</td><td>6.001133</td><td>6.239700</td><td>6.367057</td><td>6.062821</td><td>6.380545</td><td>6.317848</td><td>5.935090</td><td>6.227399</td><td>6.442183</td><td>⋯</td><td>6.231542</td><td>6.187499</td><td>6.777705</td><td>6.289393</td><td>6.146336</td><td>5.866330</td><td>5.976506</td><td>6.027865</td><td>6.445722</td><td>SLC46A1                   </td></tr>
	<tr><td>1552280_at  </td><td>3.731025</td><td>3.623076</td><td>3.693623</td><td>3.718015</td><td>4.616735</td><td>4.410913</td><td>3.899488</td><td>3.685344</td><td>3.697639</td><td>⋯</td><td>3.536507</td><td>3.777249</td><td>3.532396</td><td>3.610351</td><td>3.750931</td><td>3.943110</td><td>3.453285</td><td>3.896961</td><td>3.748689</td><td>TIMD4                     </td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>AFFX-PheX-3_at       </td><td> 5.799662</td><td> 5.361233</td><td> 5.279383</td><td> 5.552708</td><td> 4.969027</td><td> 5.343321</td><td> 7.809518</td><td> 5.631090</td><td> 5.534106</td><td>⋯</td><td> 6.816246</td><td> 5.848892</td><td> 6.308399</td><td> 5.094977</td><td> 5.226444</td><td> 5.610618</td><td> 6.835850</td><td> 5.548966</td><td> 5.533085</td><td></td></tr>
	<tr><td>AFFX-PheX-5_at       </td><td> 3.175364</td><td> 3.273273</td><td> 3.199371</td><td> 3.207001</td><td> 3.136714</td><td> 3.320854</td><td> 6.369223</td><td> 3.407609</td><td> 3.285803</td><td>⋯</td><td> 6.121487</td><td> 3.620622</td><td> 5.363175</td><td> 3.137380</td><td> 3.245167</td><td> 3.301749</td><td> 6.164385</td><td> 3.264900</td><td> 3.358127</td><td></td></tr>
	<tr><td>AFFX-PheX-M_at       </td><td> 3.258624</td><td> 3.285021</td><td> 3.181127</td><td> 3.141482</td><td> 3.052829</td><td> 3.148094</td><td> 6.147589</td><td> 3.175217</td><td> 3.309512</td><td>⋯</td><td> 6.175647</td><td> 3.319391</td><td> 5.675294</td><td> 3.123700</td><td> 3.249791</td><td> 3.283854</td><td> 6.120446</td><td> 3.257231</td><td> 3.220881</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-3_at  </td><td> 3.211913</td><td> 3.059706</td><td> 2.985300</td><td> 2.959689</td><td> 3.049319</td><td> 3.105817</td><td>10.023311</td><td> 3.003516</td><td> 3.056608</td><td>⋯</td><td> 9.894703</td><td> 3.321225</td><td> 8.904542</td><td> 3.026063</td><td> 3.174229</td><td> 3.004560</td><td> 9.713705</td><td> 2.928985</td><td> 3.056093</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-5_at  </td><td> 3.162541</td><td> 3.113690</td><td> 3.027964</td><td> 3.240127</td><td> 3.054281</td><td> 3.353850</td><td> 9.517573</td><td> 3.179267</td><td> 3.340212</td><td>⋯</td><td> 8.860095</td><td> 3.587838</td><td> 7.604186</td><td> 3.195252</td><td> 3.197339</td><td> 2.928807</td><td> 8.997940</td><td> 3.011427</td><td> 3.294881</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-M_at  </td><td> 3.157598</td><td> 2.890638</td><td> 3.136778</td><td> 2.901710</td><td> 2.929485</td><td> 2.963133</td><td> 9.817492</td><td> 3.069762</td><td> 3.277522</td><td>⋯</td><td> 9.329845</td><td> 3.188204</td><td> 8.138690</td><td> 2.965548</td><td> 3.029518</td><td> 2.946767</td><td> 9.418128</td><td> 3.241891</td><td> 3.122422</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-3_at  </td><td> 3.178294</td><td> 3.168507</td><td> 3.154042</td><td> 3.173506</td><td> 2.990508</td><td> 3.106432</td><td> 4.776025</td><td> 3.263392</td><td> 3.302927</td><td>⋯</td><td> 5.649178</td><td> 3.886967</td><td> 3.908977</td><td> 3.022060</td><td> 3.111451</td><td> 3.318567</td><td> 4.461615</td><td> 3.370078</td><td> 3.316552</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-5_at  </td><td> 3.088552</td><td> 2.986761</td><td> 2.985280</td><td> 3.113196</td><td> 2.943532</td><td> 3.131835</td><td> 5.455542</td><td> 3.179267</td><td> 3.182525</td><td>⋯</td><td> 5.814741</td><td> 3.521736</td><td> 4.057161</td><td> 3.140203</td><td> 2.966510</td><td> 3.313245</td><td> 5.983796</td><td> 3.251478</td><td> 3.308964</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-M_at  </td><td> 3.370628</td><td> 3.337118</td><td> 3.112582</td><td> 3.221729</td><td> 3.090680</td><td> 3.002346</td><td> 5.649610</td><td> 3.249833</td><td> 3.228399</td><td>⋯</td><td> 5.490851</td><td> 3.450268</td><td> 5.050364</td><td> 3.044661</td><td> 3.133026</td><td> 3.096401</td><td> 5.616145</td><td> 3.485281</td><td> 3.158677</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-3_at  </td><td> 3.643435</td><td> 3.515682</td><td> 3.489131</td><td> 3.306345</td><td> 3.429928</td><td> 3.507174</td><td> 6.941771</td><td> 3.859380</td><td> 3.595655</td><td>⋯</td><td> 6.647299</td><td> 3.788878</td><td> 4.707560</td><td> 3.377882</td><td> 3.219862</td><td> 4.018056</td><td> 6.460470</td><td> 3.521919</td><td> 3.513904</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-5_at  </td><td> 3.208177</td><td> 2.900271</td><td> 2.833011</td><td> 2.980726</td><td> 2.942104</td><td> 2.978683</td><td> 6.163290</td><td> 3.215535</td><td> 2.965791</td><td>⋯</td><td> 6.946852</td><td> 3.331943</td><td> 5.571226</td><td> 2.916105</td><td> 2.957949</td><td> 3.086792</td><td> 6.980666</td><td> 2.991682</td><td> 3.014720</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-M_at  </td><td> 3.246981</td><td> 3.197366</td><td> 3.077441</td><td> 3.308927</td><td> 2.984028</td><td> 3.171773</td><td> 6.566068</td><td> 3.127552</td><td> 3.318979</td><td>⋯</td><td> 6.421281</td><td> 3.807055</td><td> 6.011063</td><td> 3.161144</td><td> 3.018333</td><td> 3.227569</td><td> 6.673568</td><td> 3.062022</td><td> 3.216779</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-3_s_at</td><td> 3.749516</td><td> 3.831985</td><td> 3.546602</td><td> 3.703074</td><td> 3.323028</td><td> 3.604053</td><td> 9.011896</td><td> 3.670496</td><td> 3.598536</td><td>⋯</td><td> 8.102188</td><td> 4.203371</td><td> 7.438782</td><td> 3.458968</td><td> 3.674070</td><td> 3.545755</td><td> 8.231378</td><td> 3.817005</td><td> 3.853144</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-5_s_at</td><td> 3.738782</td><td> 3.494834</td><td> 3.363877</td><td> 3.489904</td><td> 3.479617</td><td> 3.717343</td><td> 8.441673</td><td> 3.770833</td><td> 3.717343</td><td>⋯</td><td> 7.898420</td><td> 3.949770</td><td> 6.769996</td><td> 3.590111</td><td> 3.689974</td><td> 3.734969</td><td> 7.933559</td><td> 3.552002</td><td> 3.614822</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-M_s_at</td><td> 3.842585</td><td> 3.893437</td><td> 3.569990</td><td> 3.400236</td><td> 3.506575</td><td> 3.817073</td><td> 8.480940</td><td> 3.722567</td><td> 3.817073</td><td>⋯</td><td> 7.765366</td><td> 3.849880</td><td> 6.840218</td><td> 3.759086</td><td> 3.899208</td><td> 3.484919</td><td> 7.939149</td><td> 3.405452</td><td> 3.844785</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-3_at </td><td> 7.767545</td><td> 7.596855</td><td> 7.807600</td><td> 7.351109</td><td> 7.681252</td><td> 7.809094</td><td> 8.296721</td><td> 8.406608</td><td> 7.853791</td><td>⋯</td><td> 7.911205</td><td> 7.574871</td><td> 8.436028</td><td> 7.108271</td><td> 7.698663</td><td> 7.776424</td><td> 7.834396</td><td> 8.494741</td><td> 7.122767</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-5_at </td><td> 7.468834</td><td> 7.259739</td><td> 7.414201</td><td> 7.146814</td><td> 7.486227</td><td> 7.626687</td><td> 7.869447</td><td> 7.963197</td><td> 7.295890</td><td>⋯</td><td> 7.759427</td><td> 7.746014</td><td> 8.219686</td><td> 6.786589</td><td> 7.495551</td><td> 7.639779</td><td> 7.601022</td><td> 8.182174</td><td> 7.408779</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-M_at </td><td> 7.434691</td><td> 7.321820</td><td> 7.682377</td><td> 7.154048</td><td> 7.785387</td><td> 7.701977</td><td> 8.108297</td><td> 8.219608</td><td> 7.643441</td><td>⋯</td><td> 7.934333</td><td> 7.800455</td><td> 8.244886</td><td> 7.095985</td><td> 7.713167</td><td> 7.854818</td><td> 7.814906</td><td> 8.254204</td><td> 7.467042</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioC-3_at </td><td> 9.479190</td><td> 9.172470</td><td> 9.506034</td><td> 9.104503</td><td> 9.419863</td><td> 9.565736</td><td> 9.918638</td><td> 9.945004</td><td> 9.470716</td><td>⋯</td><td> 9.606407</td><td> 9.579577</td><td>10.133814</td><td> 8.948487</td><td> 9.484232</td><td> 9.616063</td><td> 9.542362</td><td>10.107032</td><td> 9.296233</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioC-5_at </td><td> 9.274527</td><td> 9.146463</td><td> 9.368232</td><td> 8.807047</td><td> 9.361186</td><td> 9.346092</td><td> 9.680035</td><td> 9.707252</td><td> 9.291815</td><td>⋯</td><td> 9.507191</td><td> 9.370816</td><td>10.112901</td><td> 8.753260</td><td> 9.231800</td><td> 9.299826</td><td> 9.492625</td><td> 9.787800</td><td> 9.048998</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioD-3_at </td><td>11.962062</td><td>11.685237</td><td>11.919698</td><td>11.449102</td><td>11.878855</td><td>11.790602</td><td>11.950189</td><td>12.111121</td><td>11.973858</td><td>⋯</td><td>11.958430</td><td>11.821378</td><td>12.434823</td><td>11.285018</td><td>11.948128</td><td>11.795919</td><td>11.794224</td><td>12.240711</td><td>11.384425</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioD-5_at </td><td>11.443130</td><td>11.100516</td><td>11.340132</td><td>10.806122</td><td>11.227223</td><td>11.246998</td><td>11.573132</td><td>11.678712</td><td>11.437479</td><td>⋯</td><td>11.420393</td><td>11.357240</td><td>11.865730</td><td>10.760719</td><td>11.295197</td><td>11.349092</td><td>11.286038</td><td>11.754262</td><td>10.951819</td><td></td></tr>
	<tr><td>AFFX-r2-P1-cre-3_at  </td><td>13.069501</td><td>12.991296</td><td>13.239194</td><td>12.787590</td><td>13.057594</td><td>13.055934</td><td>13.090817</td><td>13.302660</td><td>13.198605</td><td>⋯</td><td>13.247024</td><td>13.080788</td><td>13.383656</td><td>12.781179</td><td>13.273540</td><td>13.047362</td><td>13.028270</td><td>13.323030</td><td>12.919859</td><td></td></tr>
	<tr><td>AFFX-r2-P1-cre-5_at  </td><td>13.148925</td><td>12.958498</td><td>13.129225</td><td>12.763438</td><td>12.988927</td><td>12.984767</td><td>13.126631</td><td>13.227394</td><td>13.175285</td><td>⋯</td><td>13.257020</td><td>12.954479</td><td>13.459619</td><td>12.726078</td><td>13.266652</td><td>13.025837</td><td>12.892184</td><td>13.352292</td><td>12.774202</td><td></td></tr>
	<tr><td>AFFX-ThrX-3_at       </td><td> 4.023447</td><td> 3.894434</td><td> 3.831372</td><td> 3.904646</td><td> 3.781000</td><td> 3.863492</td><td> 8.273009</td><td> 4.017678</td><td> 3.819183</td><td>⋯</td><td> 7.736063</td><td> 4.635044</td><td> 6.864505</td><td> 3.948825</td><td> 3.844685</td><td> 3.945589</td><td> 7.918078</td><td> 3.863258</td><td> 4.011895</td><td></td></tr>
	<tr><td>AFFX-ThrX-5_at       </td><td> 3.441101</td><td> 3.599296</td><td> 3.424761</td><td> 3.455952</td><td> 3.398223</td><td> 3.493427</td><td> 6.734710</td><td> 3.332242</td><td> 3.545885</td><td>⋯</td><td> 6.512621</td><td> 3.998604</td><td> 5.359918</td><td> 3.472231</td><td> 3.609189</td><td> 3.346578</td><td> 6.637332</td><td> 3.579053</td><td> 3.558303</td><td></td></tr>
	<tr><td>AFFX-ThrX-M_at       </td><td> 3.193134</td><td> 3.160905</td><td> 3.281208</td><td> 3.306119</td><td> 3.155458</td><td> 3.193315</td><td> 7.431392</td><td> 3.353211</td><td> 3.296945</td><td>⋯</td><td> 6.828242</td><td> 3.506016</td><td> 5.355181</td><td> 3.202003</td><td> 3.303350</td><td> 3.349078</td><td> 7.219454</td><td> 3.420575</td><td> 3.300860</td><td></td></tr>
	<tr><td>AFFX-TrpnX-3_at      </td><td> 3.118769</td><td> 3.111432</td><td> 3.049101</td><td> 3.144660</td><td> 2.976319</td><td> 3.133968</td><td> 3.226229</td><td> 3.112736</td><td> 3.155507</td><td>⋯</td><td> 3.062192</td><td> 3.335672</td><td> 3.117348</td><td> 3.194541</td><td> 3.013140</td><td> 3.109277</td><td> 3.031018</td><td> 3.178001</td><td> 3.199848</td><td></td></tr>
	<tr><td>AFFX-TrpnX-5_at      </td><td> 3.427693</td><td> 3.500176</td><td> 3.512857</td><td> 3.441617</td><td> 3.462671</td><td> 3.438869</td><td> 3.724412</td><td> 3.518988</td><td> 3.379785</td><td>⋯</td><td> 3.469691</td><td> 3.739438</td><td> 3.683459</td><td> 3.382735</td><td> 3.518664</td><td> 3.389215</td><td> 3.466611</td><td> 3.599959</td><td> 3.528824</td><td></td></tr>
	<tr><td>AFFX-TrpnX-M_at      </td><td> 3.325586</td><td> 3.272411</td><td> 3.287203</td><td> 3.456929</td><td> 3.133818</td><td> 3.299522</td><td> 3.779383</td><td> 3.476347</td><td> 3.241239</td><td>⋯</td><td> 3.364435</td><td> 3.711670</td><td> 3.494359</td><td> 3.385121</td><td> 3.395457</td><td> 3.301241</td><td> 3.388843</td><td> 3.263486</td><td> 3.435593</td><td></td></tr>
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
<ol class=list-inline><li>'GSM714070'</li><li>'GSM714071'</li><li>'GSM714072'</li><li>'GSM714073'</li><li>'GSM714074'</li><li>'GSM714075'</li><li>'GSM714076'</li><li>'GSM714077'</li><li>'GSM714078'</li><li>'GSM714079'</li><li>'GSM714080'</li><li>'GSM714081'</li><li>'GSM714082'</li><li>'GSM714083'</li><li>'GSM714084'</li><li>'GSM714085'</li><li>'GSM714086'</li><li>'GSM714087'</li><li>'GSM714088'</li><li>'GSM714089'</li><li>'GSM714090'</li><li>'GSM714091'</li><li>'GSM714092'</li><li>'GSM714093'</li><li>'GSM714094'</li><li>'GSM714095'</li><li>'GSM714096'</li><li>'GSM714097'</li><li>'GSM714098'</li></ol>




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


'GSM714070=mean(GSM714070),GSM714071=mean(GSM714071),GSM714072=mean(GSM714072),GSM714073=mean(GSM714073),GSM714074=mean(GSM714074),GSM714075=mean(GSM714075),GSM714076=mean(GSM714076),GSM714077=mean(GSM714077),GSM714078=mean(GSM714078),GSM714079=mean(GSM714079),GSM714080=mean(GSM714080),GSM714081=mean(GSM714081),GSM714082=mean(GSM714082),GSM714083=mean(GSM714083),GSM714084=mean(GSM714084),GSM714085=mean(GSM714085),GSM714086=mean(GSM714086),GSM714087=mean(GSM714087),GSM714088=mean(GSM714088),GSM714089=mean(GSM714089),GSM714090=mean(GSM714090),GSM714091=mean(GSM714091),GSM714092=mean(GSM714092),GSM714093=mean(GSM714093),GSM714094=mean(GSM714094),GSM714095=mean(GSM714095),GSM714096=mean(GSM714096),GSM714097=mean(GSM714097),GSM714098=mean(GSM714098),'



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
<caption>A data.frame: 6 × 30</caption>
<thead>
	<tr><th></th><th scope=col>GSM714070</th><th scope=col>GSM714071</th><th scope=col>GSM714072</th><th scope=col>GSM714073</th><th scope=col>GSM714074</th><th scope=col>GSM714075</th><th scope=col>GSM714076</th><th scope=col>GSM714077</th><th scope=col>GSM714078</th><th scope=col>GSM714079</th><th scope=col>⋯</th><th scope=col>GSM714090</th><th scope=col>GSM714091</th><th scope=col>GSM714092</th><th scope=col>GSM714093</th><th scope=col>GSM714094</th><th scope=col>GSM714095</th><th scope=col>GSM714096</th><th scope=col>GSM714097</th><th scope=col>GSM714098</th><th scope=col>symbol</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>7.476734</td><td>7.560470</td><td>7.545101</td><td>7.137021</td><td>7.313099</td><td>7.320076</td><td>6.757667</td><td>7.379789</td><td>7.545558</td><td>7.491856</td><td>⋯</td><td>7.427269</td><td>7.046685</td><td>7.362381</td><td>7.537203</td><td>7.467394</td><td>7.281994</td><td>7.340640</td><td>7.311836</td><td>7.066598</td><td>DDR1 /// MIR4640</td></tr>
	<tr><th scope=row>2</th><td>5.231710</td><td>5.571284</td><td>5.502338</td><td>5.406858</td><td>5.961533</td><td>5.719497</td><td>5.464369</td><td>5.687824</td><td>5.465359</td><td>5.275807</td><td>⋯</td><td>5.662993</td><td>5.803273</td><td>5.406180</td><td>5.713884</td><td>5.389508</td><td>5.171387</td><td>5.768445</td><td>5.735633</td><td>5.617395</td><td>RFC2            </td></tr>
	<tr><th scope=row>3</th><td>6.407551</td><td>6.456310</td><td>5.912837</td><td>6.138684</td><td>6.431826</td><td>6.274274</td><td>5.657739</td><td>5.846036</td><td>5.710486</td><td>6.053606</td><td>⋯</td><td>5.860925</td><td>6.435217</td><td>5.536453</td><td>6.086993</td><td>5.458502</td><td>5.751894</td><td>5.863148</td><td>4.882628</td><td>5.738776</td><td>HSPA6           </td></tr>
	<tr><th scope=row>4</th><td>8.086023</td><td>8.056792</td><td>7.715602</td><td>7.423653</td><td>7.747922</td><td>7.833796</td><td>7.633695</td><td>7.562925</td><td>7.985645</td><td>7.753887</td><td>⋯</td><td>7.660248</td><td>7.200921</td><td>8.062172</td><td>7.596732</td><td>7.435648</td><td>7.382841</td><td>7.493919</td><td>7.492598</td><td>7.106533</td><td>PAX8            </td></tr>
	<tr><th scope=row>5</th><td>3.334447</td><td>3.174585</td><td>3.156232</td><td>3.277893</td><td>3.197367</td><td>3.195213</td><td>3.268658</td><td>3.442809</td><td>3.269492</td><td>3.127979</td><td>⋯</td><td>3.269480</td><td>3.368173</td><td>3.337614</td><td>3.249723</td><td>3.268982</td><td>3.363065</td><td>3.181215</td><td>3.228107</td><td>3.360925</td><td>GUCA1A          </td></tr>
	<tr><th scope=row>6</th><td>7.727043</td><td>7.573688</td><td>7.457024</td><td>7.626756</td><td>7.616204</td><td>7.706446</td><td>7.274741</td><td>7.540043</td><td>7.495845</td><td>7.537099</td><td>⋯</td><td>7.439052</td><td>7.486536</td><td>8.113033</td><td>7.569026</td><td>7.275367</td><td>7.250188</td><td>7.376961</td><td>7.086117</td><td>7.022475</td><td>MIR5193 /// UBA7</td></tr>
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
              summarise(GSM714070=mean(GSM714070),GSM714071=mean(GSM714071),GSM714072=mean(GSM714072),GSM714073=mean(GSM714073),GSM714074=mean(GSM714074),GSM714075=mean(GSM714075),GSM714076=mean(GSM714076),GSM714077=mean(GSM714077),GSM714078=mean(GSM714078),GSM714079=mean(GSM714079),GSM714080=mean(GSM714080),GSM714081=mean(GSM714081),GSM714082=mean(GSM714082),GSM714083=mean(GSM714083),GSM714084=mean(GSM714084),GSM714085=mean(GSM714085),GSM714086=mean(GSM714086),GSM714087=mean(GSM714087),GSM714088=mean(GSM714088),GSM714089=mean(GSM714089),GSM714090=mean(GSM714090),GSM714091=mean(GSM714091),GSM714092=mean(GSM714092),GSM714093=mean(GSM714093),GSM714094=mean(GSM714094),GSM714095=mean(GSM714095),GSM714096=mean(GSM714096),GSM714097=mean(GSM714097),GSM714098=mean(GSM714098))  %>% 
         column_to_rownames(var ='symbol')
```


```R
gset$GSE28829_series_matrix.txt.gz@phenoData@data  %>%  head()
```


<table class="dataframe">
<caption>A data.frame: 6 × 35</caption>
<thead>
	<tr><th></th><th scope=col>title</th><th scope=col>geo_accession</th><th scope=col>status</th><th scope=col>submission_date</th><th scope=col>last_update_date</th><th scope=col>type</th><th scope=col>channel_count</th><th scope=col>source_name_ch1</th><th scope=col>organism_ch1</th><th scope=col>characteristics_ch1</th><th scope=col>⋯</th><th scope=col>contact_institute</th><th scope=col>contact_address</th><th scope=col>contact_city</th><th scope=col>contact_state</th><th scope=col>contact_zip/postal_code</th><th scope=col>contact_country</th><th scope=col>supplementary_file</th><th scope=col>data_row_count</th><th scope=col>phenotype:ch1</th><th scope=col>tissue:ch1</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM714070</th><td>advanced atherosclerotic plaque1</td><td>GSM714070</td><td>Public on Apr 25 2011</td><td>Apr 25 2011</td><td>Apr 25 2011</td><td>RNA</td><td>1</td><td>postmortem plaque</td><td>Homo sapiens</td><td>tissue: carotid artery</td><td>⋯</td><td>University of Maastricht</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>Limburg</td><td>6229  HX</td><td>Netherlands</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM714nnn/GSM714070/suppl/GSM714070_ADV1.CEL.gz</td><td>54675</td><td>advanced atherosclerotic plaque</td><td>carotid artery</td></tr>
	<tr><th scope=row>GSM714071</th><td>advanced atherosclerotic plaque2</td><td>GSM714071</td><td>Public on Apr 25 2011</td><td>Apr 25 2011</td><td>Apr 25 2011</td><td>RNA</td><td>1</td><td>postmortem plaque</td><td>Homo sapiens</td><td>tissue: carotid artery</td><td>⋯</td><td>University of Maastricht</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>Limburg</td><td>6229  HX</td><td>Netherlands</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM714nnn/GSM714071/suppl/GSM714071_ADV2.CEL.gz</td><td>54675</td><td>advanced atherosclerotic plaque</td><td>carotid artery</td></tr>
	<tr><th scope=row>GSM714072</th><td>advanced atherosclerotic plaque3</td><td>GSM714072</td><td>Public on Apr 25 2011</td><td>Apr 25 2011</td><td>Apr 25 2011</td><td>RNA</td><td>1</td><td>postmortem plaque</td><td>Homo sapiens</td><td>tissue: carotid artery</td><td>⋯</td><td>University of Maastricht</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>Limburg</td><td>6229  HX</td><td>Netherlands</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM714nnn/GSM714072/suppl/GSM714072_ADV3.CEL.gz</td><td>54675</td><td>advanced atherosclerotic plaque</td><td>carotid artery</td></tr>
	<tr><th scope=row>GSM714073</th><td>advanced atherosclerotic plaque4</td><td>GSM714073</td><td>Public on Apr 25 2011</td><td>Apr 25 2011</td><td>Apr 25 2011</td><td>RNA</td><td>1</td><td>postmortem plaque</td><td>Homo sapiens</td><td>tissue: carotid artery</td><td>⋯</td><td>University of Maastricht</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>Limburg</td><td>6229  HX</td><td>Netherlands</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM714nnn/GSM714073/suppl/GSM714073_ADV4.CEL.gz</td><td>54675</td><td>advanced atherosclerotic plaque</td><td>carotid artery</td></tr>
	<tr><th scope=row>GSM714074</th><td>advanced atherosclerotic plaque5</td><td>GSM714074</td><td>Public on Apr 25 2011</td><td>Apr 25 2011</td><td>Apr 25 2011</td><td>RNA</td><td>1</td><td>postmortem plaque</td><td>Homo sapiens</td><td>tissue: carotid artery</td><td>⋯</td><td>University of Maastricht</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>Limburg</td><td>6229  HX</td><td>Netherlands</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM714nnn/GSM714074/suppl/GSM714074_ADV5.CEL.gz</td><td>54675</td><td>advanced atherosclerotic plaque</td><td>carotid artery</td></tr>
	<tr><th scope=row>GSM714075</th><td>advanced atherosclerotic plaque6</td><td>GSM714075</td><td>Public on Apr 25 2011</td><td>Apr 25 2011</td><td>Apr 25 2011</td><td>RNA</td><td>1</td><td>postmortem plaque</td><td>Homo sapiens</td><td>tissue: carotid artery</td><td>⋯</td><td>University of Maastricht</td><td>P. Debyelaan 25</td><td>Maastricht</td><td>Limburg</td><td>6229  HX</td><td>Netherlands</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM714nnn/GSM714075/suppl/GSM714075_ADV6.CEL.gz</td><td>54675</td><td>advanced atherosclerotic plaque</td><td>carotid artery</td></tr>
</tbody>
</table>




```R
group_list  <- factor(gset$GSE28829_series_matrix.txt.gz@phenoData@data[,c('phenotype:ch1')]) 
```


```R
levels(group_list)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'advanced atherosclerotic plaque'</li><li>'early atherosclerotic plaque'</li></ol>




```R
group_list  <- fct_recode(group_list, 'EAR' = 'early atherosclerotic plaque',
                                      'ADV' = 'advanced atherosclerotic plaque')
```


```R
group_list <- relevel(group_list, ref="EAR")
```


```R
group_list
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>ADV</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li><li>EAR</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'EAR'</li><li>'ADV'</li></ol>
</details>



```R
library(ggpubr)
```


```R
selected_gene <- c('ALDOA','CREG1','LGMN','PKM')
```


```R
normalized_gset_mean[selected_gene,]
```


<table class="dataframe">
<caption>A data.frame: 4 × 29</caption>
<thead>
	<tr><th></th><th scope=col>GSM714070</th><th scope=col>GSM714071</th><th scope=col>GSM714072</th><th scope=col>GSM714073</th><th scope=col>GSM714074</th><th scope=col>GSM714075</th><th scope=col>GSM714076</th><th scope=col>GSM714077</th><th scope=col>GSM714078</th><th scope=col>GSM714079</th><th scope=col>⋯</th><th scope=col>GSM714089</th><th scope=col>GSM714090</th><th scope=col>GSM714091</th><th scope=col>GSM714092</th><th scope=col>GSM714093</th><th scope=col>GSM714094</th><th scope=col>GSM714095</th><th scope=col>GSM714096</th><th scope=col>GSM714097</th><th scope=col>GSM714098</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ALDOA</th><td> 9.893794</td><td> 9.573743</td><td> 9.716425</td><td> 9.745118</td><td> 9.830804</td><td> 9.543543</td><td> 9.575849</td><td>10.194999</td><td> 9.833627</td><td> 9.903148</td><td>⋯</td><td>10.004065</td><td> 9.794863</td><td>10.153969</td><td>10.349509</td><td> 9.629413</td><td> 8.987011</td><td> 9.826197</td><td> 9.437448</td><td>10.100124</td><td>10.161504</td></tr>
	<tr><th scope=row>CREG1</th><td>11.069110</td><td>11.330126</td><td>10.870153</td><td>10.652426</td><td>11.603152</td><td>11.340677</td><td>10.328475</td><td>10.857587</td><td>11.025890</td><td>10.565698</td><td>⋯</td><td>10.794372</td><td>10.953031</td><td>10.497119</td><td>11.138544</td><td>11.048825</td><td>10.997920</td><td>10.763623</td><td>11.118992</td><td>10.554302</td><td>10.433053</td></tr>
	<tr><th scope=row>LGMN</th><td> 9.871655</td><td>10.129602</td><td> 9.705219</td><td> 9.564114</td><td>10.355589</td><td>10.302197</td><td> 8.181282</td><td> 9.780105</td><td> 9.754146</td><td> 9.459813</td><td>⋯</td><td> 8.796401</td><td> 9.008777</td><td> 7.786350</td><td> 8.038617</td><td> 9.342329</td><td> 8.713496</td><td> 8.496280</td><td> 9.458924</td><td> 8.631696</td><td> 8.304162</td></tr>
	<tr><th scope=row>PKM</th><td> 8.043051</td><td> 7.730174</td><td> 8.176608</td><td> 8.212995</td><td> 7.797923</td><td> 7.507078</td><td> 7.932536</td><td> 8.255880</td><td> 7.875856</td><td> 8.523129</td><td>⋯</td><td> 7.836404</td><td> 7.611898</td><td> 8.310836</td><td> 8.111405</td><td> 7.676452</td><td> 7.753795</td><td> 8.227636</td><td> 7.471816</td><td> 7.791793</td><td> 8.039212</td></tr>
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
<caption>A data.frame: 29 × 6</caption>
<thead>
	<tr><th scope=col>sample</th><th scope=col>ALDOA</th><th scope=col>CREG1</th><th scope=col>LGMN</th><th scope=col>PKM</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM714070</td><td> 9.893794</td><td>11.06911</td><td> 9.871655</td><td>8.043051</td><td>ADV</td></tr>
	<tr><td>GSM714071</td><td> 9.573743</td><td>11.33013</td><td>10.129602</td><td>7.730174</td><td>ADV</td></tr>
	<tr><td>GSM714072</td><td> 9.716425</td><td>10.87015</td><td> 9.705219</td><td>8.176608</td><td>ADV</td></tr>
	<tr><td>GSM714073</td><td> 9.745118</td><td>10.65243</td><td> 9.564114</td><td>8.212995</td><td>ADV</td></tr>
	<tr><td>GSM714074</td><td> 9.830804</td><td>11.60315</td><td>10.355589</td><td>7.797923</td><td>ADV</td></tr>
	<tr><td>GSM714075</td><td> 9.543543</td><td>11.34068</td><td>10.302197</td><td>7.507078</td><td>ADV</td></tr>
	<tr><td>GSM714076</td><td> 9.575849</td><td>10.32847</td><td> 8.181282</td><td>7.932536</td><td>ADV</td></tr>
	<tr><td>GSM714077</td><td>10.194999</td><td>10.85759</td><td> 9.780105</td><td>8.255880</td><td>ADV</td></tr>
	<tr><td>GSM714078</td><td> 9.833627</td><td>11.02589</td><td> 9.754146</td><td>7.875856</td><td>ADV</td></tr>
	<tr><td>GSM714079</td><td> 9.903148</td><td>10.56570</td><td> 9.459813</td><td>8.523129</td><td>ADV</td></tr>
	<tr><td>GSM714080</td><td> 9.838939</td><td>10.56080</td><td> 9.492055</td><td>8.128461</td><td>ADV</td></tr>
	<tr><td>GSM714081</td><td>10.088661</td><td>11.03958</td><td> 9.391591</td><td>7.967070</td><td>ADV</td></tr>
	<tr><td>GSM714082</td><td> 9.873660</td><td>10.91040</td><td> 9.158378</td><td>8.001972</td><td>ADV</td></tr>
	<tr><td>GSM714083</td><td> 8.937419</td><td>10.66023</td><td> 8.272124</td><td>8.078396</td><td>ADV</td></tr>
	<tr><td>GSM714084</td><td> 9.861846</td><td>11.16510</td><td>10.231748</td><td>8.378108</td><td>ADV</td></tr>
	<tr><td>GSM714085</td><td> 9.571141</td><td>11.22439</td><td>10.194197</td><td>8.152337</td><td>ADV</td></tr>
	<tr><td>GSM714086</td><td> 9.953937</td><td>10.82712</td><td> 8.366112</td><td>8.363685</td><td>EAR</td></tr>
	<tr><td>GSM714087</td><td> 9.740244</td><td>10.79503</td><td> 8.467059</td><td>7.887106</td><td>EAR</td></tr>
	<tr><td>GSM714088</td><td> 9.595981</td><td>10.96080</td><td> 8.606146</td><td>7.519633</td><td>EAR</td></tr>
	<tr><td>GSM714089</td><td>10.004065</td><td>10.79437</td><td> 8.796401</td><td>7.836404</td><td>EAR</td></tr>
	<tr><td>GSM714090</td><td> 9.794863</td><td>10.95303</td><td> 9.008777</td><td>7.611898</td><td>EAR</td></tr>
	<tr><td>GSM714091</td><td>10.153969</td><td>10.49712</td><td> 7.786350</td><td>8.310836</td><td>EAR</td></tr>
	<tr><td>GSM714092</td><td>10.349509</td><td>11.13854</td><td> 8.038617</td><td>8.111405</td><td>EAR</td></tr>
	<tr><td>GSM714093</td><td> 9.629413</td><td>11.04883</td><td> 9.342329</td><td>7.676452</td><td>EAR</td></tr>
	<tr><td>GSM714094</td><td> 8.987011</td><td>10.99792</td><td> 8.713496</td><td>7.753795</td><td>EAR</td></tr>
	<tr><td>GSM714095</td><td> 9.826197</td><td>10.76362</td><td> 8.496280</td><td>8.227636</td><td>EAR</td></tr>
	<tr><td>GSM714096</td><td> 9.437448</td><td>11.11899</td><td> 9.458924</td><td>7.471816</td><td>EAR</td></tr>
	<tr><td>GSM714097</td><td>10.100124</td><td>10.55430</td><td> 8.631696</td><td>7.791793</td><td>EAR</td></tr>
	<tr><td>GSM714098</td><td>10.161504</td><td>10.43305</td><td> 8.304162</td><td>8.039212</td><td>EAR</td></tr>
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
         stat_compare_means(comparisons = list(c('EAR','ADV')),method = 't.test', paired = F,label = 'p.signif') 
```


    
![png](Step11.1_GSE28829_files/Step11.1_GSE28829_24_0.png)
    



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
         stat_compare_means(comparisons = list(c('EAR','ADV')),method = 't.test', paired = F,label = 'p.signif')  + 
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
         stat_compare_means(comparisons = list(c('EAR','ADV')),method = 't.test', paired = F,label = 'p.signif')  + 
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
         stat_compare_means(comparisons = list(c('EAR','ADV')),method = 't.test', paired = F,label = 'p.signif')  + 
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
         stat_compare_means(comparisons = list(c('EAR','ADV')),method = 't.test', paired = F,label = 'p.signif')  + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), axis.title.y = element_text(size=0), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p4
```


    
![png](Step11.1_GSE28829_files/Step11.1_GSE28829_26_0.png)
    



    
![png](Step11.1_GSE28829_files/Step11.1_GSE28829_26_1.png)
    



    
![png](Step11.1_GSE28829_files/Step11.1_GSE28829_26_2.png)
    



    
![png](Step11.1_GSE28829_files/Step11.1_GSE28829_26_3.png)
    



```R
p0 <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 4,guides='collect')  +  
    plot_annotation(title = 'GSE28829',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE28829_4Gene.pdf', height = 5, width = 18, units = 'cm')
```


    
![png](Step11.1_GSE28829_files/Step11.1_GSE28829_27_0.png)
    

