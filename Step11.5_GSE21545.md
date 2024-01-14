```R
########
library(GEOquery) 
library(tidyverse)
library(AnnoProbe)
#eSet <- geoChina(gse_number, destdir = '.')
```

    Loading required package: Biobase
    
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
    
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Setting options('download.file.method.GEOquery'='auto')
    
    Setting options('GEOquery.inmemory.gpl'=FALSE)
    
    Warning message in system("timedatectl", intern = TRUE):
    “running command 'timedatectl' had status 1”
    ── [1mAttaching packages[22m ─────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
    
    [32m✔[39m [34mggplot2[39m 3.4.0     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.5     [32m✔[39m [34mdplyr  [39m 1.0.7
    [32m✔[39m [34mtidyr  [39m 1.1.4     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 2.0.2     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mcombine()[39m    masks [34mBiobase[39m::combine(), [34mBiocGenerics[39m::combine()
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m     masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m        masks [34mstats[39m::lag()
    [31m✖[39m [34mggplot2[39m::[32mPosition()[39m masks [34mBiocGenerics[39m::Position(), [34mbase[39m::Position()
    
    AnnoProbe v 0.1.6  welcome to use AnnoProbe!
    If you use AnnoProbe in published research, please acknowledgements:
    We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee, for generously sharing their experience and codes.
    



```R
Sys.setenv("VROOM_CONNECTION_SIZE"= 500000 * 2)
```


```R
GSE21545  <-  GEOquery::getGEO('GSE21545', GSEMatrix = TRUE,  getGPL = F)
```

    Found 1 file(s)
    
    GSE21545_series_matrix.txt.gz
    
    [1mRows: [22m[34m54675[39m [1mColumns: [22m[34m224[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m "\t"
    [31mchr[39m   (1): ID_REF
    [32mdbl[39m (223): GSM892518, GSM892519, GSM892520, GSM892521, GSM892522, GSM892523,...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
GSE21545[[1]]
```


    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 54675 features, 223 samples 
      element names: exprs 
    protocolData: none
    phenoData
      sampleNames: GSM892518 GSM892519 ... GSM892740 (223 total)
      varLabels: title geo_accession ... sample type:ch1 (43 total)
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
      pubMedIds: 22371308
    27437576
    30444878
    30175270
    31312755
    33959644
    35494231 
    Annotation: GPL570 



```R
exprs_matrix  <- exprs(GSE21545[[1]])
```


```R
ids <- AnnoProbe::idmap('GPL570',type = 'soft')
ids   <- ids  %>% 
           mutate(ID = as.character(ID))
exprs_matrix %>% 
 as.data.frame() %>% 
  rownames_to_column(var = 'Probe') %>%  
   left_join(ids, by = c('Probe'= 'ID'))
```

    file downloaded in /home/jfckkiu/AS_HG/Final_Results/Figure7_RNA_Validation_ALDOA
    



<table class="dataframe">
<caption>A data.frame: 54675 × 225</caption>
<thead>
	<tr><th scope=col>Probe</th><th scope=col>GSM892518</th><th scope=col>GSM892519</th><th scope=col>GSM892520</th><th scope=col>GSM892521</th><th scope=col>GSM892522</th><th scope=col>GSM892523</th><th scope=col>GSM892524</th><th scope=col>GSM892525</th><th scope=col>GSM892526</th><th scope=col>⋯</th><th scope=col>GSM892732</th><th scope=col>GSM892733</th><th scope=col>GSM892734</th><th scope=col>GSM892735</th><th scope=col>GSM892736</th><th scope=col>GSM892737</th><th scope=col>GSM892738</th><th scope=col>GSM892739</th><th scope=col>GSM892740</th><th scope=col>symbol</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1007_s_at   </td><td> 6.679589</td><td> 6.550783</td><td> 7.049459</td><td> 6.613194</td><td> 6.230651</td><td> 6.865954</td><td> 6.777266</td><td> 6.994077</td><td> 6.666644</td><td>⋯</td><td>7.496145</td><td>7.163326</td><td>6.965941</td><td>7.767897</td><td>7.911350</td><td>6.759342</td><td>6.706558</td><td>7.403885</td><td>7.597200</td><td>DDR1 /// MIR4640          </td></tr>
	<tr><td>1053_at     </td><td> 7.859804</td><td> 7.199514</td><td> 7.525922</td><td> 7.619245</td><td> 7.768536</td><td> 7.516350</td><td> 7.168577</td><td> 7.487657</td><td> 7.643404</td><td>⋯</td><td>6.531912</td><td>6.437202</td><td>5.802015</td><td>6.149521</td><td>5.952339</td><td>6.115841</td><td>6.463664</td><td>6.739720</td><td>6.043171</td><td>RFC2                      </td></tr>
	<tr><td>117_at      </td><td> 8.367065</td><td> 8.425298</td><td> 8.680992</td><td> 8.384369</td><td>10.724641</td><td> 7.793695</td><td> 6.895588</td><td> 8.097322</td><td> 8.512935</td><td>⋯</td><td>8.227065</td><td>8.345329</td><td>8.003220</td><td>7.326053</td><td>7.464806</td><td>7.196055</td><td>7.680891</td><td>9.010450</td><td>7.893110</td><td>HSPA6                     </td></tr>
	<tr><td>121_at      </td><td> 7.302292</td><td> 7.097573</td><td> 6.861540</td><td> 7.152322</td><td> 7.212102</td><td> 7.019013</td><td> 7.268188</td><td> 7.153352</td><td> 7.170580</td><td>⋯</td><td>7.890883</td><td>6.993632</td><td>7.720437</td><td>7.380170</td><td>7.395961</td><td>7.360369</td><td>6.802845</td><td>7.682543</td><td>7.461491</td><td>PAX8                      </td></tr>
	<tr><td>1255_g_at   </td><td> 2.622218</td><td> 2.357024</td><td> 2.480694</td><td> 2.586321</td><td> 2.271802</td><td> 2.458667</td><td> 2.574721</td><td> 2.479205</td><td> 2.420047</td><td>⋯</td><td>2.777725</td><td>2.589518</td><td>2.806212</td><td>2.657732</td><td>2.601428</td><td>2.812067</td><td>2.617467</td><td>2.952880</td><td>2.644002</td><td>GUCA1A                    </td></tr>
	<tr><td>1294_at     </td><td> 8.382746</td><td> 8.279175</td><td> 8.614102</td><td> 8.548336</td><td> 8.656825</td><td> 8.559876</td><td> 8.961890</td><td> 8.579970</td><td> 8.852497</td><td>⋯</td><td>7.348823</td><td>7.675665</td><td>7.318562</td><td>7.590046</td><td>7.767733</td><td>7.233156</td><td>7.656409</td><td>7.424672</td><td>7.761288</td><td>MIR5193 /// UBA7          </td></tr>
	<tr><td>1316_at     </td><td> 5.538944</td><td> 5.774041</td><td> 5.520036</td><td> 5.790884</td><td> 5.692584</td><td> 5.802673</td><td> 5.840265</td><td> 5.675343</td><td> 5.560076</td><td>⋯</td><td>5.413612</td><td>5.090811</td><td>4.872417</td><td>5.293664</td><td>5.304146</td><td>4.794719</td><td>4.781214</td><td>5.249890</td><td>5.162999</td><td>THRA                      </td></tr>
	<tr><td>1320_at     </td><td> 3.735642</td><td> 3.704137</td><td> 3.702092</td><td> 3.918927</td><td> 3.780438</td><td> 3.814116</td><td> 3.569832</td><td> 3.882318</td><td> 3.792257</td><td>⋯</td><td>4.569263</td><td>4.582414</td><td>4.353643</td><td>4.670260</td><td>4.822963</td><td>5.098755</td><td>4.384210</td><td>4.171604</td><td>4.667133</td><td>PTPN21                    </td></tr>
	<tr><td>1405_i_at   </td><td>12.427214</td><td>12.379755</td><td>12.479644</td><td>11.887506</td><td>11.125330</td><td>12.460502</td><td>13.207017</td><td>12.421032</td><td>11.840692</td><td>⋯</td><td>9.236070</td><td>9.784224</td><td>9.579448</td><td>8.933905</td><td>9.134253</td><td>8.270661</td><td>8.778715</td><td>8.915736</td><td>9.499077</td><td>CCL5                      </td></tr>
	<tr><td>1431_at     </td><td> 3.493767</td><td> 3.355712</td><td> 3.245078</td><td> 3.680550</td><td> 3.279678</td><td> 3.593807</td><td> 3.726675</td><td> 3.663594</td><td> 3.303323</td><td>⋯</td><td>3.221056</td><td>3.224733</td><td>3.245364</td><td>3.282066</td><td>3.231831</td><td>2.950733</td><td>3.343796</td><td>3.443255</td><td>3.262407</td><td>CYP2E1                    </td></tr>
	<tr><td>1438_at     </td><td> 4.707948</td><td> 4.624310</td><td> 4.743726</td><td> 4.718799</td><td> 4.740787</td><td> 4.817917</td><td> 4.506859</td><td> 4.909037</td><td> 4.679351</td><td>⋯</td><td>4.731242</td><td>4.906842</td><td>4.681205</td><td>4.925192</td><td>5.019283</td><td>4.451208</td><td>4.357638</td><td>4.715933</td><td>4.811465</td><td>EPHB3                     </td></tr>
	<tr><td>1487_at     </td><td> 7.600412</td><td> 7.698399</td><td> 7.545616</td><td> 7.795813</td><td> 7.895067</td><td> 7.518234</td><td> 7.486235</td><td> 7.698046</td><td> 7.497911</td><td>⋯</td><td>7.124869</td><td>6.977318</td><td>6.329098</td><td>6.582808</td><td>6.684281</td><td>6.175827</td><td>5.588775</td><td>6.886577</td><td>6.593872</td><td>ESRRA                     </td></tr>
	<tr><td>1494_f_at   </td><td> 4.850182</td><td> 4.692186</td><td> 4.930542</td><td> 5.169379</td><td> 5.197225</td><td> 4.632590</td><td> 5.177664</td><td> 5.024212</td><td> 5.035886</td><td>⋯</td><td>4.844632</td><td>5.038067</td><td>5.003540</td><td>5.004406</td><td>5.034085</td><td>4.693081</td><td>4.731836</td><td>4.925866</td><td>4.856255</td><td>CYP2A6                    </td></tr>
	<tr><td>1552256_a_at</td><td> 5.790458</td><td> 6.411406</td><td> 6.028994</td><td> 5.528090</td><td> 5.542696</td><td> 5.993986</td><td> 5.427683</td><td> 5.789331</td><td> 6.208921</td><td>⋯</td><td>8.560696</td><td>9.035208</td><td>8.244412</td><td>7.637638</td><td>7.248802</td><td>8.481752</td><td>7.589752</td><td>8.268369</td><td>7.907737</td><td>SCARB1                    </td></tr>
	<tr><td>1552257_a_at</td><td> 7.323254</td><td> 7.406329</td><td> 7.427940</td><td> 7.394382</td><td> 7.428375</td><td> 7.423873</td><td> 6.997628</td><td> 7.312785</td><td> 7.188705</td><td>⋯</td><td>7.474827</td><td>7.202959</td><td>6.810416</td><td>6.779700</td><td>7.164774</td><td>6.705766</td><td>6.603712</td><td>6.594670</td><td>6.979144</td><td>TTLL12                    </td></tr>
	<tr><td>1552258_at  </td><td> 4.204247</td><td> 4.387105</td><td> 4.063883</td><td> 4.215286</td><td> 4.645331</td><td> 4.166324</td><td> 4.534478</td><td> 4.014600</td><td> 4.244470</td><td>⋯</td><td>4.214830</td><td>4.273041</td><td>4.208968</td><td>4.125048</td><td>4.030439</td><td>3.687653</td><td>3.700521</td><td>4.689686</td><td>4.047067</td><td>LINC00152 /// LOC101930489</td></tr>
	<tr><td>1552261_at  </td><td> 3.760853</td><td> 3.457963</td><td> 3.503154</td><td> 3.504151</td><td> 3.564234</td><td> 3.481805</td><td> 4.161457</td><td> 3.725731</td><td> 3.586947</td><td>⋯</td><td>3.893318</td><td>3.830473</td><td>4.253811</td><td>4.063046</td><td>4.163851</td><td>4.638884</td><td>3.847763</td><td>3.681185</td><td>3.800669</td><td>WFDC2                     </td></tr>
	<tr><td>1552263_at  </td><td> 8.611584</td><td> 8.635892</td><td> 8.715976</td><td> 8.573162</td><td> 9.569571</td><td> 8.398983</td><td> 8.270812</td><td> 8.423107</td><td> 8.780900</td><td>⋯</td><td>6.944120</td><td>6.681771</td><td>6.701175</td><td>5.962845</td><td>5.976677</td><td>6.238839</td><td>5.910938</td><td>6.873009</td><td>6.146677</td><td>MAPK1                     </td></tr>
	<tr><td>1552264_a_at</td><td> 9.141372</td><td> 8.715761</td><td> 9.089015</td><td> 8.897100</td><td> 9.933612</td><td> 8.882382</td><td> 8.285977</td><td> 9.046600</td><td> 9.183162</td><td>⋯</td><td>7.298019</td><td>6.731071</td><td>7.195554</td><td>6.976353</td><td>6.821065</td><td>7.626919</td><td>7.277828</td><td>7.069859</td><td>7.268221</td><td>MAPK1                     </td></tr>
	<tr><td>1552266_at  </td><td> 2.358678</td><td> 2.229809</td><td> 2.199685</td><td> 2.322233</td><td> 2.603621</td><td> 2.210992</td><td> 2.465447</td><td> 2.243770</td><td> 2.477069</td><td>⋯</td><td>2.922386</td><td>3.376340</td><td>3.648905</td><td>3.007819</td><td>3.112032</td><td>3.298305</td><td>2.934632</td><td>3.002818</td><td>3.091407</td><td>ADAM32                    </td></tr>
	<tr><td>1552269_at  </td><td> 2.363987</td><td> 2.251266</td><td> 2.468833</td><td> 2.626568</td><td> 2.281775</td><td> 2.555919</td><td> 2.864503</td><td> 2.522353</td><td> 2.591598</td><td>⋯</td><td>2.702444</td><td>2.761779</td><td>3.019255</td><td>2.855972</td><td>2.795409</td><td>3.095035</td><td>3.080215</td><td>2.607486</td><td>3.033505</td><td>SPATA17                   </td></tr>
	<tr><td>1552271_at  </td><td> 4.640230</td><td> 4.449571</td><td> 4.762516</td><td> 4.568005</td><td> 4.692367</td><td> 4.598604</td><td> 4.824459</td><td> 4.875617</td><td> 4.828250</td><td>⋯</td><td>4.829769</td><td>4.782295</td><td>4.759207</td><td>4.513405</td><td>4.491487</td><td>4.767900</td><td>4.213869</td><td>4.421339</td><td>4.500796</td><td>PRR22                     </td></tr>
	<tr><td>1552272_a_at</td><td> 4.455310</td><td> 4.100657</td><td> 4.421236</td><td> 4.463934</td><td> 4.239524</td><td> 4.418479</td><td> 4.105108</td><td> 4.546553</td><td> 4.235814</td><td>⋯</td><td>4.669299</td><td>4.560938</td><td>4.382450</td><td>4.555901</td><td>4.547502</td><td>4.132515</td><td>3.816676</td><td>4.348214</td><td>4.467834</td><td>PRR22                     </td></tr>
	<tr><td>1552274_at  </td><td> 6.753764</td><td> 6.824846</td><td> 6.951324</td><td> 6.653160</td><td> 7.957189</td><td> 7.784609</td><td> 6.695751</td><td> 7.300474</td><td> 7.613547</td><td>⋯</td><td>5.780412</td><td>5.434913</td><td>6.171067</td><td>5.956352</td><td>6.297948</td><td>7.056348</td><td>7.758583</td><td>5.671466</td><td>6.452730</td><td>PXK                       </td></tr>
	<tr><td>1552275_s_at</td><td> 7.042399</td><td> 7.188016</td><td> 6.952320</td><td> 6.711789</td><td> 7.982206</td><td> 6.944154</td><td> 5.789157</td><td> 7.018306</td><td> 7.074936</td><td>⋯</td><td>5.938226</td><td>4.990320</td><td>5.824650</td><td>5.584332</td><td>5.688073</td><td>6.408547</td><td>6.227653</td><td>5.959026</td><td>6.106590</td><td>PXK                       </td></tr>
	<tr><td>1552276_a_at</td><td> 5.116415</td><td> 5.275616</td><td> 5.207492</td><td> 5.067102</td><td> 5.351849</td><td> 4.855468</td><td> 5.368724</td><td> 5.175319</td><td> 5.107496</td><td>⋯</td><td>5.043140</td><td>4.412130</td><td>4.484251</td><td>4.693415</td><td>4.403857</td><td>5.052961</td><td>4.254601</td><td>4.473880</td><td>4.671905</td><td>VPS18                     </td></tr>
	<tr><td>1552277_a_at</td><td> 7.070547</td><td> 6.632322</td><td> 7.192575</td><td> 6.424360</td><td> 6.295723</td><td> 7.148551</td><td> 8.685769</td><td> 7.344805</td><td> 7.128713</td><td>⋯</td><td>8.709568</td><td>7.921903</td><td>8.357311</td><td>7.620396</td><td>7.187924</td><td>7.612945</td><td>7.092512</td><td>8.410635</td><td>7.611760</td><td>MSANTD3                   </td></tr>
	<tr><td>1552278_a_at</td><td> 3.731603</td><td> 3.171154</td><td> 3.300317</td><td> 3.128599</td><td> 3.207171</td><td> 3.450796</td><td> 3.043935</td><td> 3.192074</td><td> 3.379181</td><td>⋯</td><td>5.479267</td><td>5.168351</td><td>5.117080</td><td>5.322355</td><td>5.026272</td><td>5.272848</td><td>5.501095</td><td>6.038411</td><td>5.526573</td><td>SLC46A1                   </td></tr>
	<tr><td>1552279_a_at</td><td> 5.128919</td><td> 4.803682</td><td> 5.021896</td><td> 5.139371</td><td> 5.003249</td><td> 4.976917</td><td> 4.658185</td><td> 4.727988</td><td> 5.027161</td><td>⋯</td><td>6.601545</td><td>6.217399</td><td>6.342857</td><td>6.298828</td><td>6.069605</td><td>6.266139</td><td>5.791994</td><td>6.499609</td><td>6.477843</td><td>SLC46A1                   </td></tr>
	<tr><td>1552280_at  </td><td> 3.248613</td><td> 3.165785</td><td> 3.610717</td><td> 3.604242</td><td> 3.407655</td><td> 3.834008</td><td> 3.288295</td><td> 2.985274</td><td> 3.970349</td><td>⋯</td><td>4.874487</td><td>3.139437</td><td>6.392161</td><td>3.229272</td><td>3.430259</td><td>3.729117</td><td>3.924164</td><td>6.850711</td><td>5.377627</td><td>TIMD4                     </td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>AFFX-PheX-3_at       </td><td> 4.598850</td><td> 4.593038</td><td> 4.939596</td><td> 4.716156</td><td> 5.212725</td><td> 4.871588</td><td> 6.587520</td><td> 4.559078</td><td> 4.873095</td><td>⋯</td><td> 7.671242</td><td> 7.095251</td><td> 6.932433</td><td> 6.253713</td><td> 7.402154</td><td> 5.391137</td><td> 7.796067</td><td> 4.755294</td><td> 6.938544</td><td></td></tr>
	<tr><td>AFFX-PheX-5_at       </td><td> 2.948022</td><td> 3.289071</td><td> 3.871543</td><td> 3.170529</td><td> 3.906662</td><td> 3.308499</td><td> 5.714208</td><td> 3.467126</td><td> 3.404479</td><td>⋯</td><td> 6.093308</td><td> 4.888843</td><td> 5.172395</td><td> 4.830196</td><td> 5.902344</td><td> 3.632517</td><td> 7.359934</td><td> 3.451628</td><td> 5.945526</td><td></td></tr>
	<tr><td>AFFX-PheX-M_at       </td><td> 3.512788</td><td> 3.767741</td><td> 4.072215</td><td> 3.749433</td><td> 4.449909</td><td> 4.019325</td><td> 6.119774</td><td> 3.770393</td><td> 4.079216</td><td>⋯</td><td> 7.347003</td><td> 6.615943</td><td> 6.654928</td><td> 5.903176</td><td> 7.046821</td><td> 4.177118</td><td> 8.187283</td><td> 3.989763</td><td> 7.214057</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-3_at  </td><td> 5.500505</td><td> 6.533205</td><td> 7.157503</td><td> 6.514233</td><td> 8.045747</td><td> 6.848733</td><td> 9.981649</td><td> 6.700129</td><td> 7.074279</td><td>⋯</td><td>11.183459</td><td>10.851567</td><td>10.385841</td><td> 9.687287</td><td>10.844909</td><td> 6.359231</td><td>11.025456</td><td> 7.019251</td><td>11.005642</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-5_at  </td><td> 3.649117</td><td> 4.985278</td><td> 5.304532</td><td> 4.024629</td><td> 5.293309</td><td> 5.345589</td><td> 7.907466</td><td> 4.798766</td><td> 4.880565</td><td>⋯</td><td> 8.815082</td><td> 8.051313</td><td> 9.339799</td><td> 7.545550</td><td> 8.704691</td><td> 5.436083</td><td>10.291169</td><td> 4.798968</td><td> 9.918831</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-dap-M_at  </td><td> 5.039503</td><td> 5.893046</td><td> 6.376059</td><td> 5.588172</td><td> 7.148634</td><td> 6.219545</td><td> 9.144369</td><td> 5.967737</td><td> 6.398823</td><td>⋯</td><td> 9.915791</td><td> 9.436178</td><td> 9.587664</td><td> 8.408367</td><td> 9.560547</td><td> 5.563287</td><td>10.358116</td><td> 5.798771</td><td>10.298945</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-3_at  </td><td> 3.069867</td><td> 3.094153</td><td> 3.259531</td><td> 3.070309</td><td> 3.644247</td><td> 3.162961</td><td> 5.535428</td><td> 3.051491</td><td> 3.485043</td><td>⋯</td><td> 7.477420</td><td> 6.957379</td><td> 5.852352</td><td> 5.719437</td><td> 6.866384</td><td> 3.560889</td><td> 7.303582</td><td> 3.587459</td><td> 6.515841</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-5_at  </td><td> 2.730032</td><td> 2.845381</td><td> 3.478408</td><td> 2.852971</td><td> 3.227108</td><td> 3.160855</td><td> 4.930513</td><td> 2.728919</td><td> 2.983155</td><td>⋯</td><td> 6.531420</td><td> 5.461736</td><td> 5.085634</td><td> 4.662033</td><td> 6.054897</td><td> 3.035848</td><td> 7.625180</td><td> 2.861149</td><td> 6.135659</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-lys-M_at  </td><td> 2.744625</td><td> 3.111330</td><td> 3.090429</td><td> 2.637663</td><td> 3.439011</td><td> 2.991128</td><td> 5.224392</td><td> 2.896202</td><td> 3.224233</td><td>⋯</td><td> 7.138676</td><td> 6.209415</td><td> 5.965075</td><td> 5.375557</td><td> 6.767990</td><td> 3.289086</td><td> 7.256461</td><td> 2.767038</td><td> 6.627020</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-3_at  </td><td> 4.758886</td><td> 5.229571</td><td> 5.832445</td><td> 5.309914</td><td> 6.253840</td><td> 5.727051</td><td> 7.955780</td><td> 5.543867</td><td> 5.516492</td><td>⋯</td><td> 9.149036</td><td> 8.585205</td><td> 8.269863</td><td> 7.710047</td><td> 8.710980</td><td> 5.045581</td><td> 8.814096</td><td> 5.201465</td><td> 8.597860</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-5_at  </td><td> 3.173830</td><td> 3.565787</td><td> 4.074979</td><td> 3.299200</td><td> 4.484291</td><td> 3.569928</td><td> 6.418922</td><td> 3.777348</td><td> 3.829655</td><td>⋯</td><td> 6.900966</td><td> 5.980139</td><td> 6.068337</td><td> 5.473513</td><td> 6.709771</td><td> 4.355112</td><td> 8.470416</td><td> 3.762183</td><td> 7.047038</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-phe-M_at  </td><td> 3.030041</td><td> 3.333375</td><td> 4.062231</td><td> 3.167991</td><td> 4.369677</td><td> 3.524606</td><td> 6.412157</td><td> 3.021550</td><td> 3.766191</td><td>⋯</td><td> 7.525586</td><td> 6.782962</td><td> 6.525617</td><td> 5.828970</td><td> 7.082039</td><td> 3.644525</td><td> 7.944444</td><td> 3.480250</td><td> 7.013636</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-3_s_at</td><td> 4.614092</td><td> 5.536753</td><td> 6.055873</td><td> 5.528090</td><td> 6.571493</td><td> 5.638018</td><td> 8.393967</td><td> 5.398113</td><td> 5.764516</td><td>⋯</td><td> 9.458451</td><td> 9.176939</td><td> 8.978787</td><td> 7.755788</td><td> 9.028556</td><td> 5.539396</td><td> 9.556980</td><td> 5.296652</td><td> 9.472929</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-5_s_at</td><td> 3.038524</td><td> 2.991084</td><td> 3.438014</td><td> 3.014909</td><td> 3.975421</td><td> 3.240991</td><td> 5.675220</td><td> 3.288321</td><td> 3.274580</td><td>⋯</td><td> 6.769257</td><td> 5.733184</td><td> 7.644040</td><td> 5.355193</td><td> 6.391538</td><td> 4.413013</td><td> 8.575562</td><td> 3.340780</td><td> 7.938341</td><td></td></tr>
	<tr><td>AFFX-r2-Bs-thr-M_s_at</td><td> 4.587094</td><td> 4.335601</td><td> 4.808331</td><td> 4.862180</td><td> 5.508714</td><td> 4.943562</td><td> 7.422636</td><td> 5.075566</td><td> 5.280243</td><td>⋯</td><td> 8.133500</td><td> 7.480804</td><td> 8.289995</td><td> 6.490385</td><td> 7.567460</td><td> 5.073923</td><td> 9.209276</td><td> 4.537547</td><td> 8.597494</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-3_at </td><td> 7.782976</td><td> 8.220614</td><td> 8.238296</td><td> 7.916081</td><td> 8.738632</td><td> 8.086689</td><td> 9.798906</td><td> 8.207686</td><td> 8.053735</td><td>⋯</td><td> 8.028776</td><td> 7.828738</td><td> 9.212185</td><td> 7.479221</td><td> 7.685919</td><td> 8.020703</td><td> 8.595616</td><td> 8.121945</td><td> 8.242413</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-5_at </td><td> 7.761238</td><td> 8.141143</td><td> 8.083486</td><td> 7.887460</td><td> 8.511340</td><td> 8.068363</td><td> 9.541846</td><td> 8.240755</td><td> 8.110680</td><td>⋯</td><td> 7.918144</td><td> 7.905457</td><td> 9.021824</td><td> 7.460960</td><td> 7.568240</td><td> 8.009838</td><td> 8.800904</td><td> 7.925460</td><td> 8.213657</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioB-M_at </td><td> 7.881876</td><td> 8.326139</td><td> 8.379622</td><td> 7.947517</td><td> 8.778000</td><td> 8.194189</td><td> 9.733321</td><td> 8.328685</td><td> 8.160250</td><td>⋯</td><td> 8.085122</td><td> 8.043810</td><td> 9.252537</td><td> 7.623616</td><td> 7.786479</td><td> 8.064274</td><td> 8.801022</td><td> 8.106188</td><td> 8.350028</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioC-3_at </td><td> 9.431962</td><td> 9.851560</td><td> 9.820634</td><td> 9.527379</td><td>10.084927</td><td> 9.705727</td><td>11.183300</td><td> 9.819315</td><td> 9.778932</td><td>⋯</td><td> 9.391100</td><td> 9.368980</td><td>10.480904</td><td> 9.162528</td><td> 9.168831</td><td> 9.533878</td><td>10.131610</td><td> 9.602517</td><td> 9.958054</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioC-5_at </td><td> 9.206812</td><td> 9.574285</td><td> 9.635055</td><td> 9.314941</td><td> 9.892033</td><td> 9.463784</td><td>11.063276</td><td> 9.590261</td><td> 9.581481</td><td>⋯</td><td> 9.277653</td><td> 9.121820</td><td>10.476738</td><td> 8.957124</td><td> 8.838986</td><td> 9.445741</td><td>10.091344</td><td> 9.542138</td><td> 9.778421</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioD-3_at </td><td>12.162478</td><td>12.573793</td><td>12.512616</td><td>12.359484</td><td>12.618239</td><td>12.420189</td><td>13.420640</td><td>12.614016</td><td>12.419661</td><td>⋯</td><td>11.848646</td><td>11.749284</td><td>12.656261</td><td>11.749962</td><td>11.689466</td><td>11.747099</td><td>12.310210</td><td>12.251810</td><td>12.317773</td><td></td></tr>
	<tr><td>AFFX-r2-Ec-bioD-5_at </td><td>11.331308</td><td>11.971818</td><td>11.833563</td><td>11.651573</td><td>12.006638</td><td>11.749531</td><td>12.938845</td><td>11.890554</td><td>11.768396</td><td>⋯</td><td>11.127818</td><td>11.110795</td><td>12.115626</td><td>10.965360</td><td>10.946355</td><td>11.252270</td><td>11.951325</td><td>11.528131</td><td>11.541979</td><td></td></tr>
	<tr><td>AFFX-r2-P1-cre-3_at  </td><td>13.356565</td><td>13.742737</td><td>13.620602</td><td>13.627280</td><td>13.871845</td><td>13.617009</td><td>14.053436</td><td>13.672234</td><td>13.717266</td><td>⋯</td><td>13.187270</td><td>13.222558</td><td>13.313591</td><td>13.062346</td><td>13.151017</td><td>13.191377</td><td>13.551493</td><td>13.515646</td><td>13.289804</td><td></td></tr>
	<tr><td>AFFX-r2-P1-cre-5_at  </td><td>13.203354</td><td>13.624773</td><td>13.467918</td><td>13.442908</td><td>13.722678</td><td>13.442813</td><td>13.980241</td><td>13.500278</td><td>13.493535</td><td>⋯</td><td>13.121904</td><td>12.984264</td><td>13.399894</td><td>13.030261</td><td>12.992787</td><td>13.166106</td><td>13.540354</td><td>13.417514</td><td>13.270340</td><td></td></tr>
	<tr><td>AFFX-ThrX-3_at       </td><td> 4.251202</td><td> 4.535661</td><td> 5.099850</td><td> 4.707099</td><td> 5.458428</td><td> 4.832445</td><td> 7.603611</td><td> 4.769020</td><td> 4.949368</td><td>⋯</td><td> 8.679110</td><td> 8.378815</td><td> 8.133379</td><td> 7.057436</td><td> 8.169890</td><td> 5.066096</td><td> 9.050793</td><td> 4.792563</td><td> 8.631628</td><td></td></tr>
	<tr><td>AFFX-ThrX-5_at       </td><td> 3.655921</td><td> 3.456087</td><td> 3.649160</td><td> 3.574803</td><td> 3.704276</td><td> 3.675291</td><td> 5.081594</td><td> 3.639642</td><td> 3.663042</td><td>⋯</td><td> 5.792305</td><td> 4.747768</td><td> 6.315986</td><td> 4.296591</td><td> 5.529841</td><td> 4.077782</td><td> 7.513978</td><td> 3.328826</td><td> 6.844849</td><td></td></tr>
	<tr><td>AFFX-ThrX-M_at       </td><td> 3.546140</td><td> 3.509505</td><td> 3.805052</td><td> 3.667717</td><td> 4.454777</td><td> 3.823343</td><td> 6.184042</td><td> 3.668257</td><td> 3.787996</td><td>⋯</td><td> 7.448852</td><td> 6.643659</td><td> 7.512262</td><td> 5.488834</td><td> 6.906802</td><td> 3.857880</td><td> 8.220250</td><td> 3.480863</td><td> 7.995241</td><td></td></tr>
	<tr><td>AFFX-TrpnX-3_at      </td><td> 2.472267</td><td> 2.427352</td><td> 2.286778</td><td> 2.396934</td><td> 2.471298</td><td> 2.348854</td><td> 2.566027</td><td> 2.479157</td><td> 2.496064</td><td>⋯</td><td> 2.606253</td><td> 2.564303</td><td> 2.438315</td><td> 2.517019</td><td> 2.476587</td><td> 2.662852</td><td> 2.389300</td><td> 2.476874</td><td> 2.544633</td><td></td></tr>
	<tr><td>AFFX-TrpnX-5_at      </td><td> 3.065608</td><td> 2.796632</td><td> 2.802991</td><td> 2.848606</td><td> 3.183681</td><td> 2.928234</td><td> 2.842875</td><td> 2.750526</td><td> 2.850283</td><td>⋯</td><td> 2.980219</td><td> 2.857410</td><td> 2.843217</td><td> 2.900981</td><td> 2.945798</td><td> 2.897270</td><td> 2.834124</td><td> 3.033088</td><td> 2.978618</td><td></td></tr>
	<tr><td>AFFX-TrpnX-M_at      </td><td> 2.553665</td><td> 2.488507</td><td> 2.506789</td><td> 2.506172</td><td> 2.578257</td><td> 2.535920</td><td> 2.535563</td><td> 2.678342</td><td> 2.505362</td><td>⋯</td><td> 2.852714</td><td> 2.903490</td><td> 2.917171</td><td> 2.784664</td><td> 2.872477</td><td> 3.197873</td><td> 2.716546</td><td> 2.867979</td><td> 2.797712</td><td></td></tr>
</tbody>
</table>




```R
pheno_group  <- phenoData(GSE21545[[1]])@data[,c(2,38:42)]  %>% as_tibble()  %>% magrittr::set_colnames(c('Samples','Ages','Type','Patients','Outcome','Time'))  %>% 
                   mutate(Type = ifelse(Type == 'peripheral blood mononuclear cells', 'PBMC','Plaques'))
```


```R
pheno_group
```


<table class="dataframe">
<caption>A tibble: 223 × 6</caption>
<thead>
	<tr><th scope=col>Samples</th><th scope=col>Ages</th><th scope=col>Type</th><th scope=col>Patients</th><th scope=col>Outcome</th><th scope=col>Time</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM892518</td><td>68</td><td>PBMC</td><td>279</td><td>FALSE</td><td>1738</td></tr>
	<tr><td>GSM892519</td><td>78</td><td>PBMC</td><td>281</td><td>FALSE</td><td>697 </td></tr>
	<tr><td>GSM892520</td><td>84</td><td>PBMC</td><td>282</td><td>TRUE </td><td>2   </td></tr>
	<tr><td>GSM892521</td><td>80</td><td>PBMC</td><td>284</td><td>FALSE</td><td>1767</td></tr>
	<tr><td>GSM892522</td><td>81</td><td>PBMC</td><td>285</td><td>FALSE</td><td>1326</td></tr>
	<tr><td>GSM892523</td><td>50</td><td>PBMC</td><td>286</td><td>FALSE</td><td>1427</td></tr>
	<tr><td>GSM892524</td><td>75</td><td>PBMC</td><td>288</td><td>TRUE </td><td>7   </td></tr>
	<tr><td>GSM892525</td><td>68</td><td>PBMC</td><td>292</td><td>FALSE</td><td>1422</td></tr>
	<tr><td>GSM892526</td><td>73</td><td>PBMC</td><td>296</td><td>FALSE</td><td>1320</td></tr>
	<tr><td>GSM892527</td><td>75</td><td>PBMC</td><td>298</td><td>FALSE</td><td>689 </td></tr>
	<tr><td>GSM892528</td><td>64</td><td>PBMC</td><td>301</td><td>FALSE</td><td>676 </td></tr>
	<tr><td>GSM892529</td><td>66</td><td>PBMC</td><td>302</td><td>FALSE</td><td>1710</td></tr>
	<tr><td>GSM892530</td><td>71</td><td>PBMC</td><td>306</td><td>FALSE</td><td>1688</td></tr>
	<tr><td>GSM892531</td><td>56</td><td>PBMC</td><td>311</td><td>FALSE</td><td>1711</td></tr>
	<tr><td>GSM892532</td><td>76</td><td>PBMC</td><td>313</td><td>FALSE</td><td>751 </td></tr>
	<tr><td>GSM892533</td><td>78</td><td>PBMC</td><td>318</td><td>FALSE</td><td>1680</td></tr>
	<tr><td>GSM892534</td><td>69</td><td>PBMC</td><td>327</td><td>FALSE</td><td>1876</td></tr>
	<tr><td>GSM892535</td><td>78</td><td>PBMC</td><td>333</td><td>FALSE</td><td>627 </td></tr>
	<tr><td>GSM892536</td><td>66</td><td>PBMC</td><td>334</td><td>FALSE</td><td>1779</td></tr>
	<tr><td>GSM892537</td><td>84</td><td>PBMC</td><td>337</td><td>FALSE</td><td>1542</td></tr>
	<tr><td>GSM892538</td><td>79</td><td>PBMC</td><td>340</td><td>FALSE</td><td>777 </td></tr>
	<tr><td>GSM892539</td><td>78</td><td>PBMC</td><td>345</td><td>FALSE</td><td>687 </td></tr>
	<tr><td>GSM892540</td><td>65</td><td>PBMC</td><td>355</td><td>FALSE</td><td>739 </td></tr>
	<tr><td>GSM892541</td><td>77</td><td>PBMC</td><td>356</td><td>TRUE </td><td>2   </td></tr>
	<tr><td>GSM892542</td><td>79</td><td>PBMC</td><td>359</td><td>FALSE</td><td>1674</td></tr>
	<tr><td>GSM892543</td><td>72</td><td>PBMC</td><td>363</td><td>FALSE</td><td>1649</td></tr>
	<tr><td>GSM892544</td><td>63</td><td>PBMC</td><td>364</td><td>FALSE</td><td>1842</td></tr>
	<tr><td>GSM892545</td><td>81</td><td>PBMC</td><td>365</td><td>TRUE </td><td>922 </td></tr>
	<tr><td>GSM892546</td><td>59</td><td>PBMC</td><td>366</td><td>TRUE </td><td>1   </td></tr>
	<tr><td>GSM892547</td><td>69</td><td>PBMC</td><td>367</td><td>FALSE</td><td>1886</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>GSM892711</td><td>61</td><td>Plaques</td><td>587</td><td>FALSE</td><td>2083</td></tr>
	<tr><td>GSM892712</td><td>48</td><td>Plaques</td><td>588</td><td>FALSE</td><td>1704</td></tr>
	<tr><td>GSM892713</td><td>77</td><td>Plaques</td><td>591</td><td>TRUE </td><td>21  </td></tr>
	<tr><td>GSM892714</td><td>66</td><td>Plaques</td><td>597</td><td>FALSE</td><td>2115</td></tr>
	<tr><td>GSM892715</td><td>73</td><td>Plaques</td><td>603</td><td>FALSE</td><td>1306</td></tr>
	<tr><td>GSM892716</td><td>77</td><td>Plaques</td><td>605</td><td>FALSE</td><td>1870</td></tr>
	<tr><td>GSM892717</td><td>62</td><td>Plaques</td><td>606</td><td>FALSE</td><td>1173</td></tr>
	<tr><td>GSM892718</td><td>82</td><td>Plaques</td><td>618</td><td>FALSE</td><td>395 </td></tr>
	<tr><td>GSM892719</td><td>80</td><td>Plaques</td><td>619</td><td>FALSE</td><td>2619</td></tr>
	<tr><td>GSM892720</td><td>68</td><td>Plaques</td><td>624</td><td>FALSE</td><td>1752</td></tr>
	<tr><td>GSM892721</td><td>70</td><td>Plaques</td><td>625</td><td>FALSE</td><td>1207</td></tr>
	<tr><td>GSM892722</td><td>70</td><td>Plaques</td><td>627</td><td>FALSE</td><td>757 </td></tr>
	<tr><td>GSM892723</td><td>69</td><td>Plaques</td><td>629</td><td>FALSE</td><td>2888</td></tr>
	<tr><td>GSM892724</td><td>71</td><td>Plaques</td><td>630</td><td>TRUE </td><td>423 </td></tr>
	<tr><td>GSM892725</td><td>68</td><td>Plaques</td><td>635</td><td>TRUE </td><td>122 </td></tr>
	<tr><td>GSM892726</td><td>49</td><td>Plaques</td><td>638</td><td>FALSE</td><td>1890</td></tr>
	<tr><td>GSM892727</td><td>76</td><td>Plaques</td><td>642</td><td>FALSE</td><td>1907</td></tr>
	<tr><td>GSM892728</td><td>79</td><td>Plaques</td><td>645</td><td>FALSE</td><td>2397</td></tr>
	<tr><td>GSM892729</td><td>82</td><td>Plaques</td><td>646</td><td>FALSE</td><td>686 </td></tr>
	<tr><td>GSM892730</td><td>81</td><td>Plaques</td><td>651</td><td>FALSE</td><td>2997</td></tr>
	<tr><td>GSM892731</td><td>NA</td><td>Plaques</td><td>654</td><td>FALSE</td><td>652 </td></tr>
	<tr><td>GSM892732</td><td>81</td><td>Plaques</td><td>659</td><td>TRUE </td><td>796 </td></tr>
	<tr><td>GSM892733</td><td>51</td><td>Plaques</td><td>664</td><td>FALSE</td><td>1186</td></tr>
	<tr><td>GSM892734</td><td>82</td><td>Plaques</td><td>669</td><td>TRUE </td><td>516 </td></tr>
	<tr><td>GSM892735</td><td>56</td><td>Plaques</td><td>675</td><td>FALSE</td><td>2223</td></tr>
	<tr><td>GSM892736</td><td>60</td><td>Plaques</td><td>677</td><td>FALSE</td><td>1347</td></tr>
	<tr><td>GSM892737</td><td>79</td><td>Plaques</td><td>678</td><td>FALSE</td><td>2123</td></tr>
	<tr><td>GSM892738</td><td>69</td><td>Plaques</td><td>680</td><td>FALSE</td><td>1792</td></tr>
	<tr><td>GSM892739</td><td>85</td><td>Plaques</td><td>681</td><td>TRUE </td><td>103 </td></tr>
	<tr><td>GSM892740</td><td>73</td><td>Plaques</td><td>686</td><td>FALSE</td><td>1479</td></tr>
</tbody>
</table>




```R
sample_name <- colnames(exprs_matrix %>% 
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
<ol class=list-inline><li>'GSM892518'</li><li>'GSM892519'</li><li>'GSM892520'</li><li>'GSM892521'</li><li>'GSM892522'</li><li>'GSM892523'</li><li>'GSM892524'</li><li>'GSM892525'</li><li>'GSM892526'</li><li>'GSM892527'</li><li>'GSM892528'</li><li>'GSM892529'</li><li>'GSM892530'</li><li>'GSM892531'</li><li>'GSM892532'</li><li>'GSM892533'</li><li>'GSM892534'</li><li>'GSM892535'</li><li>'GSM892536'</li><li>'GSM892537'</li><li>'GSM892538'</li><li>'GSM892539'</li><li>'GSM892540'</li><li>'GSM892541'</li><li>'GSM892542'</li><li>'GSM892543'</li><li>'GSM892544'</li><li>'GSM892545'</li><li>'GSM892546'</li><li>'GSM892547'</li><li>'GSM892548'</li><li>'GSM892549'</li><li>'GSM892550'</li><li>'GSM892551'</li><li>'GSM892552'</li><li>'GSM892553'</li><li>'GSM892554'</li><li>'GSM892555'</li><li>'GSM892556'</li><li>'GSM892557'</li><li>'GSM892558'</li><li>'GSM892559'</li><li>'GSM892560'</li><li>'GSM892561'</li><li>'GSM892562'</li><li>'GSM892563'</li><li>'GSM892564'</li><li>'GSM892565'</li><li>'GSM892566'</li><li>'GSM892567'</li><li>'GSM892568'</li><li>'GSM892569'</li><li>'GSM892570'</li><li>'GSM892571'</li><li>'GSM892572'</li><li>'GSM892573'</li><li>'GSM892574'</li><li>'GSM892575'</li><li>'GSM892576'</li><li>'GSM892577'</li><li>'GSM892578'</li><li>'GSM892579'</li><li>'GSM892580'</li><li>'GSM892581'</li><li>'GSM892582'</li><li>'GSM892583'</li><li>'GSM892584'</li><li>'GSM892585'</li><li>'GSM892586'</li><li>'GSM892587'</li><li>'GSM892588'</li><li>'GSM892589'</li><li>'GSM892590'</li><li>'GSM892591'</li><li>'GSM892592'</li><li>'GSM892593'</li><li>'GSM892594'</li><li>'GSM892595'</li><li>'GSM892596'</li><li>'GSM892597'</li><li>'GSM892598'</li><li>'GSM892599'</li><li>'GSM892600'</li><li>'GSM892601'</li><li>'GSM892602'</li><li>'GSM892603'</li><li>'GSM892604'</li><li>'GSM892605'</li><li>'GSM892606'</li><li>'GSM892607'</li><li>'GSM892608'</li><li>'GSM892609'</li><li>'GSM892610'</li><li>'GSM892611'</li><li>'GSM892612'</li><li>'GSM892613'</li><li>'GSM892614'</li><li>'GSM892615'</li><li>'GSM892616'</li><li>'GSM892617'</li><li>'GSM892618'</li><li>'GSM892619'</li><li>'GSM892620'</li><li>'GSM892621'</li><li>'GSM892622'</li><li>'GSM892623'</li><li>'GSM892624'</li><li>'GSM892625'</li><li>'GSM892626'</li><li>'GSM892627'</li><li>'GSM892628'</li><li>'GSM892629'</li><li>'GSM892630'</li><li>'GSM892631'</li><li>'GSM892632'</li><li>'GSM892633'</li><li>'GSM892634'</li><li>'GSM892635'</li><li>'GSM892636'</li><li>'GSM892637'</li><li>'GSM892638'</li><li>'GSM892639'</li><li>'GSM892640'</li><li>'GSM892641'</li><li>'GSM892642'</li><li>'GSM892643'</li><li>'GSM892644'</li><li>'GSM892645'</li><li>'GSM892646'</li><li>'GSM892647'</li><li>'GSM892648'</li><li>'GSM892649'</li><li>'GSM892650'</li><li>'GSM892651'</li><li>'GSM892652'</li><li>'GSM892653'</li><li>'GSM892654'</li><li>'GSM892655'</li><li>'GSM892656'</li><li>'GSM892657'</li><li>'GSM892658'</li><li>'GSM892659'</li><li>'GSM892660'</li><li>'GSM892661'</li><li>'GSM892662'</li><li>'GSM892663'</li><li>'GSM892664'</li><li>'GSM892665'</li><li>'GSM892666'</li><li>'GSM892667'</li><li>'GSM892668'</li><li>'GSM892669'</li><li>'GSM892670'</li><li>'GSM892671'</li><li>'GSM892672'</li><li>'GSM892673'</li><li>'GSM892674'</li><li>'GSM892675'</li><li>'GSM892676'</li><li>'GSM892677'</li><li>'GSM892678'</li><li>'GSM892679'</li><li>'GSM892680'</li><li>'GSM892681'</li><li>'GSM892682'</li><li>'GSM892683'</li><li>'GSM892684'</li><li>'GSM892685'</li><li>'GSM892686'</li><li>'GSM892687'</li><li>'GSM892688'</li><li>'GSM892689'</li><li>'GSM892690'</li><li>'GSM892691'</li><li>'GSM892692'</li><li>'GSM892693'</li><li>'GSM892694'</li><li>'GSM892695'</li><li>'GSM892696'</li><li>'GSM892697'</li><li>'GSM892698'</li><li>'GSM892699'</li><li>'GSM892700'</li><li>'GSM892701'</li><li>'GSM892702'</li><li>'GSM892703'</li><li>'GSM892704'</li><li>'GSM892705'</li><li>'GSM892706'</li><li>'GSM892707'</li><li>'GSM892708'</li><li>'GSM892709'</li><li>'GSM892710'</li><li>'GSM892711'</li><li>'GSM892712'</li><li>'GSM892713'</li><li>'GSM892714'</li><li>'GSM892715'</li><li>'GSM892716'</li><li>'GSM892717'</li><li>'GSM892718'</li><li>'GSM892719'</li><li>'GSM892720'</li><li>'GSM892721'</li><li>'GSM892722'</li><li>'GSM892723'</li><li>'GSM892724'</li><li>'GSM892725'</li><li>'GSM892726'</li><li>'GSM892727'</li><li>'GSM892728'</li><li>'GSM892729'</li><li>'GSM892730'</li><li>'GSM892731'</li><li>'GSM892732'</li><li>'GSM892733'</li><li>'GSM892734'</li><li>'GSM892735'</li><li>'GSM892736'</li><li>'GSM892737'</li><li>'GSM892738'</li><li>'GSM892739'</li><li>'GSM892740'</li></ol>




```R
exprs_matrix_before  <-  exprs_matrix  %>% 
 as.data.frame()  %>% 
  rownames_to_column(var = 'Probe')  %>%  
   left_join(ids, by = c('Probe'= 'ID')) %>% 
    filter(!is.na(symbol))   %>% 
     filter(symbol != '')   %>% 
     .[,2:ncol(.)]   %>% 
       group_by(symbol) 
```


```R
####get the formula
Formula  <- exprs_matrix_before  %>% 
   .[,1:(ncol(.)-1)]  %>% 
     colnames(.) %>% 
      tibble() %>% 
       magrittr::set_colnames(c('Samples')) %>% 
        mutate(Samples2 = Samples) %>% 
         mutate(Formula = paste0(Samples2, ' = mean(', Samples,'), '))  %>% 
           pull(Formula)
```


```R
Formula
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'GSM892518 = mean(GSM892518), '</li><li>'GSM892519 = mean(GSM892519), '</li><li>'GSM892520 = mean(GSM892520), '</li><li>'GSM892521 = mean(GSM892521), '</li><li>'GSM892522 = mean(GSM892522), '</li><li>'GSM892523 = mean(GSM892523), '</li><li>'GSM892524 = mean(GSM892524), '</li><li>'GSM892525 = mean(GSM892525), '</li><li>'GSM892526 = mean(GSM892526), '</li><li>'GSM892527 = mean(GSM892527), '</li><li>'GSM892528 = mean(GSM892528), '</li><li>'GSM892529 = mean(GSM892529), '</li><li>'GSM892530 = mean(GSM892530), '</li><li>'GSM892531 = mean(GSM892531), '</li><li>'GSM892532 = mean(GSM892532), '</li><li>'GSM892533 = mean(GSM892533), '</li><li>'GSM892534 = mean(GSM892534), '</li><li>'GSM892535 = mean(GSM892535), '</li><li>'GSM892536 = mean(GSM892536), '</li><li>'GSM892537 = mean(GSM892537), '</li><li>'GSM892538 = mean(GSM892538), '</li><li>'GSM892539 = mean(GSM892539), '</li><li>'GSM892540 = mean(GSM892540), '</li><li>'GSM892541 = mean(GSM892541), '</li><li>'GSM892542 = mean(GSM892542), '</li><li>'GSM892543 = mean(GSM892543), '</li><li>'GSM892544 = mean(GSM892544), '</li><li>'GSM892545 = mean(GSM892545), '</li><li>'GSM892546 = mean(GSM892546), '</li><li>'GSM892547 = mean(GSM892547), '</li><li>'GSM892548 = mean(GSM892548), '</li><li>'GSM892549 = mean(GSM892549), '</li><li>'GSM892550 = mean(GSM892550), '</li><li>'GSM892551 = mean(GSM892551), '</li><li>'GSM892552 = mean(GSM892552), '</li><li>'GSM892553 = mean(GSM892553), '</li><li>'GSM892554 = mean(GSM892554), '</li><li>'GSM892555 = mean(GSM892555), '</li><li>'GSM892556 = mean(GSM892556), '</li><li>'GSM892557 = mean(GSM892557), '</li><li>'GSM892558 = mean(GSM892558), '</li><li>'GSM892559 = mean(GSM892559), '</li><li>'GSM892560 = mean(GSM892560), '</li><li>'GSM892561 = mean(GSM892561), '</li><li>'GSM892562 = mean(GSM892562), '</li><li>'GSM892563 = mean(GSM892563), '</li><li>'GSM892564 = mean(GSM892564), '</li><li>'GSM892565 = mean(GSM892565), '</li><li>'GSM892566 = mean(GSM892566), '</li><li>'GSM892567 = mean(GSM892567), '</li><li>'GSM892568 = mean(GSM892568), '</li><li>'GSM892569 = mean(GSM892569), '</li><li>'GSM892570 = mean(GSM892570), '</li><li>'GSM892571 = mean(GSM892571), '</li><li>'GSM892572 = mean(GSM892572), '</li><li>'GSM892573 = mean(GSM892573), '</li><li>'GSM892574 = mean(GSM892574), '</li><li>'GSM892575 = mean(GSM892575), '</li><li>'GSM892576 = mean(GSM892576), '</li><li>'GSM892577 = mean(GSM892577), '</li><li>'GSM892578 = mean(GSM892578), '</li><li>'GSM892579 = mean(GSM892579), '</li><li>'GSM892580 = mean(GSM892580), '</li><li>'GSM892581 = mean(GSM892581), '</li><li>'GSM892582 = mean(GSM892582), '</li><li>'GSM892583 = mean(GSM892583), '</li><li>'GSM892584 = mean(GSM892584), '</li><li>'GSM892585 = mean(GSM892585), '</li><li>'GSM892586 = mean(GSM892586), '</li><li>'GSM892587 = mean(GSM892587), '</li><li>'GSM892588 = mean(GSM892588), '</li><li>'GSM892589 = mean(GSM892589), '</li><li>'GSM892590 = mean(GSM892590), '</li><li>'GSM892591 = mean(GSM892591), '</li><li>'GSM892592 = mean(GSM892592), '</li><li>'GSM892593 = mean(GSM892593), '</li><li>'GSM892594 = mean(GSM892594), '</li><li>'GSM892595 = mean(GSM892595), '</li><li>'GSM892596 = mean(GSM892596), '</li><li>'GSM892597 = mean(GSM892597), '</li><li>'GSM892598 = mean(GSM892598), '</li><li>'GSM892599 = mean(GSM892599), '</li><li>'GSM892600 = mean(GSM892600), '</li><li>'GSM892601 = mean(GSM892601), '</li><li>'GSM892602 = mean(GSM892602), '</li><li>'GSM892603 = mean(GSM892603), '</li><li>'GSM892604 = mean(GSM892604), '</li><li>'GSM892605 = mean(GSM892605), '</li><li>'GSM892606 = mean(GSM892606), '</li><li>'GSM892607 = mean(GSM892607), '</li><li>'GSM892608 = mean(GSM892608), '</li><li>'GSM892609 = mean(GSM892609), '</li><li>'GSM892610 = mean(GSM892610), '</li><li>'GSM892611 = mean(GSM892611), '</li><li>'GSM892612 = mean(GSM892612), '</li><li>'GSM892613 = mean(GSM892613), '</li><li>'GSM892614 = mean(GSM892614), '</li><li>'GSM892615 = mean(GSM892615), '</li><li>'GSM892616 = mean(GSM892616), '</li><li>'GSM892617 = mean(GSM892617), '</li><li>'GSM892618 = mean(GSM892618), '</li><li>'GSM892619 = mean(GSM892619), '</li><li>'GSM892620 = mean(GSM892620), '</li><li>'GSM892621 = mean(GSM892621), '</li><li>'GSM892622 = mean(GSM892622), '</li><li>'GSM892623 = mean(GSM892623), '</li><li>'GSM892624 = mean(GSM892624), '</li><li>'GSM892625 = mean(GSM892625), '</li><li>'GSM892626 = mean(GSM892626), '</li><li>'GSM892627 = mean(GSM892627), '</li><li>'GSM892628 = mean(GSM892628), '</li><li>'GSM892629 = mean(GSM892629), '</li><li>'GSM892630 = mean(GSM892630), '</li><li>'GSM892631 = mean(GSM892631), '</li><li>'GSM892632 = mean(GSM892632), '</li><li>'GSM892633 = mean(GSM892633), '</li><li>'GSM892634 = mean(GSM892634), '</li><li>'GSM892635 = mean(GSM892635), '</li><li>'GSM892636 = mean(GSM892636), '</li><li>'GSM892637 = mean(GSM892637), '</li><li>'GSM892638 = mean(GSM892638), '</li><li>'GSM892639 = mean(GSM892639), '</li><li>'GSM892640 = mean(GSM892640), '</li><li>'GSM892641 = mean(GSM892641), '</li><li>'GSM892642 = mean(GSM892642), '</li><li>'GSM892643 = mean(GSM892643), '</li><li>'GSM892644 = mean(GSM892644), '</li><li>'GSM892645 = mean(GSM892645), '</li><li>'GSM892646 = mean(GSM892646), '</li><li>'GSM892647 = mean(GSM892647), '</li><li>'GSM892648 = mean(GSM892648), '</li><li>'GSM892649 = mean(GSM892649), '</li><li>'GSM892650 = mean(GSM892650), '</li><li>'GSM892651 = mean(GSM892651), '</li><li>'GSM892652 = mean(GSM892652), '</li><li>'GSM892653 = mean(GSM892653), '</li><li>'GSM892654 = mean(GSM892654), '</li><li>'GSM892655 = mean(GSM892655), '</li><li>'GSM892656 = mean(GSM892656), '</li><li>'GSM892657 = mean(GSM892657), '</li><li>'GSM892658 = mean(GSM892658), '</li><li>'GSM892659 = mean(GSM892659), '</li><li>'GSM892660 = mean(GSM892660), '</li><li>'GSM892661 = mean(GSM892661), '</li><li>'GSM892662 = mean(GSM892662), '</li><li>'GSM892663 = mean(GSM892663), '</li><li>'GSM892664 = mean(GSM892664), '</li><li>'GSM892665 = mean(GSM892665), '</li><li>'GSM892666 = mean(GSM892666), '</li><li>'GSM892667 = mean(GSM892667), '</li><li>'GSM892668 = mean(GSM892668), '</li><li>'GSM892669 = mean(GSM892669), '</li><li>'GSM892670 = mean(GSM892670), '</li><li>'GSM892671 = mean(GSM892671), '</li><li>'GSM892672 = mean(GSM892672), '</li><li>'GSM892673 = mean(GSM892673), '</li><li>'GSM892674 = mean(GSM892674), '</li><li>'GSM892675 = mean(GSM892675), '</li><li>'GSM892676 = mean(GSM892676), '</li><li>'GSM892677 = mean(GSM892677), '</li><li>'GSM892678 = mean(GSM892678), '</li><li>'GSM892679 = mean(GSM892679), '</li><li>'GSM892680 = mean(GSM892680), '</li><li>'GSM892681 = mean(GSM892681), '</li><li>'GSM892682 = mean(GSM892682), '</li><li>'GSM892683 = mean(GSM892683), '</li><li>'GSM892684 = mean(GSM892684), '</li><li>'GSM892685 = mean(GSM892685), '</li><li>'GSM892686 = mean(GSM892686), '</li><li>'GSM892687 = mean(GSM892687), '</li><li>'GSM892688 = mean(GSM892688), '</li><li>'GSM892689 = mean(GSM892689), '</li><li>'GSM892690 = mean(GSM892690), '</li><li>'GSM892691 = mean(GSM892691), '</li><li>'GSM892692 = mean(GSM892692), '</li><li>'GSM892693 = mean(GSM892693), '</li><li>'GSM892694 = mean(GSM892694), '</li><li>'GSM892695 = mean(GSM892695), '</li><li>'GSM892696 = mean(GSM892696), '</li><li>'GSM892697 = mean(GSM892697), '</li><li>'GSM892698 = mean(GSM892698), '</li><li>'GSM892699 = mean(GSM892699), '</li><li>'GSM892700 = mean(GSM892700), '</li><li>'GSM892701 = mean(GSM892701), '</li><li>'GSM892702 = mean(GSM892702), '</li><li>'GSM892703 = mean(GSM892703), '</li><li>'GSM892704 = mean(GSM892704), '</li><li>'GSM892705 = mean(GSM892705), '</li><li>'GSM892706 = mean(GSM892706), '</li><li>'GSM892707 = mean(GSM892707), '</li><li>'GSM892708 = mean(GSM892708), '</li><li>'GSM892709 = mean(GSM892709), '</li><li>'GSM892710 = mean(GSM892710), '</li><li>'GSM892711 = mean(GSM892711), '</li><li>'GSM892712 = mean(GSM892712), '</li><li>'GSM892713 = mean(GSM892713), '</li><li>'GSM892714 = mean(GSM892714), '</li><li>'GSM892715 = mean(GSM892715), '</li><li>'GSM892716 = mean(GSM892716), '</li><li>'GSM892717 = mean(GSM892717), '</li><li>'GSM892718 = mean(GSM892718), '</li><li>'GSM892719 = mean(GSM892719), '</li><li>'GSM892720 = mean(GSM892720), '</li><li>'GSM892721 = mean(GSM892721), '</li><li>'GSM892722 = mean(GSM892722), '</li><li>'GSM892723 = mean(GSM892723), '</li><li>'GSM892724 = mean(GSM892724), '</li><li>'GSM892725 = mean(GSM892725), '</li><li>'GSM892726 = mean(GSM892726), '</li><li>'GSM892727 = mean(GSM892727), '</li><li>'GSM892728 = mean(GSM892728), '</li><li>'GSM892729 = mean(GSM892729), '</li><li>'GSM892730 = mean(GSM892730), '</li><li>'GSM892731 = mean(GSM892731), '</li><li>'GSM892732 = mean(GSM892732), '</li><li>'GSM892733 = mean(GSM892733), '</li><li>'GSM892734 = mean(GSM892734), '</li><li>'GSM892735 = mean(GSM892735), '</li><li>'GSM892736 = mean(GSM892736), '</li><li>'GSM892737 = mean(GSM892737), '</li><li>'GSM892738 = mean(GSM892738), '</li><li>'GSM892739 = mean(GSM892739), '</li><li>'GSM892740 = mean(GSM892740), '</li></ol>




```R
normalized_gset_mean <- exprs_matrix_before  %>% 
              summarise(
              GSM892518 = mean(GSM892518), GSM892519 = mean(GSM892519), GSM892520 = mean(GSM892520), GSM892521 = mean(GSM892521), 
              GSM892522 = mean(GSM892522), GSM892523 = mean(GSM892523), GSM892524 = mean(GSM892524), GSM892525 = mean(GSM892525), 
              GSM892526 = mean(GSM892526), GSM892527 = mean(GSM892527), GSM892528 = mean(GSM892528), GSM892529 = mean(GSM892529), 
              GSM892530 = mean(GSM892530), GSM892531 = mean(GSM892531), GSM892532 = mean(GSM892532), GSM892533 = mean(GSM892533), 
              GSM892534 = mean(GSM892534), GSM892535 = mean(GSM892535), GSM892536 = mean(GSM892536), GSM892537 = mean(GSM892537), 
              GSM892538 = mean(GSM892538), GSM892539 = mean(GSM892539), GSM892540 = mean(GSM892540), GSM892541 = mean(GSM892541),
              GSM892542 = mean(GSM892542), GSM892543 = mean(GSM892543), GSM892544 = mean(GSM892544), GSM892545 = mean(GSM892545), 
              GSM892546 = mean(GSM892546), GSM892547 = mean(GSM892547), GSM892548 = mean(GSM892548), GSM892549 = mean(GSM892549), 
              GSM892550 = mean(GSM892550), GSM892551 = mean(GSM892551), GSM892552 = mean(GSM892552), GSM892553 = mean(GSM892553), 
              GSM892554 = mean(GSM892554), GSM892555 = mean(GSM892555), GSM892556 = mean(GSM892556), GSM892557 = mean(GSM892557), 
              GSM892558 = mean(GSM892558), GSM892559 = mean(GSM892559), GSM892560 = mean(GSM892560), GSM892561 = mean(GSM892561), 
              GSM892562 = mean(GSM892562), GSM892563 = mean(GSM892563), GSM892564 = mean(GSM892564), GSM892565 = mean(GSM892565), 
              GSM892566 = mean(GSM892566), GSM892567 = mean(GSM892567), GSM892568 = mean(GSM892568), GSM892569 = mean(GSM892569),
              GSM892570 = mean(GSM892570), GSM892571 = mean(GSM892571), GSM892572 = mean(GSM892572), GSM892573 = mean(GSM892573), 
              GSM892574 = mean(GSM892574), GSM892575 = mean(GSM892575), GSM892576 = mean(GSM892576), GSM892577 = mean(GSM892577), 
              GSM892578 = mean(GSM892578), GSM892579 = mean(GSM892579), GSM892580 = mean(GSM892580), GSM892581 = mean(GSM892581),
              GSM892582 = mean(GSM892582), GSM892583 = mean(GSM892583), GSM892584 = mean(GSM892584), GSM892585 = mean(GSM892585),
              GSM892586 = mean(GSM892586), GSM892587 = mean(GSM892587), GSM892588 = mean(GSM892588), GSM892589 = mean(GSM892589), 
              GSM892590 = mean(GSM892590), GSM892591 = mean(GSM892591), GSM892592 = mean(GSM892592), GSM892593 = mean(GSM892593),
              GSM892594 = mean(GSM892594), GSM892595 = mean(GSM892595), GSM892596 = mean(GSM892596), GSM892597 = mean(GSM892597), 
              GSM892598 = mean(GSM892598), GSM892599 = mean(GSM892599), GSM892600 = mean(GSM892600), GSM892601 = mean(GSM892601),
              GSM892602 = mean(GSM892602), GSM892603 = mean(GSM892603), GSM892604 = mean(GSM892604), GSM892605 = mean(GSM892605), 
              GSM892606 = mean(GSM892606), GSM892607 = mean(GSM892607), GSM892608 = mean(GSM892608), GSM892609 = mean(GSM892609), 
              GSM892610 = mean(GSM892610), GSM892611 = mean(GSM892611), GSM892612 = mean(GSM892612), GSM892613 = mean(GSM892613),
              GSM892614 = mean(GSM892614), GSM892615 = mean(GSM892615), GSM892616 = mean(GSM892616), GSM892617 = mean(GSM892617), 
              GSM892618 = mean(GSM892618), GSM892619 = mean(GSM892619), GSM892620 = mean(GSM892620), GSM892621 = mean(GSM892621), 
              GSM892622 = mean(GSM892622), GSM892623 = mean(GSM892623), GSM892624 = mean(GSM892624), GSM892625 = mean(GSM892625), 
              GSM892626 = mean(GSM892626), GSM892627 = mean(GSM892627), GSM892628 = mean(GSM892628), GSM892629 = mean(GSM892629), 
              GSM892630 = mean(GSM892630), GSM892631 = mean(GSM892631), GSM892632 = mean(GSM892632), GSM892633 = mean(GSM892633), 
              GSM892634 = mean(GSM892634), GSM892635 = mean(GSM892635), GSM892636 = mean(GSM892636), GSM892637 = mean(GSM892637), 
              GSM892638 = mean(GSM892638), GSM892639 = mean(GSM892639), GSM892640 = mean(GSM892640), GSM892641 = mean(GSM892641), 
              GSM892642 = mean(GSM892642), GSM892643 = mean(GSM892643), GSM892644 = mean(GSM892644), GSM892645 = mean(GSM892645), 
              GSM892646 = mean(GSM892646), GSM892647 = mean(GSM892647), GSM892648 = mean(GSM892648), GSM892649 = mean(GSM892649), 
              GSM892650 = mean(GSM892650), GSM892651 = mean(GSM892651), GSM892652 = mean(GSM892652), GSM892653 = mean(GSM892653), 
              GSM892654 = mean(GSM892654), GSM892655 = mean(GSM892655), GSM892656 = mean(GSM892656), GSM892657 = mean(GSM892657), 
              GSM892658 = mean(GSM892658), GSM892659 = mean(GSM892659), GSM892660 = mean(GSM892660), GSM892661 = mean(GSM892661), 
              GSM892662 = mean(GSM892662), GSM892663 = mean(GSM892663), GSM892664 = mean(GSM892664), GSM892665 = mean(GSM892665), 
              GSM892666 = mean(GSM892666), GSM892667 = mean(GSM892667), GSM892668 = mean(GSM892668), GSM892669 = mean(GSM892669), 
              GSM892670 = mean(GSM892670), GSM892671 = mean(GSM892671), GSM892672 = mean(GSM892672), GSM892673 = mean(GSM892673),
              GSM892674 = mean(GSM892674), GSM892675 = mean(GSM892675), GSM892676 = mean(GSM892676), GSM892677 = mean(GSM892677), 
              GSM892678 = mean(GSM892678), GSM892679 = mean(GSM892679), GSM892680 = mean(GSM892680), GSM892681 = mean(GSM892681), 
              GSM892682 = mean(GSM892682), GSM892683 = mean(GSM892683), GSM892684 = mean(GSM892684), GSM892685 = mean(GSM892685), 
              GSM892686 = mean(GSM892686), GSM892687 = mean(GSM892687), GSM892688 = mean(GSM892688), GSM892689 = mean(GSM892689), 
              GSM892690 = mean(GSM892690), GSM892691 = mean(GSM892691), GSM892692 = mean(GSM892692), GSM892693 = mean(GSM892693), 
              GSM892694 = mean(GSM892694), GSM892695 = mean(GSM892695), GSM892696 = mean(GSM892696), GSM892697 = mean(GSM892697), 
              GSM892698 = mean(GSM892698), GSM892699 = mean(GSM892699), GSM892700 = mean(GSM892700), GSM892701 = mean(GSM892701), 
              GSM892702 = mean(GSM892702), GSM892703 = mean(GSM892703), GSM892704 = mean(GSM892704), GSM892705 = mean(GSM892705), 
              GSM892706 = mean(GSM892706), GSM892707 = mean(GSM892707), GSM892708 = mean(GSM892708), GSM892709 = mean(GSM892709), 
              GSM892710 = mean(GSM892710), GSM892711 = mean(GSM892711), GSM892712 = mean(GSM892712), GSM892713 = mean(GSM892713), 
              GSM892714 = mean(GSM892714), GSM892715 = mean(GSM892715), GSM892716 = mean(GSM892716), GSM892717 = mean(GSM892717), 
              GSM892718 = mean(GSM892718), GSM892719 = mean(GSM892719), GSM892720 = mean(GSM892720), GSM892721 = mean(GSM892721), 
              GSM892722 = mean(GSM892722), GSM892723 = mean(GSM892723), GSM892724 = mean(GSM892724), GSM892725 = mean(GSM892725), 
              GSM892726 = mean(GSM892726), GSM892727 = mean(GSM892727), GSM892728 = mean(GSM892728), GSM892729 = mean(GSM892729), 
              GSM892730 = mean(GSM892730), GSM892731 = mean(GSM892731), GSM892732 = mean(GSM892732), GSM892733 = mean(GSM892733), 
              GSM892734 = mean(GSM892734), GSM892735 = mean(GSM892735), GSM892736 = mean(GSM892736), GSM892737 = mean(GSM892737), 
              GSM892738 = mean(GSM892738), GSM892739 = mean(GSM892739), GSM892740 = mean(GSM892740))  %>% 
         column_to_rownames(var ='symbol')
```


```R
library(ggpubr)
```


```R
selected_gene <- c('ALDOA','CREG1','LGMN','PKM')
```


```R
normalized_gset_mean[selected_gene,] %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(Yellow_Module,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Yellow_Module", y = "Time", facet.by = 'Type',
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson")
```


<table class="dataframe">
<caption>A data.frame: 4 × 223</caption>
<thead>
	<tr><th></th><th scope=col>GSM892518</th><th scope=col>GSM892519</th><th scope=col>GSM892520</th><th scope=col>GSM892521</th><th scope=col>GSM892522</th><th scope=col>GSM892523</th><th scope=col>GSM892524</th><th scope=col>GSM892525</th><th scope=col>GSM892526</th><th scope=col>GSM892527</th><th scope=col>⋯</th><th scope=col>GSM892731</th><th scope=col>GSM892732</th><th scope=col>GSM892733</th><th scope=col>GSM892734</th><th scope=col>GSM892735</th><th scope=col>GSM892736</th><th scope=col>GSM892737</th><th scope=col>GSM892738</th><th scope=col>GSM892739</th><th scope=col>GSM892740</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ALDOA</th><td>10.886218</td><td>10.734034</td><td>10.829143</td><td>10.735125</td><td>11.052067</td><td>10.584845</td><td>10.623859</td><td>10.756301</td><td>10.789925</td><td>10.275724</td><td>⋯</td><td>10.534107</td><td>11.176459</td><td>10.392333</td><td>10.634499</td><td>10.755860</td><td>10.446059</td><td>10.606249</td><td> 9.932506</td><td>10.967463</td><td>10.792791</td></tr>
	<tr><th scope=row>CREG1</th><td>10.364084</td><td>10.461196</td><td>10.480717</td><td>10.027528</td><td>10.002447</td><td>10.025216</td><td> 9.633990</td><td>10.411462</td><td>10.433035</td><td> 9.931594</td><td>⋯</td><td>11.440414</td><td>12.107784</td><td>11.671501</td><td>11.767611</td><td>11.428542</td><td>11.049492</td><td>11.560494</td><td>11.193442</td><td>11.995219</td><td>11.929663</td></tr>
	<tr><th scope=row>LGMN</th><td> 6.135511</td><td> 6.653872</td><td> 6.514739</td><td> 6.428221</td><td> 5.929549</td><td> 7.215479</td><td> 6.835788</td><td> 6.536581</td><td> 6.590901</td><td> 6.717306</td><td>⋯</td><td>11.868023</td><td>12.253235</td><td>11.432726</td><td>11.382816</td><td>11.239292</td><td>10.671554</td><td>11.435031</td><td>10.837853</td><td>11.944222</td><td>11.630108</td></tr>
	<tr><th scope=row>PKM</th><td> 6.685186</td><td> 7.000903</td><td> 6.911282</td><td> 6.785602</td><td> 6.366404</td><td> 6.512294</td><td> 7.441863</td><td> 6.741107</td><td> 6.862693</td><td> 6.769872</td><td>⋯</td><td> 7.107402</td><td> 9.744618</td><td> 9.916361</td><td> 9.227988</td><td> 9.108028</td><td> 8.927357</td><td> 8.815684</td><td> 7.970550</td><td> 7.613127</td><td> 8.930595</td></tr>
</tbody>
</table>




```R
selected_col  <- c('#DA4C35','#FFB300')
# Only colour strips in x-direction
strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = selected_col))
```


```R
p1 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(ALDOA,Type,Ages) %>% 
        mutate(Ages = as.integer(Ages)) %>% 
           ggscatter(x = "Ages", y = "ALDOA", xlab = 'Ages',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p1

p2 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(CREG1,Type,Ages) %>% 
        mutate(Ages = as.integer(Ages)) %>% 
           ggscatter(x = "Ages", y = "CREG1", xlab = 'Ages',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p2


p3 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(LGMN,Type,Ages) %>% 
        mutate(Ages = as.integer(Ages)) %>% 
           ggscatter(x = "Ages", y = "LGMN", xlab = 'Ages',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p3

p4 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(PKM,Type,Ages) %>% 
        mutate(Ages = as.integer(Ages)) %>% 
           ggscatter(x = "Ages", y = "PKM", xlab = 'Ages',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p4

p0 <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 2,guides='collect')  +  
    plot_annotation(title = 'GSE21545',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE21545_4Gene_Ages.pdf', height = 10, width = 18, units = 'cm')
```

    Warning message in mask$eval_all_mutate(quo):
    “NAs introduced by coercion”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”
    Warning message in mask$eval_all_mutate(quo):
    “NAs introduced by coercion”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_17_1.png)
    


    Warning message in mask$eval_all_mutate(quo):
    “NAs introduced by coercion”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_17_3.png)
    


    Warning message in mask$eval_all_mutate(quo):
    “NAs introduced by coercion”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_17_5.png)
    


    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_17_7.png)
    


    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_smooth()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing non-finite values (`stat_cor()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values (`geom_point()`).”



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_17_9.png)
    



```R
library(patchwork)
```


```R
p1  <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(ALDOA,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "ALDOA", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p2  <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(CREG1,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "CREG1", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p3  <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(LGMN,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "LGMN", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p4  <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(PKM,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "PKM", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p0 <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 4,guides='collect')  +  
    plot_annotation(title = 'GSE21545',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE21545_4Gene_Days.pdf', height = 5, width = 16, units = 'cm')
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_19_1.png)
    



```R

p1.1 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(ALDOA,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "ALDOA", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p1.1 



p2.1 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(CREG1,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "CREG1", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p2.1 

p3.1 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(LGMN,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "LGMN", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p3.1 

p4.1 <- normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(PKM,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "PKM", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2,label.x = 800)+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p4.1 


p0 <- p1.1 + p2.1 + p3.1 + p4.1 + 
  plot_layout(ncol = 2,guides='collect')  +  
    plot_annotation(title = 'GSE21545',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE21545_4Gene_Days_2_Sites.pdf', height = 10, width = 17, units = 'cm')
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_20_1.png)
    


    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_20_3.png)
    


    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_20_5.png)
    


    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_20_7.png)
    


    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_20_9.png)
    



```R
normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(CREG1,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "CREG1", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson")+
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))

normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(CREG1,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "CREG1", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson")+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))

```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_21_1.png)
    



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_21_2.png)
    



```R

normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(LGMN,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "LGMN", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson")+
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))



normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(LGMN,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "LGMN", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson")+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_22_1.png)
    



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_22_2.png)
    



```R

normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(PKM,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "PKM", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson")+
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))

normalized_gset_mean[selected_gene,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() %>% 
        select(PKM,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "PKM", xlab = 'Days',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson")+ ggh4x::facet_wrap2(~Type, scales = "free", strip = strip) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'
    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_23_1.png)
    



    
![png](Step11.5_GSE21545_files/Step11.5_GSE21545_23_2.png)
    

