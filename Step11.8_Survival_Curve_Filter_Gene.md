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
head(Formula)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'GSM892518 = mean(GSM892518), '</li><li>'GSM892519 = mean(GSM892519), '</li><li>'GSM892520 = mean(GSM892520), '</li><li>'GSM892521 = mean(GSM892521), '</li><li>'GSM892522 = mean(GSM892522), '</li><li>'GSM892523 = mean(GSM892523), '</li></ol>




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
selected_gene <- c('CCL23','DUSP5','IL31RA','UPP1')
```


```R
normalized_gset_mean_tibble <- normalized_gset_mean    %>% 
   rownames_to_column(var = 'Gene') %>% 
   pivot_longer(-Gene)
```


```R
normalized_gset_mean_tibble  %>% 
 group_by(Gene) %>% 
   summarise(value_mean = mean(value))  %>% 
    arrange(-value_mean) %>% 
     pull(value_mean) %>% 
      quantile()
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>0%</dt><dd>2.21321157070852</dd><dt>25%</dt><dd>3.78574013311659</dd><dt>50%</dt><dd>5.03609850201794</dd><dt>75%</dt><dd>6.52484728306745</dd><dt>100%</dt><dd>13.9613926591031</dd></dl>




```R
selected_gene  <- normalized_gset_mean_tibble  %>% 
 group_by(Gene) %>% 
   summarise(value_mean = mean(value)) %>% 
    filter(value_mean > 6.5) %>% 
      pull(Gene)
```


```R
normalized_gset_mean_high_expression  <- normalized_gset_mean_tibble %>% 
  filter(Gene %in% selected_gene) %>% 
    pivot_wider(names_from = 'name', values_from = 'value') %>% 
     column_to_rownames(var = 'Gene')
```


```R
normalized_gset_mean_high_expression %>% 
  t()  %>% 
   as.data.frame()  
```


<table class="dataframe">
<caption>A data.frame: 223 × 5958</caption>
<thead>
	<tr><th></th><th scope=col>AACS</th><th scope=col>AAED1</th><th scope=col>AAMP</th><th scope=col>AAR2</th><th scope=col>AARS</th><th scope=col>AARSD1 /// PTGES3L /// PTGES3L-AARSD1</th><th scope=col>AASDH</th><th scope=col>AATF</th><th scope=col>AATK</th><th scope=col>ABCB10</th><th scope=col>⋯</th><th scope=col>ZSCAN21</th><th scope=col>ZSCAN29</th><th scope=col>ZSWIM6</th><th scope=col>ZSWIM7</th><th scope=col>ZSWIM8</th><th scope=col>ZUFSP</th><th scope=col>ZW10</th><th scope=col>ZYG11B</th><th scope=col>ZYX</th><th scope=col>ZZZ3</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM892518</th><td>6.744481</td><td>8.098326</td><td>8.845298</td><td>8.027744</td><td>8.133400</td><td>6.779367</td><td>7.192553</td><td>8.788231</td><td>7.296900</td><td>8.386076</td><td>⋯</td><td>7.154543</td><td>7.512841</td><td> 9.404165</td><td>7.772632</td><td>6.273566</td><td>7.424264</td><td>7.443882</td><td>7.309862</td><td>7.839783</td><td>7.745338</td></tr>
	<tr><th scope=row>GSM892519</th><td>6.606717</td><td>7.951058</td><td>8.353823</td><td>7.704787</td><td>7.991139</td><td>6.713854</td><td>7.712399</td><td>8.498978</td><td>6.748701</td><td>8.571861</td><td>⋯</td><td>7.078088</td><td>7.612578</td><td> 9.580369</td><td>7.768381</td><td>6.577481</td><td>7.568137</td><td>7.431788</td><td>7.383383</td><td>7.792338</td><td>7.728572</td></tr>
	<tr><th scope=row>GSM892520</th><td>6.703049</td><td>8.253258</td><td>8.821917</td><td>7.935994</td><td>8.121270</td><td>6.987739</td><td>7.483351</td><td>8.788333</td><td>7.220110</td><td>8.360093</td><td>⋯</td><td>7.195625</td><td>7.729012</td><td> 9.391700</td><td>7.845753</td><td>6.349280</td><td>7.342488</td><td>7.437652</td><td>7.422046</td><td>7.996883</td><td>7.661320</td></tr>
	<tr><th scope=row>GSM892521</th><td>6.817530</td><td>8.024947</td><td>8.754522</td><td>7.956263</td><td>8.256909</td><td>6.871659</td><td>7.263732</td><td>8.747234</td><td>7.395759</td><td>8.485493</td><td>⋯</td><td>6.971817</td><td>7.460289</td><td> 9.213891</td><td>7.530890</td><td>6.449211</td><td>7.261644</td><td>7.395750</td><td>7.281587</td><td>7.620354</td><td>7.538805</td></tr>
	<tr><th scope=row>GSM892522</th><td>6.587242</td><td>7.161036</td><td>8.411995</td><td>7.162934</td><td>7.625582</td><td>6.343630</td><td>6.524414</td><td>8.697660</td><td>8.946331</td><td>7.628403</td><td>⋯</td><td>6.314725</td><td>7.075525</td><td> 9.849266</td><td>6.988772</td><td>6.730376</td><td>6.431582</td><td>6.985614</td><td>8.140622</td><td>8.675458</td><td>6.631154</td></tr>
	<tr><th scope=row>GSM892523</th><td>6.791856</td><td>8.253478</td><td>8.745131</td><td>7.927454</td><td>8.174660</td><td>6.963125</td><td>7.517947</td><td>8.762901</td><td>7.163435</td><td>8.739387</td><td>⋯</td><td>7.188488</td><td>7.622161</td><td> 9.148234</td><td>7.903907</td><td>6.306595</td><td>7.325780</td><td>7.334377</td><td>7.511323</td><td>7.650941</td><td>7.820282</td></tr>
	<tr><th scope=row>GSM892524</th><td>6.245855</td><td>7.671536</td><td>8.251558</td><td>7.391992</td><td>8.000434</td><td>6.286422</td><td>6.818797</td><td>7.967802</td><td>7.539410</td><td>7.715247</td><td>⋯</td><td>6.061974</td><td>7.027863</td><td> 8.259271</td><td>6.847237</td><td>6.849651</td><td>6.762690</td><td>6.710881</td><td>6.816398</td><td>8.360958</td><td>6.858813</td></tr>
	<tr><th scope=row>GSM892525</th><td>6.596868</td><td>8.229281</td><td>8.623323</td><td>7.784096</td><td>8.086842</td><td>6.771321</td><td>7.363387</td><td>8.629280</td><td>7.370437</td><td>8.418601</td><td>⋯</td><td>7.047557</td><td>7.618737</td><td> 9.329595</td><td>7.689948</td><td>6.369818</td><td>7.225717</td><td>7.376693</td><td>7.504508</td><td>7.956391</td><td>7.679804</td></tr>
	<tr><th scope=row>GSM892526</th><td>6.880399</td><td>7.963736</td><td>8.823275</td><td>7.667458</td><td>8.079904</td><td>7.004932</td><td>7.637911</td><td>8.700804</td><td>7.325166</td><td>8.490112</td><td>⋯</td><td>7.046551</td><td>7.717766</td><td> 9.615519</td><td>8.097353</td><td>6.651920</td><td>7.569444</td><td>7.341047</td><td>7.328069</td><td>7.704761</td><td>7.629483</td></tr>
	<tr><th scope=row>GSM892527</th><td>6.591960</td><td>8.003476</td><td>7.922650</td><td>7.644931</td><td>7.917782</td><td>6.790382</td><td>7.965053</td><td>8.175806</td><td>6.849212</td><td>8.580217</td><td>⋯</td><td>7.012370</td><td>7.580656</td><td> 9.392378</td><td>7.807005</td><td>6.274718</td><td>7.569269</td><td>7.284206</td><td>7.592794</td><td>7.356051</td><td>7.989371</td></tr>
	<tr><th scope=row>GSM892528</th><td>6.434061</td><td>8.431501</td><td>8.215267</td><td>7.856309</td><td>7.846105</td><td>6.780654</td><td>7.873020</td><td>8.416400</td><td>6.566155</td><td>8.516838</td><td>⋯</td><td>7.056252</td><td>7.580622</td><td> 9.186873</td><td>7.735456</td><td>6.483284</td><td>7.388772</td><td>7.445995</td><td>7.542868</td><td>8.729524</td><td>7.787358</td></tr>
	<tr><th scope=row>GSM892529</th><td>6.801279</td><td>8.228248</td><td>8.419544</td><td>8.283660</td><td>8.045290</td><td>6.386840</td><td>7.166886</td><td>7.812879</td><td>7.414223</td><td>9.056989</td><td>⋯</td><td>7.345189</td><td>7.454230</td><td> 9.698171</td><td>7.735707</td><td>6.920172</td><td>7.556130</td><td>7.316671</td><td>7.392102</td><td>8.336837</td><td>7.180302</td></tr>
	<tr><th scope=row>GSM892530</th><td>6.261183</td><td>7.501167</td><td>8.191854</td><td>7.966121</td><td>7.940615</td><td>6.552336</td><td>6.743088</td><td>8.778617</td><td>7.204549</td><td>8.153302</td><td>⋯</td><td>7.083876</td><td>7.224174</td><td> 8.806175</td><td>7.395989</td><td>6.474138</td><td>7.725389</td><td>7.216404</td><td>7.188727</td><td>9.718680</td><td>7.055392</td></tr>
	<tr><th scope=row>GSM892531</th><td>6.219873</td><td>7.582889</td><td>7.774610</td><td>7.567474</td><td>7.593643</td><td>6.403376</td><td>5.354727</td><td>7.809812</td><td>7.442859</td><td>8.586112</td><td>⋯</td><td>6.746093</td><td>6.685073</td><td> 9.000582</td><td>7.220321</td><td>6.792826</td><td>6.660123</td><td>6.984694</td><td>6.906403</td><td>9.461871</td><td>6.881501</td></tr>
	<tr><th scope=row>GSM892532</th><td>6.742950</td><td>8.320929</td><td>8.839143</td><td>7.924190</td><td>8.291375</td><td>6.569814</td><td>7.301869</td><td>8.796987</td><td>6.963654</td><td>8.385561</td><td>⋯</td><td>6.956621</td><td>7.521725</td><td> 9.533637</td><td>7.814331</td><td>6.430978</td><td>7.297665</td><td>7.430950</td><td>7.495326</td><td>7.720845</td><td>7.554511</td></tr>
	<tr><th scope=row>GSM892533</th><td>6.704944</td><td>7.412264</td><td>8.473529</td><td>7.956586</td><td>8.178099</td><td>6.089149</td><td>6.577089</td><td>8.721642</td><td>7.251669</td><td>7.476402</td><td>⋯</td><td>7.025184</td><td>6.514719</td><td> 8.166519</td><td>7.036061</td><td>6.757614</td><td>7.235078</td><td>6.754703</td><td>6.925736</td><td>9.670314</td><td>7.166199</td></tr>
	<tr><th scope=row>GSM892534</th><td>6.182084</td><td>7.144035</td><td>8.281567</td><td>7.870729</td><td>7.794427</td><td>6.127934</td><td>7.122769</td><td>8.834673</td><td>7.151682</td><td>7.486488</td><td>⋯</td><td>6.964867</td><td>7.488915</td><td> 9.595377</td><td>7.060084</td><td>6.836105</td><td>7.046702</td><td>6.761914</td><td>7.947575</td><td>9.083170</td><td>7.039949</td></tr>
	<tr><th scope=row>GSM892535</th><td>6.337192</td><td>8.366446</td><td>8.088356</td><td>7.671704</td><td>7.699893</td><td>6.773937</td><td>7.650213</td><td>8.284341</td><td>6.828798</td><td>8.681060</td><td>⋯</td><td>7.041633</td><td>7.485423</td><td> 9.461841</td><td>7.379815</td><td>6.425944</td><td>7.716786</td><td>7.287202</td><td>7.624202</td><td>7.607068</td><td>7.762243</td></tr>
	<tr><th scope=row>GSM892536</th><td>5.688162</td><td>6.933059</td><td>8.123178</td><td>7.463102</td><td>8.057786</td><td>6.333349</td><td>7.435997</td><td>8.847742</td><td>7.009881</td><td>6.966679</td><td>⋯</td><td>6.782929</td><td>6.560395</td><td> 9.354076</td><td>7.195035</td><td>6.671343</td><td>7.516456</td><td>6.767039</td><td>7.374235</td><td>9.384887</td><td>6.928326</td></tr>
	<tr><th scope=row>GSM892537</th><td>6.668096</td><td>8.540159</td><td>8.773367</td><td>7.833906</td><td>7.979246</td><td>6.788792</td><td>7.573157</td><td>8.660890</td><td>7.478607</td><td>8.410810</td><td>⋯</td><td>7.060252</td><td>7.422690</td><td> 9.650551</td><td>7.966881</td><td>6.454921</td><td>7.442702</td><td>7.277000</td><td>7.768285</td><td>7.903778</td><td>7.652281</td></tr>
	<tr><th scope=row>GSM892538</th><td>6.367278</td><td>7.982309</td><td>8.250587</td><td>7.747187</td><td>7.768878</td><td>6.531466</td><td>7.755633</td><td>8.219508</td><td>7.173153</td><td>8.561161</td><td>⋯</td><td>6.811171</td><td>7.615087</td><td> 9.919866</td><td>7.807143</td><td>6.402823</td><td>7.630336</td><td>7.411405</td><td>7.337271</td><td>7.679841</td><td>7.537877</td></tr>
	<tr><th scope=row>GSM892539</th><td>6.507528</td><td>7.889120</td><td>8.132082</td><td>7.773882</td><td>7.781864</td><td>6.812221</td><td>7.906423</td><td>8.213667</td><td>6.627808</td><td>8.710067</td><td>⋯</td><td>6.867453</td><td>7.584055</td><td> 8.958087</td><td>7.979214</td><td>6.504291</td><td>7.524379</td><td>7.231299</td><td>7.448077</td><td>7.491265</td><td>7.864270</td></tr>
	<tr><th scope=row>GSM892540</th><td>6.294595</td><td>7.948546</td><td>8.502375</td><td>7.807183</td><td>8.037376</td><td>6.929834</td><td>7.934246</td><td>8.358458</td><td>7.204202</td><td>8.306314</td><td>⋯</td><td>7.119266</td><td>7.703995</td><td> 9.813183</td><td>7.941755</td><td>6.279432</td><td>7.298962</td><td>7.355470</td><td>7.663083</td><td>7.579487</td><td>7.685918</td></tr>
	<tr><th scope=row>GSM892541</th><td>6.534798</td><td>8.048214</td><td>8.701014</td><td>7.776260</td><td>7.823811</td><td>6.490875</td><td>7.246580</td><td>8.262712</td><td>6.955950</td><td>8.693712</td><td>⋯</td><td>6.914460</td><td>7.524863</td><td>10.074511</td><td>7.741070</td><td>6.386677</td><td>7.474729</td><td>7.429227</td><td>7.379764</td><td>8.019258</td><td>7.365580</td></tr>
	<tr><th scope=row>GSM892542</th><td>6.184547</td><td>7.297091</td><td>8.421844</td><td>8.065336</td><td>8.315868</td><td>6.636063</td><td>6.882514</td><td>8.836912</td><td>7.158297</td><td>7.770093</td><td>⋯</td><td>6.968772</td><td>7.086351</td><td> 9.177625</td><td>7.350160</td><td>6.711161</td><td>7.376510</td><td>7.210573</td><td>7.689020</td><td>9.682822</td><td>7.192851</td></tr>
	<tr><th scope=row>GSM892543</th><td>6.682148</td><td>7.950826</td><td>8.690258</td><td>7.901859</td><td>8.012897</td><td>6.792236</td><td>7.206095</td><td>8.718302</td><td>7.237448</td><td>8.492662</td><td>⋯</td><td>6.900959</td><td>7.464481</td><td> 9.442640</td><td>7.874393</td><td>6.585328</td><td>7.254298</td><td>7.244149</td><td>7.140478</td><td>7.975937</td><td>7.541750</td></tr>
	<tr><th scope=row>GSM892544</th><td>6.733925</td><td>8.062635</td><td>8.722910</td><td>7.800925</td><td>8.062221</td><td>6.682119</td><td>7.104700</td><td>8.629489</td><td>7.041196</td><td>8.922461</td><td>⋯</td><td>7.214763</td><td>7.402581</td><td> 9.238843</td><td>7.562664</td><td>6.471253</td><td>7.149058</td><td>7.246975</td><td>7.499716</td><td>7.967075</td><td>7.733221</td></tr>
	<tr><th scope=row>GSM892545</th><td>6.696357</td><td>8.007845</td><td>8.620403</td><td>7.979154</td><td>8.280297</td><td>6.727042</td><td>7.497030</td><td>8.630631</td><td>7.153650</td><td>8.476046</td><td>⋯</td><td>6.992558</td><td>7.599650</td><td> 9.566730</td><td>7.820617</td><td>6.389318</td><td>7.240853</td><td>7.387185</td><td>7.413491</td><td>7.819232</td><td>7.543667</td></tr>
	<tr><th scope=row>GSM892546</th><td>6.601588</td><td>7.878466</td><td>8.451771</td><td>7.831509</td><td>7.983359</td><td>6.661006</td><td>7.841769</td><td>8.233061</td><td>6.618159</td><td>8.413212</td><td>⋯</td><td>6.855004</td><td>7.468985</td><td> 9.409230</td><td>7.853540</td><td>6.398046</td><td>7.021363</td><td>7.325438</td><td>7.700557</td><td>7.671356</td><td>7.845263</td></tr>
	<tr><th scope=row>GSM892547</th><td>6.651946</td><td>7.434621</td><td>8.308435</td><td>7.983728</td><td>8.212196</td><td>5.414503</td><td>6.959527</td><td>8.632683</td><td>6.911231</td><td>7.349661</td><td>⋯</td><td>6.765987</td><td>7.413592</td><td> 9.472469</td><td>6.938442</td><td>7.085260</td><td>6.797690</td><td>6.738146</td><td>7.718725</td><td>9.098681</td><td>7.301284</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>GSM892711</th><td>6.476200</td><td>8.363302</td><td>8.127312</td><td>6.978710</td><td>8.151245</td><td>6.465365</td><td>5.171469</td><td>7.381079</td><td>6.208254</td><td>5.508167</td><td>⋯</td><td>6.826936</td><td>5.407265</td><td>7.618136</td><td>6.846758</td><td>7.431044</td><td>4.863038</td><td>6.485421</td><td>6.043969</td><td> 8.580735</td><td>6.389012</td></tr>
	<tr><th scope=row>GSM892712</th><td>6.289872</td><td>7.490396</td><td>7.613947</td><td>7.325866</td><td>8.625121</td><td>6.873631</td><td>6.088767</td><td>7.830076</td><td>5.972325</td><td>7.047833</td><td>⋯</td><td>7.083781</td><td>6.897591</td><td>8.618126</td><td>6.575474</td><td>6.417783</td><td>6.575952</td><td>7.227973</td><td>7.687075</td><td> 9.358290</td><td>7.200833</td></tr>
	<tr><th scope=row>GSM892713</th><td>6.568658</td><td>7.914156</td><td>7.760917</td><td>7.092881</td><td>8.585764</td><td>6.907888</td><td>6.837888</td><td>7.927918</td><td>6.287703</td><td>6.591545</td><td>⋯</td><td>7.070480</td><td>7.000971</td><td>9.271294</td><td>6.947144</td><td>6.701812</td><td>6.288006</td><td>5.968029</td><td>7.576178</td><td> 8.985936</td><td>6.940660</td></tr>
	<tr><th scope=row>GSM892714</th><td>6.260210</td><td>7.432401</td><td>7.981251</td><td>7.065626</td><td>8.443269</td><td>6.734447</td><td>6.098279</td><td>7.850232</td><td>6.082650</td><td>6.865432</td><td>⋯</td><td>7.381078</td><td>6.824001</td><td>8.539393</td><td>6.414655</td><td>6.943558</td><td>6.156067</td><td>7.020572</td><td>8.010377</td><td> 8.809008</td><td>7.363817</td></tr>
	<tr><th scope=row>GSM892715</th><td>6.634492</td><td>6.957463</td><td>7.636013</td><td>7.313244</td><td>8.402642</td><td>6.522632</td><td>5.733629</td><td>7.960610</td><td>6.249646</td><td>6.506355</td><td>⋯</td><td>6.743300</td><td>6.730268</td><td>8.261791</td><td>6.178705</td><td>6.727753</td><td>5.388038</td><td>6.122999</td><td>7.332406</td><td> 9.623936</td><td>6.773690</td></tr>
	<tr><th scope=row>GSM892716</th><td>6.786530</td><td>7.646512</td><td>7.666652</td><td>7.172274</td><td>8.670431</td><td>6.564646</td><td>6.254124</td><td>7.691596</td><td>6.029386</td><td>6.677923</td><td>⋯</td><td>6.118135</td><td>6.929292</td><td>9.247521</td><td>6.224404</td><td>6.811507</td><td>6.310932</td><td>6.539631</td><td>7.927077</td><td> 9.073853</td><td>7.028067</td></tr>
	<tr><th scope=row>GSM892717</th><td>6.592317</td><td>7.368495</td><td>7.838677</td><td>7.605901</td><td>8.409795</td><td>6.940597</td><td>5.609481</td><td>7.956074</td><td>6.353457</td><td>7.303067</td><td>⋯</td><td>6.981124</td><td>5.722455</td><td>8.106521</td><td>6.724852</td><td>6.402136</td><td>6.617191</td><td>6.825132</td><td>6.513575</td><td> 9.657655</td><td>7.325079</td></tr>
	<tr><th scope=row>GSM892718</th><td>6.768997</td><td>7.508028</td><td>7.630184</td><td>6.978505</td><td>8.792552</td><td>6.415695</td><td>6.269389</td><td>7.983391</td><td>6.211742</td><td>6.891727</td><td>⋯</td><td>7.065066</td><td>7.068208</td><td>9.388699</td><td>6.831353</td><td>6.813834</td><td>6.423029</td><td>6.690049</td><td>7.769320</td><td> 9.251383</td><td>6.997929</td></tr>
	<tr><th scope=row>GSM892719</th><td>6.137274</td><td>8.175573</td><td>7.801282</td><td>7.166578</td><td>8.116656</td><td>6.208456</td><td>6.940733</td><td>7.695558</td><td>5.893956</td><td>7.040804</td><td>⋯</td><td>6.530052</td><td>6.432623</td><td>9.373218</td><td>6.391285</td><td>6.543695</td><td>6.445889</td><td>6.984509</td><td>7.989008</td><td> 9.096490</td><td>6.993594</td></tr>
	<tr><th scope=row>GSM892720</th><td>6.161883</td><td>7.579085</td><td>7.689805</td><td>7.206432</td><td>8.669464</td><td>6.052609</td><td>5.501467</td><td>7.721230</td><td>6.794644</td><td>5.944503</td><td>⋯</td><td>6.317202</td><td>6.165129</td><td>8.786184</td><td>6.454978</td><td>6.557341</td><td>5.479841</td><td>6.484029</td><td>6.937832</td><td> 9.177315</td><td>6.756525</td></tr>
	<tr><th scope=row>GSM892721</th><td>6.936574</td><td>7.887124</td><td>7.895321</td><td>7.279579</td><td>8.645885</td><td>6.510732</td><td>6.423610</td><td>7.879193</td><td>5.922117</td><td>7.442526</td><td>⋯</td><td>7.009328</td><td>6.659047</td><td>8.699425</td><td>6.770520</td><td>6.471355</td><td>6.326404</td><td>6.797455</td><td>7.409474</td><td> 9.768182</td><td>7.322568</td></tr>
	<tr><th scope=row>GSM892722</th><td>6.971671</td><td>7.809692</td><td>7.915753</td><td>7.543959</td><td>8.447548</td><td>7.173433</td><td>7.029512</td><td>7.729862</td><td>5.969891</td><td>7.465948</td><td>⋯</td><td>7.354816</td><td>6.945015</td><td>8.739375</td><td>6.976216</td><td>6.543381</td><td>6.567784</td><td>6.927991</td><td>7.636797</td><td> 7.478889</td><td>7.593385</td></tr>
	<tr><th scope=row>GSM892723</th><td>5.985059</td><td>7.757325</td><td>7.717067</td><td>6.649571</td><td>7.908695</td><td>6.427296</td><td>5.216108</td><td>7.698404</td><td>6.671056</td><td>5.161368</td><td>⋯</td><td>4.862884</td><td>6.089914</td><td>9.171755</td><td>5.619530</td><td>7.008153</td><td>4.875154</td><td>5.629678</td><td>7.004463</td><td> 9.087836</td><td>5.865254</td></tr>
	<tr><th scope=row>GSM892724</th><td>6.506324</td><td>7.006040</td><td>7.831571</td><td>7.516133</td><td>8.455084</td><td>6.612935</td><td>5.514891</td><td>7.436909</td><td>6.393400</td><td>7.150989</td><td>⋯</td><td>6.910452</td><td>6.560502</td><td>8.281061</td><td>6.025772</td><td>6.736728</td><td>5.294410</td><td>6.641917</td><td>7.198413</td><td>10.231938</td><td>6.875560</td></tr>
	<tr><th scope=row>GSM892725</th><td>6.652975</td><td>7.659634</td><td>7.798400</td><td>7.000333</td><td>8.426011</td><td>6.391030</td><td>5.312481</td><td>7.291015</td><td>6.407956</td><td>6.162020</td><td>⋯</td><td>6.546185</td><td>5.910178</td><td>8.350952</td><td>6.298148</td><td>6.620398</td><td>5.696531</td><td>6.408603</td><td>6.890606</td><td> 8.721910</td><td>6.641521</td></tr>
	<tr><th scope=row>GSM892726</th><td>6.461473</td><td>7.817090</td><td>7.769957</td><td>7.375247</td><td>8.326782</td><td>6.366425</td><td>5.738586</td><td>8.004114</td><td>6.113817</td><td>7.552871</td><td>⋯</td><td>7.066131</td><td>6.865552</td><td>9.047524</td><td>6.708812</td><td>6.383058</td><td>6.532482</td><td>6.840590</td><td>7.625847</td><td> 9.777910</td><td>7.469105</td></tr>
	<tr><th scope=row>GSM892727</th><td>6.792632</td><td>7.382763</td><td>7.766917</td><td>7.208697</td><td>8.617717</td><td>6.863750</td><td>6.445228</td><td>7.949494</td><td>5.516677</td><td>6.813826</td><td>⋯</td><td>7.369951</td><td>7.008676</td><td>9.045257</td><td>6.182111</td><td>6.959245</td><td>6.167560</td><td>6.555311</td><td>7.949245</td><td> 8.958348</td><td>6.625944</td></tr>
	<tr><th scope=row>GSM892728</th><td>6.489961</td><td>8.013562</td><td>7.666948</td><td>7.024002</td><td>8.189994</td><td>6.037001</td><td>6.772912</td><td>7.672432</td><td>5.749248</td><td>7.261598</td><td>⋯</td><td>6.676088</td><td>6.701813</td><td>9.295983</td><td>6.268745</td><td>6.889879</td><td>6.199361</td><td>6.614349</td><td>8.226501</td><td> 8.897003</td><td>6.940988</td></tr>
	<tr><th scope=row>GSM892729</th><td>6.545018</td><td>8.164848</td><td>7.761506</td><td>7.121277</td><td>8.182803</td><td>6.789162</td><td>6.954836</td><td>7.498273</td><td>5.957108</td><td>7.735294</td><td>⋯</td><td>6.837937</td><td>6.855791</td><td>9.321107</td><td>7.010450</td><td>6.159609</td><td>6.830631</td><td>7.143977</td><td>7.602906</td><td> 7.395089</td><td>7.589764</td></tr>
	<tr><th scope=row>GSM892730</th><td>6.044620</td><td>8.188033</td><td>7.755512</td><td>6.982878</td><td>8.267727</td><td>6.048832</td><td>6.574394</td><td>7.617016</td><td>6.013510</td><td>7.151277</td><td>⋯</td><td>6.566443</td><td>6.537270</td><td>9.455387</td><td>6.349760</td><td>6.613317</td><td>6.182800</td><td>6.735401</td><td>8.260799</td><td> 9.235865</td><td>6.998553</td></tr>
	<tr><th scope=row>GSM892731</th><td>6.903474</td><td>8.426809</td><td>7.934202</td><td>7.236329</td><td>8.252199</td><td>6.443834</td><td>7.058370</td><td>7.553470</td><td>5.928988</td><td>7.647825</td><td>⋯</td><td>7.023655</td><td>6.625362</td><td>9.343620</td><td>6.810409</td><td>6.288356</td><td>6.420266</td><td>7.184801</td><td>7.625527</td><td> 7.338900</td><td>7.614750</td></tr>
	<tr><th scope=row>GSM892732</th><td>6.232739</td><td>8.138495</td><td>7.851660</td><td>7.349124</td><td>8.364646</td><td>6.412122</td><td>5.993294</td><td>8.037159</td><td>6.694239</td><td>7.992721</td><td>⋯</td><td>6.454157</td><td>6.340876</td><td>8.846103</td><td>6.864891</td><td>6.397387</td><td>6.445142</td><td>6.852658</td><td>7.347422</td><td>10.269789</td><td>7.473535</td></tr>
	<tr><th scope=row>GSM892733</th><td>6.655518</td><td>7.966378</td><td>7.698371</td><td>7.400564</td><td>8.435264</td><td>6.497104</td><td>6.140997</td><td>8.035092</td><td>5.929801</td><td>7.638686</td><td>⋯</td><td>7.028635</td><td>6.681626</td><td>8.824372</td><td>6.899025</td><td>6.420495</td><td>6.717390</td><td>6.863758</td><td>7.543075</td><td>10.238442</td><td>7.583005</td></tr>
	<tr><th scope=row>GSM892734</th><td>6.699976</td><td>8.217901</td><td>6.900384</td><td>7.340885</td><td>8.357698</td><td>5.696248</td><td>5.470826</td><td>7.486252</td><td>5.895026</td><td>7.442163</td><td>⋯</td><td>6.995989</td><td>6.367713</td><td>9.026293</td><td>6.163627</td><td>6.833714</td><td>6.542262</td><td>6.907784</td><td>7.632727</td><td> 9.983087</td><td>6.828129</td></tr>
	<tr><th scope=row>GSM892735</th><td>6.635070</td><td>7.519029</td><td>7.757553</td><td>7.251901</td><td>8.443706</td><td>6.477393</td><td>5.618210</td><td>7.735360</td><td>6.281759</td><td>7.398695</td><td>⋯</td><td>6.645776</td><td>6.326635</td><td>9.226130</td><td>6.521542</td><td>6.420993</td><td>6.214451</td><td>6.824885</td><td>7.371379</td><td> 9.940597</td><td>6.977336</td></tr>
	<tr><th scope=row>GSM892736</th><td>6.688627</td><td>7.559205</td><td>7.848844</td><td>7.423926</td><td>8.409471</td><td>6.352663</td><td>5.792368</td><td>7.918395</td><td>5.788065</td><td>7.345189</td><td>⋯</td><td>6.885424</td><td>6.600315</td><td>8.845849</td><td>6.349378</td><td>6.653121</td><td>6.177181</td><td>6.458729</td><td>7.380818</td><td> 9.878802</td><td>6.996492</td></tr>
	<tr><th scope=row>GSM892737</th><td>6.461210</td><td>8.076686</td><td>8.114321</td><td>6.987587</td><td>8.770922</td><td>6.336274</td><td>6.229084</td><td>7.698160</td><td>5.769223</td><td>7.038824</td><td>⋯</td><td>7.052092</td><td>6.336557</td><td>9.200771</td><td>6.769130</td><td>6.694464</td><td>6.213077</td><td>6.658322</td><td>8.189443</td><td> 9.201405</td><td>7.541182</td></tr>
	<tr><th scope=row>GSM892738</th><td>6.840229</td><td>7.750591</td><td>7.479951</td><td>7.036618</td><td>8.441403</td><td>6.394248</td><td>6.980082</td><td>7.760172</td><td>5.634933</td><td>7.293573</td><td>⋯</td><td>7.277475</td><td>6.773417</td><td>9.119492</td><td>6.533849</td><td>6.587136</td><td>6.510829</td><td>6.877735</td><td>8.362280</td><td> 8.578688</td><td>7.402224</td></tr>
	<tr><th scope=row>GSM892739</th><td>6.334797</td><td>8.555989</td><td>7.770727</td><td>7.252122</td><td>8.212847</td><td>6.438389</td><td>7.001591</td><td>7.465521</td><td>5.946766</td><td>7.722824</td><td>⋯</td><td>6.981306</td><td>6.687297</td><td>9.478243</td><td>7.030569</td><td>6.234560</td><td>6.907458</td><td>7.112260</td><td>7.348922</td><td> 7.519751</td><td>7.430321</td></tr>
	<tr><th scope=row>GSM892740</th><td>6.122527</td><td>7.977229</td><td>7.734951</td><td>7.289281</td><td>8.308845</td><td>6.052669</td><td>5.713552</td><td>7.948864</td><td>6.569940</td><td>7.291691</td><td>⋯</td><td>6.709204</td><td>6.423148</td><td>9.210131</td><td>6.427857</td><td>6.663105</td><td>6.144444</td><td>6.756438</td><td>7.566873</td><td>10.409470</td><td>7.066382</td></tr>
</tbody>
</table>




```R
gene_time  <- normalized_gset_mean_high_expression %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')  %>% 
      as_tibble() 
gene_time
```


<table class="dataframe">
<caption>A tibble: 223 × 5964</caption>
<thead>
	<tr><th scope=col>Samples</th><th scope=col>AACS</th><th scope=col>AAED1</th><th scope=col>AAMP</th><th scope=col>AAR2</th><th scope=col>AARS</th><th scope=col>AARSD1 /// PTGES3L /// PTGES3L-AARSD1</th><th scope=col>AASDH</th><th scope=col>AATF</th><th scope=col>AATK</th><th scope=col>⋯</th><th scope=col>ZUFSP</th><th scope=col>ZW10</th><th scope=col>ZYG11B</th><th scope=col>ZYX</th><th scope=col>ZZZ3</th><th scope=col>Ages</th><th scope=col>Type</th><th scope=col>Patients</th><th scope=col>Outcome</th><th scope=col>Time</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM892518</td><td>6.744481</td><td>8.098326</td><td>8.845298</td><td>8.027744</td><td>8.133400</td><td>6.779367</td><td>7.192553</td><td>8.788231</td><td>7.296900</td><td>⋯</td><td>7.424264</td><td>7.443882</td><td>7.309862</td><td>7.839783</td><td>7.745338</td><td>68</td><td>PBMC</td><td>279</td><td>FALSE</td><td>1738</td></tr>
	<tr><td>GSM892519</td><td>6.606717</td><td>7.951058</td><td>8.353823</td><td>7.704787</td><td>7.991139</td><td>6.713854</td><td>7.712399</td><td>8.498978</td><td>6.748701</td><td>⋯</td><td>7.568137</td><td>7.431788</td><td>7.383383</td><td>7.792338</td><td>7.728572</td><td>78</td><td>PBMC</td><td>281</td><td>FALSE</td><td>697 </td></tr>
	<tr><td>GSM892520</td><td>6.703049</td><td>8.253258</td><td>8.821917</td><td>7.935994</td><td>8.121270</td><td>6.987739</td><td>7.483351</td><td>8.788333</td><td>7.220110</td><td>⋯</td><td>7.342488</td><td>7.437652</td><td>7.422046</td><td>7.996883</td><td>7.661320</td><td>84</td><td>PBMC</td><td>282</td><td>TRUE </td><td>2   </td></tr>
	<tr><td>GSM892521</td><td>6.817530</td><td>8.024947</td><td>8.754522</td><td>7.956263</td><td>8.256909</td><td>6.871659</td><td>7.263732</td><td>8.747234</td><td>7.395759</td><td>⋯</td><td>7.261644</td><td>7.395750</td><td>7.281587</td><td>7.620354</td><td>7.538805</td><td>80</td><td>PBMC</td><td>284</td><td>FALSE</td><td>1767</td></tr>
	<tr><td>GSM892522</td><td>6.587242</td><td>7.161036</td><td>8.411995</td><td>7.162934</td><td>7.625582</td><td>6.343630</td><td>6.524414</td><td>8.697660</td><td>8.946331</td><td>⋯</td><td>6.431582</td><td>6.985614</td><td>8.140622</td><td>8.675458</td><td>6.631154</td><td>81</td><td>PBMC</td><td>285</td><td>FALSE</td><td>1326</td></tr>
	<tr><td>GSM892523</td><td>6.791856</td><td>8.253478</td><td>8.745131</td><td>7.927454</td><td>8.174660</td><td>6.963125</td><td>7.517947</td><td>8.762901</td><td>7.163435</td><td>⋯</td><td>7.325780</td><td>7.334377</td><td>7.511323</td><td>7.650941</td><td>7.820282</td><td>50</td><td>PBMC</td><td>286</td><td>FALSE</td><td>1427</td></tr>
	<tr><td>GSM892524</td><td>6.245855</td><td>7.671536</td><td>8.251558</td><td>7.391992</td><td>8.000434</td><td>6.286422</td><td>6.818797</td><td>7.967802</td><td>7.539410</td><td>⋯</td><td>6.762690</td><td>6.710881</td><td>6.816398</td><td>8.360958</td><td>6.858813</td><td>75</td><td>PBMC</td><td>288</td><td>TRUE </td><td>7   </td></tr>
	<tr><td>GSM892525</td><td>6.596868</td><td>8.229281</td><td>8.623323</td><td>7.784096</td><td>8.086842</td><td>6.771321</td><td>7.363387</td><td>8.629280</td><td>7.370437</td><td>⋯</td><td>7.225717</td><td>7.376693</td><td>7.504508</td><td>7.956391</td><td>7.679804</td><td>68</td><td>PBMC</td><td>292</td><td>FALSE</td><td>1422</td></tr>
	<tr><td>GSM892526</td><td>6.880399</td><td>7.963736</td><td>8.823275</td><td>7.667458</td><td>8.079904</td><td>7.004932</td><td>7.637911</td><td>8.700804</td><td>7.325166</td><td>⋯</td><td>7.569444</td><td>7.341047</td><td>7.328069</td><td>7.704761</td><td>7.629483</td><td>73</td><td>PBMC</td><td>296</td><td>FALSE</td><td>1320</td></tr>
	<tr><td>GSM892527</td><td>6.591960</td><td>8.003476</td><td>7.922650</td><td>7.644931</td><td>7.917782</td><td>6.790382</td><td>7.965053</td><td>8.175806</td><td>6.849212</td><td>⋯</td><td>7.569269</td><td>7.284206</td><td>7.592794</td><td>7.356051</td><td>7.989371</td><td>75</td><td>PBMC</td><td>298</td><td>FALSE</td><td>689 </td></tr>
	<tr><td>GSM892528</td><td>6.434061</td><td>8.431501</td><td>8.215267</td><td>7.856309</td><td>7.846105</td><td>6.780654</td><td>7.873020</td><td>8.416400</td><td>6.566155</td><td>⋯</td><td>7.388772</td><td>7.445995</td><td>7.542868</td><td>8.729524</td><td>7.787358</td><td>64</td><td>PBMC</td><td>301</td><td>FALSE</td><td>676 </td></tr>
	<tr><td>GSM892529</td><td>6.801279</td><td>8.228248</td><td>8.419544</td><td>8.283660</td><td>8.045290</td><td>6.386840</td><td>7.166886</td><td>7.812879</td><td>7.414223</td><td>⋯</td><td>7.556130</td><td>7.316671</td><td>7.392102</td><td>8.336837</td><td>7.180302</td><td>66</td><td>PBMC</td><td>302</td><td>FALSE</td><td>1710</td></tr>
	<tr><td>GSM892530</td><td>6.261183</td><td>7.501167</td><td>8.191854</td><td>7.966121</td><td>7.940615</td><td>6.552336</td><td>6.743088</td><td>8.778617</td><td>7.204549</td><td>⋯</td><td>7.725389</td><td>7.216404</td><td>7.188727</td><td>9.718680</td><td>7.055392</td><td>71</td><td>PBMC</td><td>306</td><td>FALSE</td><td>1688</td></tr>
	<tr><td>GSM892531</td><td>6.219873</td><td>7.582889</td><td>7.774610</td><td>7.567474</td><td>7.593643</td><td>6.403376</td><td>5.354727</td><td>7.809812</td><td>7.442859</td><td>⋯</td><td>6.660123</td><td>6.984694</td><td>6.906403</td><td>9.461871</td><td>6.881501</td><td>56</td><td>PBMC</td><td>311</td><td>FALSE</td><td>1711</td></tr>
	<tr><td>GSM892532</td><td>6.742950</td><td>8.320929</td><td>8.839143</td><td>7.924190</td><td>8.291375</td><td>6.569814</td><td>7.301869</td><td>8.796987</td><td>6.963654</td><td>⋯</td><td>7.297665</td><td>7.430950</td><td>7.495326</td><td>7.720845</td><td>7.554511</td><td>76</td><td>PBMC</td><td>313</td><td>FALSE</td><td>751 </td></tr>
	<tr><td>GSM892533</td><td>6.704944</td><td>7.412264</td><td>8.473529</td><td>7.956586</td><td>8.178099</td><td>6.089149</td><td>6.577089</td><td>8.721642</td><td>7.251669</td><td>⋯</td><td>7.235078</td><td>6.754703</td><td>6.925736</td><td>9.670314</td><td>7.166199</td><td>78</td><td>PBMC</td><td>318</td><td>FALSE</td><td>1680</td></tr>
	<tr><td>GSM892534</td><td>6.182084</td><td>7.144035</td><td>8.281567</td><td>7.870729</td><td>7.794427</td><td>6.127934</td><td>7.122769</td><td>8.834673</td><td>7.151682</td><td>⋯</td><td>7.046702</td><td>6.761914</td><td>7.947575</td><td>9.083170</td><td>7.039949</td><td>69</td><td>PBMC</td><td>327</td><td>FALSE</td><td>1876</td></tr>
	<tr><td>GSM892535</td><td>6.337192</td><td>8.366446</td><td>8.088356</td><td>7.671704</td><td>7.699893</td><td>6.773937</td><td>7.650213</td><td>8.284341</td><td>6.828798</td><td>⋯</td><td>7.716786</td><td>7.287202</td><td>7.624202</td><td>7.607068</td><td>7.762243</td><td>78</td><td>PBMC</td><td>333</td><td>FALSE</td><td>627 </td></tr>
	<tr><td>GSM892536</td><td>5.688162</td><td>6.933059</td><td>8.123178</td><td>7.463102</td><td>8.057786</td><td>6.333349</td><td>7.435997</td><td>8.847742</td><td>7.009881</td><td>⋯</td><td>7.516456</td><td>6.767039</td><td>7.374235</td><td>9.384887</td><td>6.928326</td><td>66</td><td>PBMC</td><td>334</td><td>FALSE</td><td>1779</td></tr>
	<tr><td>GSM892537</td><td>6.668096</td><td>8.540159</td><td>8.773367</td><td>7.833906</td><td>7.979246</td><td>6.788792</td><td>7.573157</td><td>8.660890</td><td>7.478607</td><td>⋯</td><td>7.442702</td><td>7.277000</td><td>7.768285</td><td>7.903778</td><td>7.652281</td><td>84</td><td>PBMC</td><td>337</td><td>FALSE</td><td>1542</td></tr>
	<tr><td>GSM892538</td><td>6.367278</td><td>7.982309</td><td>8.250587</td><td>7.747187</td><td>7.768878</td><td>6.531466</td><td>7.755633</td><td>8.219508</td><td>7.173153</td><td>⋯</td><td>7.630336</td><td>7.411405</td><td>7.337271</td><td>7.679841</td><td>7.537877</td><td>79</td><td>PBMC</td><td>340</td><td>FALSE</td><td>777 </td></tr>
	<tr><td>GSM892539</td><td>6.507528</td><td>7.889120</td><td>8.132082</td><td>7.773882</td><td>7.781864</td><td>6.812221</td><td>7.906423</td><td>8.213667</td><td>6.627808</td><td>⋯</td><td>7.524379</td><td>7.231299</td><td>7.448077</td><td>7.491265</td><td>7.864270</td><td>78</td><td>PBMC</td><td>345</td><td>FALSE</td><td>687 </td></tr>
	<tr><td>GSM892540</td><td>6.294595</td><td>7.948546</td><td>8.502375</td><td>7.807183</td><td>8.037376</td><td>6.929834</td><td>7.934246</td><td>8.358458</td><td>7.204202</td><td>⋯</td><td>7.298962</td><td>7.355470</td><td>7.663083</td><td>7.579487</td><td>7.685918</td><td>65</td><td>PBMC</td><td>355</td><td>FALSE</td><td>739 </td></tr>
	<tr><td>GSM892541</td><td>6.534798</td><td>8.048214</td><td>8.701014</td><td>7.776260</td><td>7.823811</td><td>6.490875</td><td>7.246580</td><td>8.262712</td><td>6.955950</td><td>⋯</td><td>7.474729</td><td>7.429227</td><td>7.379764</td><td>8.019258</td><td>7.365580</td><td>77</td><td>PBMC</td><td>356</td><td>TRUE </td><td>2   </td></tr>
	<tr><td>GSM892542</td><td>6.184547</td><td>7.297091</td><td>8.421844</td><td>8.065336</td><td>8.315868</td><td>6.636063</td><td>6.882514</td><td>8.836912</td><td>7.158297</td><td>⋯</td><td>7.376510</td><td>7.210573</td><td>7.689020</td><td>9.682822</td><td>7.192851</td><td>79</td><td>PBMC</td><td>359</td><td>FALSE</td><td>1674</td></tr>
	<tr><td>GSM892543</td><td>6.682148</td><td>7.950826</td><td>8.690258</td><td>7.901859</td><td>8.012897</td><td>6.792236</td><td>7.206095</td><td>8.718302</td><td>7.237448</td><td>⋯</td><td>7.254298</td><td>7.244149</td><td>7.140478</td><td>7.975937</td><td>7.541750</td><td>72</td><td>PBMC</td><td>363</td><td>FALSE</td><td>1649</td></tr>
	<tr><td>GSM892544</td><td>6.733925</td><td>8.062635</td><td>8.722910</td><td>7.800925</td><td>8.062221</td><td>6.682119</td><td>7.104700</td><td>8.629489</td><td>7.041196</td><td>⋯</td><td>7.149058</td><td>7.246975</td><td>7.499716</td><td>7.967075</td><td>7.733221</td><td>63</td><td>PBMC</td><td>364</td><td>FALSE</td><td>1842</td></tr>
	<tr><td>GSM892545</td><td>6.696357</td><td>8.007845</td><td>8.620403</td><td>7.979154</td><td>8.280297</td><td>6.727042</td><td>7.497030</td><td>8.630631</td><td>7.153650</td><td>⋯</td><td>7.240853</td><td>7.387185</td><td>7.413491</td><td>7.819232</td><td>7.543667</td><td>81</td><td>PBMC</td><td>365</td><td>TRUE </td><td>922 </td></tr>
	<tr><td>GSM892546</td><td>6.601588</td><td>7.878466</td><td>8.451771</td><td>7.831509</td><td>7.983359</td><td>6.661006</td><td>7.841769</td><td>8.233061</td><td>6.618159</td><td>⋯</td><td>7.021363</td><td>7.325438</td><td>7.700557</td><td>7.671356</td><td>7.845263</td><td>59</td><td>PBMC</td><td>366</td><td>TRUE </td><td>1   </td></tr>
	<tr><td>GSM892547</td><td>6.651946</td><td>7.434621</td><td>8.308435</td><td>7.983728</td><td>8.212196</td><td>5.414503</td><td>6.959527</td><td>8.632683</td><td>6.911231</td><td>⋯</td><td>6.797690</td><td>6.738146</td><td>7.718725</td><td>9.098681</td><td>7.301284</td><td>69</td><td>PBMC</td><td>367</td><td>FALSE</td><td>1886</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>GSM892711</td><td>6.476200</td><td>8.363302</td><td>8.127312</td><td>6.978710</td><td>8.151245</td><td>6.465365</td><td>5.171469</td><td>7.381079</td><td>6.208254</td><td>⋯</td><td>4.863038</td><td>6.485421</td><td>6.043969</td><td> 8.580735</td><td>6.389012</td><td>61</td><td>Plaques</td><td>587</td><td>FALSE</td><td>2083</td></tr>
	<tr><td>GSM892712</td><td>6.289872</td><td>7.490396</td><td>7.613947</td><td>7.325866</td><td>8.625121</td><td>6.873631</td><td>6.088767</td><td>7.830076</td><td>5.972325</td><td>⋯</td><td>6.575952</td><td>7.227973</td><td>7.687075</td><td> 9.358290</td><td>7.200833</td><td>48</td><td>Plaques</td><td>588</td><td>FALSE</td><td>1704</td></tr>
	<tr><td>GSM892713</td><td>6.568658</td><td>7.914156</td><td>7.760917</td><td>7.092881</td><td>8.585764</td><td>6.907888</td><td>6.837888</td><td>7.927918</td><td>6.287703</td><td>⋯</td><td>6.288006</td><td>5.968029</td><td>7.576178</td><td> 8.985936</td><td>6.940660</td><td>77</td><td>Plaques</td><td>591</td><td>TRUE </td><td>21  </td></tr>
	<tr><td>GSM892714</td><td>6.260210</td><td>7.432401</td><td>7.981251</td><td>7.065626</td><td>8.443269</td><td>6.734447</td><td>6.098279</td><td>7.850232</td><td>6.082650</td><td>⋯</td><td>6.156067</td><td>7.020572</td><td>8.010377</td><td> 8.809008</td><td>7.363817</td><td>66</td><td>Plaques</td><td>597</td><td>FALSE</td><td>2115</td></tr>
	<tr><td>GSM892715</td><td>6.634492</td><td>6.957463</td><td>7.636013</td><td>7.313244</td><td>8.402642</td><td>6.522632</td><td>5.733629</td><td>7.960610</td><td>6.249646</td><td>⋯</td><td>5.388038</td><td>6.122999</td><td>7.332406</td><td> 9.623936</td><td>6.773690</td><td>73</td><td>Plaques</td><td>603</td><td>FALSE</td><td>1306</td></tr>
	<tr><td>GSM892716</td><td>6.786530</td><td>7.646512</td><td>7.666652</td><td>7.172274</td><td>8.670431</td><td>6.564646</td><td>6.254124</td><td>7.691596</td><td>6.029386</td><td>⋯</td><td>6.310932</td><td>6.539631</td><td>7.927077</td><td> 9.073853</td><td>7.028067</td><td>77</td><td>Plaques</td><td>605</td><td>FALSE</td><td>1870</td></tr>
	<tr><td>GSM892717</td><td>6.592317</td><td>7.368495</td><td>7.838677</td><td>7.605901</td><td>8.409795</td><td>6.940597</td><td>5.609481</td><td>7.956074</td><td>6.353457</td><td>⋯</td><td>6.617191</td><td>6.825132</td><td>6.513575</td><td> 9.657655</td><td>7.325079</td><td>62</td><td>Plaques</td><td>606</td><td>FALSE</td><td>1173</td></tr>
	<tr><td>GSM892718</td><td>6.768997</td><td>7.508028</td><td>7.630184</td><td>6.978505</td><td>8.792552</td><td>6.415695</td><td>6.269389</td><td>7.983391</td><td>6.211742</td><td>⋯</td><td>6.423029</td><td>6.690049</td><td>7.769320</td><td> 9.251383</td><td>6.997929</td><td>82</td><td>Plaques</td><td>618</td><td>FALSE</td><td>395 </td></tr>
	<tr><td>GSM892719</td><td>6.137274</td><td>8.175573</td><td>7.801282</td><td>7.166578</td><td>8.116656</td><td>6.208456</td><td>6.940733</td><td>7.695558</td><td>5.893956</td><td>⋯</td><td>6.445889</td><td>6.984509</td><td>7.989008</td><td> 9.096490</td><td>6.993594</td><td>80</td><td>Plaques</td><td>619</td><td>FALSE</td><td>2619</td></tr>
	<tr><td>GSM892720</td><td>6.161883</td><td>7.579085</td><td>7.689805</td><td>7.206432</td><td>8.669464</td><td>6.052609</td><td>5.501467</td><td>7.721230</td><td>6.794644</td><td>⋯</td><td>5.479841</td><td>6.484029</td><td>6.937832</td><td> 9.177315</td><td>6.756525</td><td>68</td><td>Plaques</td><td>624</td><td>FALSE</td><td>1752</td></tr>
	<tr><td>GSM892721</td><td>6.936574</td><td>7.887124</td><td>7.895321</td><td>7.279579</td><td>8.645885</td><td>6.510732</td><td>6.423610</td><td>7.879193</td><td>5.922117</td><td>⋯</td><td>6.326404</td><td>6.797455</td><td>7.409474</td><td> 9.768182</td><td>7.322568</td><td>70</td><td>Plaques</td><td>625</td><td>FALSE</td><td>1207</td></tr>
	<tr><td>GSM892722</td><td>6.971671</td><td>7.809692</td><td>7.915753</td><td>7.543959</td><td>8.447548</td><td>7.173433</td><td>7.029512</td><td>7.729862</td><td>5.969891</td><td>⋯</td><td>6.567784</td><td>6.927991</td><td>7.636797</td><td> 7.478889</td><td>7.593385</td><td>70</td><td>Plaques</td><td>627</td><td>FALSE</td><td>757 </td></tr>
	<tr><td>GSM892723</td><td>5.985059</td><td>7.757325</td><td>7.717067</td><td>6.649571</td><td>7.908695</td><td>6.427296</td><td>5.216108</td><td>7.698404</td><td>6.671056</td><td>⋯</td><td>4.875154</td><td>5.629678</td><td>7.004463</td><td> 9.087836</td><td>5.865254</td><td>69</td><td>Plaques</td><td>629</td><td>FALSE</td><td>2888</td></tr>
	<tr><td>GSM892724</td><td>6.506324</td><td>7.006040</td><td>7.831571</td><td>7.516133</td><td>8.455084</td><td>6.612935</td><td>5.514891</td><td>7.436909</td><td>6.393400</td><td>⋯</td><td>5.294410</td><td>6.641917</td><td>7.198413</td><td>10.231938</td><td>6.875560</td><td>71</td><td>Plaques</td><td>630</td><td>TRUE </td><td>423 </td></tr>
	<tr><td>GSM892725</td><td>6.652975</td><td>7.659634</td><td>7.798400</td><td>7.000333</td><td>8.426011</td><td>6.391030</td><td>5.312481</td><td>7.291015</td><td>6.407956</td><td>⋯</td><td>5.696531</td><td>6.408603</td><td>6.890606</td><td> 8.721910</td><td>6.641521</td><td>68</td><td>Plaques</td><td>635</td><td>TRUE </td><td>122 </td></tr>
	<tr><td>GSM892726</td><td>6.461473</td><td>7.817090</td><td>7.769957</td><td>7.375247</td><td>8.326782</td><td>6.366425</td><td>5.738586</td><td>8.004114</td><td>6.113817</td><td>⋯</td><td>6.532482</td><td>6.840590</td><td>7.625847</td><td> 9.777910</td><td>7.469105</td><td>49</td><td>Plaques</td><td>638</td><td>FALSE</td><td>1890</td></tr>
	<tr><td>GSM892727</td><td>6.792632</td><td>7.382763</td><td>7.766917</td><td>7.208697</td><td>8.617717</td><td>6.863750</td><td>6.445228</td><td>7.949494</td><td>5.516677</td><td>⋯</td><td>6.167560</td><td>6.555311</td><td>7.949245</td><td> 8.958348</td><td>6.625944</td><td>76</td><td>Plaques</td><td>642</td><td>FALSE</td><td>1907</td></tr>
	<tr><td>GSM892728</td><td>6.489961</td><td>8.013562</td><td>7.666948</td><td>7.024002</td><td>8.189994</td><td>6.037001</td><td>6.772912</td><td>7.672432</td><td>5.749248</td><td>⋯</td><td>6.199361</td><td>6.614349</td><td>8.226501</td><td> 8.897003</td><td>6.940988</td><td>79</td><td>Plaques</td><td>645</td><td>FALSE</td><td>2397</td></tr>
	<tr><td>GSM892729</td><td>6.545018</td><td>8.164848</td><td>7.761506</td><td>7.121277</td><td>8.182803</td><td>6.789162</td><td>6.954836</td><td>7.498273</td><td>5.957108</td><td>⋯</td><td>6.830631</td><td>7.143977</td><td>7.602906</td><td> 7.395089</td><td>7.589764</td><td>82</td><td>Plaques</td><td>646</td><td>FALSE</td><td>686 </td></tr>
	<tr><td>GSM892730</td><td>6.044620</td><td>8.188033</td><td>7.755512</td><td>6.982878</td><td>8.267727</td><td>6.048832</td><td>6.574394</td><td>7.617016</td><td>6.013510</td><td>⋯</td><td>6.182800</td><td>6.735401</td><td>8.260799</td><td> 9.235865</td><td>6.998553</td><td>81</td><td>Plaques</td><td>651</td><td>FALSE</td><td>2997</td></tr>
	<tr><td>GSM892731</td><td>6.903474</td><td>8.426809</td><td>7.934202</td><td>7.236329</td><td>8.252199</td><td>6.443834</td><td>7.058370</td><td>7.553470</td><td>5.928988</td><td>⋯</td><td>6.420266</td><td>7.184801</td><td>7.625527</td><td> 7.338900</td><td>7.614750</td><td>NA</td><td>Plaques</td><td>654</td><td>FALSE</td><td>652 </td></tr>
	<tr><td>GSM892732</td><td>6.232739</td><td>8.138495</td><td>7.851660</td><td>7.349124</td><td>8.364646</td><td>6.412122</td><td>5.993294</td><td>8.037159</td><td>6.694239</td><td>⋯</td><td>6.445142</td><td>6.852658</td><td>7.347422</td><td>10.269789</td><td>7.473535</td><td>81</td><td>Plaques</td><td>659</td><td>TRUE </td><td>796 </td></tr>
	<tr><td>GSM892733</td><td>6.655518</td><td>7.966378</td><td>7.698371</td><td>7.400564</td><td>8.435264</td><td>6.497104</td><td>6.140997</td><td>8.035092</td><td>5.929801</td><td>⋯</td><td>6.717390</td><td>6.863758</td><td>7.543075</td><td>10.238442</td><td>7.583005</td><td>51</td><td>Plaques</td><td>664</td><td>FALSE</td><td>1186</td></tr>
	<tr><td>GSM892734</td><td>6.699976</td><td>8.217901</td><td>6.900384</td><td>7.340885</td><td>8.357698</td><td>5.696248</td><td>5.470826</td><td>7.486252</td><td>5.895026</td><td>⋯</td><td>6.542262</td><td>6.907784</td><td>7.632727</td><td> 9.983087</td><td>6.828129</td><td>82</td><td>Plaques</td><td>669</td><td>TRUE </td><td>516 </td></tr>
	<tr><td>GSM892735</td><td>6.635070</td><td>7.519029</td><td>7.757553</td><td>7.251901</td><td>8.443706</td><td>6.477393</td><td>5.618210</td><td>7.735360</td><td>6.281759</td><td>⋯</td><td>6.214451</td><td>6.824885</td><td>7.371379</td><td> 9.940597</td><td>6.977336</td><td>56</td><td>Plaques</td><td>675</td><td>FALSE</td><td>2223</td></tr>
	<tr><td>GSM892736</td><td>6.688627</td><td>7.559205</td><td>7.848844</td><td>7.423926</td><td>8.409471</td><td>6.352663</td><td>5.792368</td><td>7.918395</td><td>5.788065</td><td>⋯</td><td>6.177181</td><td>6.458729</td><td>7.380818</td><td> 9.878802</td><td>6.996492</td><td>60</td><td>Plaques</td><td>677</td><td>FALSE</td><td>1347</td></tr>
	<tr><td>GSM892737</td><td>6.461210</td><td>8.076686</td><td>8.114321</td><td>6.987587</td><td>8.770922</td><td>6.336274</td><td>6.229084</td><td>7.698160</td><td>5.769223</td><td>⋯</td><td>6.213077</td><td>6.658322</td><td>8.189443</td><td> 9.201405</td><td>7.541182</td><td>79</td><td>Plaques</td><td>678</td><td>FALSE</td><td>2123</td></tr>
	<tr><td>GSM892738</td><td>6.840229</td><td>7.750591</td><td>7.479951</td><td>7.036618</td><td>8.441403</td><td>6.394248</td><td>6.980082</td><td>7.760172</td><td>5.634933</td><td>⋯</td><td>6.510829</td><td>6.877735</td><td>8.362280</td><td> 8.578688</td><td>7.402224</td><td>69</td><td>Plaques</td><td>680</td><td>FALSE</td><td>1792</td></tr>
	<tr><td>GSM892739</td><td>6.334797</td><td>8.555989</td><td>7.770727</td><td>7.252122</td><td>8.212847</td><td>6.438389</td><td>7.001591</td><td>7.465521</td><td>5.946766</td><td>⋯</td><td>6.907458</td><td>7.112260</td><td>7.348922</td><td> 7.519751</td><td>7.430321</td><td>85</td><td>Plaques</td><td>681</td><td>TRUE </td><td>103 </td></tr>
	<tr><td>GSM892740</td><td>6.122527</td><td>7.977229</td><td>7.734951</td><td>7.289281</td><td>8.308845</td><td>6.052669</td><td>5.713552</td><td>7.948864</td><td>6.569940</td><td>⋯</td><td>6.144444</td><td>6.756438</td><td>7.566873</td><td>10.409470</td><td>7.066382</td><td>73</td><td>Plaques</td><td>686</td><td>FALSE</td><td>1479</td></tr>
</tbody>
</table>




```R
 cor.test(gene_time$AARS,as.numeric(gene_time$Time)) 
cor.test(gene_time$AARS,as.numeric(gene_time$Time))  %>% 
  str()
demo  <- cor.test(gene_time$AARS,as.numeric(gene_time$Time)) 
```


    
    	Pearson's product-moment correlation
    
    data:  gene_time$AARS and as.numeric(gene_time$Time)
    t = 3.8513, df = 221, p-value = 0.0001539
    alternative hypothesis: true correlation is not equal to 0
    95 percent confidence interval:
     0.1234761 0.3699722
    sample estimates:
          cor 
    0.2507848 



    List of 9
     $ statistic  : Named num 3.85
      ..- attr(*, "names")= chr "t"
     $ parameter  : Named int 221
      ..- attr(*, "names")= chr "df"
     $ p.value    : num 0.000154
     $ estimate   : Named num 0.251
      ..- attr(*, "names")= chr "cor"
     $ null.value : Named num 0
      ..- attr(*, "names")= chr "correlation"
     $ alternative: chr "two.sided"
     $ method     : chr "Pearson's product-moment correlation"
     $ data.name  : chr "gene_time$AARS and as.numeric(gene_time$Time)"
     $ conf.int   : num [1:2] 0.123 0.37
      ..- attr(*, "conf.level")= num 0.95
     - attr(*, "class")= chr "htest"



```R
demo[['p.value']]
```


0.000153937933091306



```R

pearson_tibble  <- tibble(Gene = 'AACS',pearson_index = cor(gene_time$AACS,as.numeric(gene_time$Time), method = 'pearson'), 
pearson_value = cor.test(gene_time$AACS,as.numeric(gene_time$Time))$p.value)

```


```R

for (i in 2:length(selected_gene)){

tibble_gene  <- tibble(Gene = selected_gene[i],
                       pearson_index = cor(gene_time[,selected_gene[i]] %>% pull(),as.numeric(gene_time$Time), method = 'pearson'),
                       pearson_value = cor.test(gene_time[,selected_gene[i]]%>% pull(),as.numeric(gene_time$Time))$p.value)
 pearson_tibble  <- bind_rows(pearson_tibble,tibble_gene)   
}
```


```R
pearson_tibble  %>% 
  filter(pearson_value < 0.05)  %>% 
    arrange(pearson_index)
```


<table class="dataframe">
<caption>A tibble: 2970 × 3</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>pearson_index</th><th scope=col>pearson_value</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>KIAA0753                                </td><td>-0.4086887</td><td>2.179628e-10</td></tr>
	<tr><td>EI24                                    </td><td>-0.3745263</td><td>7.790535e-09</td></tr>
	<tr><td>ID2 /// ID2B                            </td><td>-0.3568507</td><td>4.241849e-08</td></tr>
	<tr><td>SUDS3                                   </td><td>-0.3558234</td><td>4.666316e-08</td></tr>
	<tr><td>ZNF419                                  </td><td>-0.3488254</td><td>8.856138e-08</td></tr>
	<tr><td>ZHX1                                    </td><td>-0.3481078</td><td>9.449307e-08</td></tr>
	<tr><td>SUGP1                                   </td><td>-0.3476511</td><td>9.846453e-08</td></tr>
	<tr><td>FHOD1                                   </td><td>-0.3474296</td><td>1.004488e-07</td></tr>
	<tr><td>CAPZA2                                  </td><td>-0.3450190</td><td>1.246771e-07</td></tr>
	<tr><td>UGP2                                    </td><td>-0.3440390</td><td>1.360543e-07</td></tr>
	<tr><td>BCL2L2-PABPN1 /// PABPN1                </td><td>-0.3423130</td><td>1.585597e-07</td></tr>
	<tr><td>COPG1                                   </td><td>-0.3417881</td><td>1.660839e-07</td></tr>
	<tr><td>MIRLET7D                                </td><td>-0.3413357</td><td>1.728450e-07</td></tr>
	<tr><td>LOC150776 /// SMPD4                     </td><td>-0.3412929</td><td>1.734981e-07</td></tr>
	<tr><td>HNRNPU-AS1                              </td><td>-0.3408608</td><td>1.802243e-07</td></tr>
	<tr><td>KPNA6                                   </td><td>-0.3401516</td><td>1.918118e-07</td></tr>
	<tr><td>PPIL4                                   </td><td>-0.3396256</td><td>2.008648e-07</td></tr>
	<tr><td>CTNNB1                                  </td><td>-0.3385975</td><td>2.197565e-07</td></tr>
	<tr><td>TAF1D                                   </td><td>-0.3365936</td><td>2.615953e-07</td></tr>
	<tr><td>AMFR                                    </td><td>-0.3345207</td><td>3.128671e-07</td></tr>
	<tr><td>LOC220729 /// SDHA /// SDHAP1 /// SDHAP2</td><td>-0.3343838</td><td>3.165725e-07</td></tr>
	<tr><td>GOLGA1                                  </td><td>-0.3338451</td><td>3.315682e-07</td></tr>
	<tr><td>EAF1                                    </td><td>-0.3318437</td><td>3.934643e-07</td></tr>
	<tr><td>HERC2P2 /// HERC2P9                     </td><td>-0.3290204</td><td>4.998855e-07</td></tr>
	<tr><td>CISD1                                   </td><td>-0.3282547</td><td>5.331957e-07</td></tr>
	<tr><td>RLF                                     </td><td>-0.3280828</td><td>5.409593e-07</td></tr>
	<tr><td>SWI5                                    </td><td>-0.3261869</td><td>6.341072e-07</td></tr>
	<tr><td>BAZ2B                                   </td><td>-0.3253241</td><td>6.814089e-07</td></tr>
	<tr><td>GNE                                     </td><td>-0.3251192</td><td>6.931282e-07</td></tr>
	<tr><td>PGAP2                                   </td><td>-0.3246960</td><td>7.179530e-07</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>PFKL         </td><td>0.3042432</td><td>3.692719e-06</td></tr>
	<tr><td>PPP2R1A      </td><td>0.3043567</td><td>3.660521e-06</td></tr>
	<tr><td>ARMC5        </td><td>0.3063851</td><td>3.128409e-06</td></tr>
	<tr><td>CALM1        </td><td>0.3108287</td><td>2.208465e-06</td></tr>
	<tr><td>SETD8        </td><td>0.3114557</td><td>2.101599e-06</td></tr>
	<tr><td>RAD23A       </td><td>0.3116660</td><td>2.066890e-06</td></tr>
	<tr><td>DNAJC3       </td><td>0.3150232</td><td>1.581466e-06</td></tr>
	<tr><td>PNPLA2       </td><td>0.3173938</td><td>1.306511e-06</td></tr>
	<tr><td>NRIP1        </td><td>0.3218751</td><td>9.064667e-07</td></tr>
	<tr><td>YWHAE        </td><td>0.3219642</td><td>8.998504e-07</td></tr>
	<tr><td>IRF2BP1      </td><td>0.3236908</td><td>7.803623e-07</td></tr>
	<tr><td>UBA1         </td><td>0.3254337</td><td>6.752221e-07</td></tr>
	<tr><td>VMP1         </td><td>0.3258950</td><td>6.497454e-07</td></tr>
	<tr><td>KLHDC3       </td><td>0.3294664</td><td>4.814081e-07</td></tr>
	<tr><td>DRAP1        </td><td>0.3301479</td><td>4.544371e-07</td></tr>
	<tr><td>CD74         </td><td>0.3306902</td><td>4.340173e-07</td></tr>
	<tr><td>AP2M1        </td><td>0.3321197</td><td>3.843133e-07</td></tr>
	<tr><td>PRKCSH       </td><td>0.3365815</td><td>2.618681e-07</td></tr>
	<tr><td>PPIA         </td><td>0.3384597</td><td>2.224145e-07</td></tr>
	<tr><td>TOMM40       </td><td>0.3385098</td><td>2.214452e-07</td></tr>
	<tr><td>SEC22B       </td><td>0.3429202</td><td>1.502619e-07</td></tr>
	<tr><td>TMEM214      </td><td>0.3473144</td><td>1.014949e-07</td></tr>
	<tr><td>CTC-425F1.4  </td><td>0.3480781</td><td>9.474601e-08</td></tr>
	<tr><td>TRAPPC1      </td><td>0.3515733</td><td>6.898916e-08</td></tr>
	<tr><td>CFL1         </td><td>0.3556685</td><td>4.733778e-08</td></tr>
	<tr><td>PCGF1        </td><td>0.3669461</td><td>1.631720e-08</td></tr>
	<tr><td>BANF1        </td><td>0.3684243</td><td>1.414764e-08</td></tr>
	<tr><td>RP11-403P17.3</td><td>0.3813622</td><td>3.934151e-09</td></tr>
	<tr><td>Y16709       </td><td>0.3858911</td><td>2.480114e-09</td></tr>
	<tr><td>MOSPD3       </td><td>0.4026712</td><td>4.216687e-10</td></tr>
</tbody>
</table>




```R
normalized_gset_mean['KIAA0753',]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples')
```


<table class="dataframe">
<caption>A data.frame: 223 × 7</caption>
<thead>
	<tr><th scope=col>Samples</th><th scope=col>KIAA0753</th><th scope=col>Ages</th><th scope=col>Type</th><th scope=col>Patients</th><th scope=col>Outcome</th><th scope=col>Time</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM892518</td><td>7.061467</td><td>68</td><td>PBMC</td><td>279</td><td>FALSE</td><td>1738</td></tr>
	<tr><td>GSM892519</td><td>7.197194</td><td>78</td><td>PBMC</td><td>281</td><td>FALSE</td><td>697 </td></tr>
	<tr><td>GSM892520</td><td>7.029200</td><td>84</td><td>PBMC</td><td>282</td><td>TRUE </td><td>2   </td></tr>
	<tr><td>GSM892521</td><td>7.007012</td><td>80</td><td>PBMC</td><td>284</td><td>FALSE</td><td>1767</td></tr>
	<tr><td>GSM892522</td><td>7.376876</td><td>81</td><td>PBMC</td><td>285</td><td>FALSE</td><td>1326</td></tr>
	<tr><td>GSM892523</td><td>6.918401</td><td>50</td><td>PBMC</td><td>286</td><td>FALSE</td><td>1427</td></tr>
	<tr><td>GSM892524</td><td>6.146062</td><td>75</td><td>PBMC</td><td>288</td><td>TRUE </td><td>7   </td></tr>
	<tr><td>GSM892525</td><td>7.025291</td><td>68</td><td>PBMC</td><td>292</td><td>FALSE</td><td>1422</td></tr>
	<tr><td>GSM892526</td><td>6.933379</td><td>73</td><td>PBMC</td><td>296</td><td>FALSE</td><td>1320</td></tr>
	<tr><td>GSM892527</td><td>7.228920</td><td>75</td><td>PBMC</td><td>298</td><td>FALSE</td><td>689 </td></tr>
	<tr><td>GSM892528</td><td>7.108720</td><td>64</td><td>PBMC</td><td>301</td><td>FALSE</td><td>676 </td></tr>
	<tr><td>GSM892529</td><td>7.008571</td><td>66</td><td>PBMC</td><td>302</td><td>FALSE</td><td>1710</td></tr>
	<tr><td>GSM892530</td><td>6.716637</td><td>71</td><td>PBMC</td><td>306</td><td>FALSE</td><td>1688</td></tr>
	<tr><td>GSM892531</td><td>6.923172</td><td>56</td><td>PBMC</td><td>311</td><td>FALSE</td><td>1711</td></tr>
	<tr><td>GSM892532</td><td>6.906010</td><td>76</td><td>PBMC</td><td>313</td><td>FALSE</td><td>751 </td></tr>
	<tr><td>GSM892533</td><td>6.401335</td><td>78</td><td>PBMC</td><td>318</td><td>FALSE</td><td>1680</td></tr>
	<tr><td>GSM892534</td><td>6.121687</td><td>69</td><td>PBMC</td><td>327</td><td>FALSE</td><td>1876</td></tr>
	<tr><td>GSM892535</td><td>7.057700</td><td>78</td><td>PBMC</td><td>333</td><td>FALSE</td><td>627 </td></tr>
	<tr><td>GSM892536</td><td>6.368621</td><td>66</td><td>PBMC</td><td>334</td><td>FALSE</td><td>1779</td></tr>
	<tr><td>GSM892537</td><td>6.944781</td><td>84</td><td>PBMC</td><td>337</td><td>FALSE</td><td>1542</td></tr>
	<tr><td>GSM892538</td><td>7.081838</td><td>79</td><td>PBMC</td><td>340</td><td>FALSE</td><td>777 </td></tr>
	<tr><td>GSM892539</td><td>7.221115</td><td>78</td><td>PBMC</td><td>345</td><td>FALSE</td><td>687 </td></tr>
	<tr><td>GSM892540</td><td>7.074890</td><td>65</td><td>PBMC</td><td>355</td><td>FALSE</td><td>739 </td></tr>
	<tr><td>GSM892541</td><td>7.021139</td><td>77</td><td>PBMC</td><td>356</td><td>TRUE </td><td>2   </td></tr>
	<tr><td>GSM892542</td><td>6.574691</td><td>79</td><td>PBMC</td><td>359</td><td>FALSE</td><td>1674</td></tr>
	<tr><td>GSM892543</td><td>7.072705</td><td>72</td><td>PBMC</td><td>363</td><td>FALSE</td><td>1649</td></tr>
	<tr><td>GSM892544</td><td>6.787240</td><td>63</td><td>PBMC</td><td>364</td><td>FALSE</td><td>1842</td></tr>
	<tr><td>GSM892545</td><td>6.918430</td><td>81</td><td>PBMC</td><td>365</td><td>TRUE </td><td>922 </td></tr>
	<tr><td>GSM892546</td><td>7.051250</td><td>59</td><td>PBMC</td><td>366</td><td>TRUE </td><td>1   </td></tr>
	<tr><td>GSM892547</td><td>6.747977</td><td>69</td><td>PBMC</td><td>367</td><td>FALSE</td><td>1886</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>GSM892711</td><td>6.188320</td><td>61</td><td>Plaques</td><td>587</td><td>FALSE</td><td>2083</td></tr>
	<tr><td>GSM892712</td><td>6.919713</td><td>48</td><td>Plaques</td><td>588</td><td>FALSE</td><td>1704</td></tr>
	<tr><td>GSM892713</td><td>6.850619</td><td>77</td><td>Plaques</td><td>591</td><td>TRUE </td><td>21  </td></tr>
	<tr><td>GSM892714</td><td>6.351322</td><td>66</td><td>Plaques</td><td>597</td><td>FALSE</td><td>2115</td></tr>
	<tr><td>GSM892715</td><td>7.250307</td><td>73</td><td>Plaques</td><td>603</td><td>FALSE</td><td>1306</td></tr>
	<tr><td>GSM892716</td><td>6.706782</td><td>77</td><td>Plaques</td><td>605</td><td>FALSE</td><td>1870</td></tr>
	<tr><td>GSM892717</td><td>7.060458</td><td>62</td><td>Plaques</td><td>606</td><td>FALSE</td><td>1173</td></tr>
	<tr><td>GSM892718</td><td>6.567417</td><td>82</td><td>Plaques</td><td>618</td><td>FALSE</td><td>395 </td></tr>
	<tr><td>GSM892719</td><td>6.364502</td><td>80</td><td>Plaques</td><td>619</td><td>FALSE</td><td>2619</td></tr>
	<tr><td>GSM892720</td><td>6.846550</td><td>68</td><td>Plaques</td><td>624</td><td>FALSE</td><td>1752</td></tr>
	<tr><td>GSM892721</td><td>7.335131</td><td>70</td><td>Plaques</td><td>625</td><td>FALSE</td><td>1207</td></tr>
	<tr><td>GSM892722</td><td>7.225115</td><td>70</td><td>Plaques</td><td>627</td><td>FALSE</td><td>757 </td></tr>
	<tr><td>GSM892723</td><td>5.988878</td><td>69</td><td>Plaques</td><td>629</td><td>FALSE</td><td>2888</td></tr>
	<tr><td>GSM892724</td><td>7.145711</td><td>71</td><td>Plaques</td><td>630</td><td>TRUE </td><td>423 </td></tr>
	<tr><td>GSM892725</td><td>6.552802</td><td>68</td><td>Plaques</td><td>635</td><td>TRUE </td><td>122 </td></tr>
	<tr><td>GSM892726</td><td>7.030119</td><td>49</td><td>Plaques</td><td>638</td><td>FALSE</td><td>1890</td></tr>
	<tr><td>GSM892727</td><td>6.653690</td><td>76</td><td>Plaques</td><td>642</td><td>FALSE</td><td>1907</td></tr>
	<tr><td>GSM892728</td><td>6.352171</td><td>79</td><td>Plaques</td><td>645</td><td>FALSE</td><td>2397</td></tr>
	<tr><td>GSM892729</td><td>7.355985</td><td>82</td><td>Plaques</td><td>646</td><td>FALSE</td><td>686 </td></tr>
	<tr><td>GSM892730</td><td>6.245331</td><td>81</td><td>Plaques</td><td>651</td><td>FALSE</td><td>2997</td></tr>
	<tr><td>GSM892731</td><td>7.308389</td><td>NA</td><td>Plaques</td><td>654</td><td>FALSE</td><td>652 </td></tr>
	<tr><td>GSM892732</td><td>7.145004</td><td>81</td><td>Plaques</td><td>659</td><td>TRUE </td><td>796 </td></tr>
	<tr><td>GSM892733</td><td>7.144392</td><td>51</td><td>Plaques</td><td>664</td><td>FALSE</td><td>1186</td></tr>
	<tr><td>GSM892734</td><td>6.827092</td><td>82</td><td>Plaques</td><td>669</td><td>TRUE </td><td>516 </td></tr>
	<tr><td>GSM892735</td><td>6.831446</td><td>56</td><td>Plaques</td><td>675</td><td>FALSE</td><td>2223</td></tr>
	<tr><td>GSM892736</td><td>7.018951</td><td>60</td><td>Plaques</td><td>677</td><td>FALSE</td><td>1347</td></tr>
	<tr><td>GSM892737</td><td>6.417925</td><td>79</td><td>Plaques</td><td>678</td><td>FALSE</td><td>2123</td></tr>
	<tr><td>GSM892738</td><td>6.811850</td><td>69</td><td>Plaques</td><td>680</td><td>FALSE</td><td>1792</td></tr>
	<tr><td>GSM892739</td><td>7.482747</td><td>85</td><td>Plaques</td><td>681</td><td>TRUE </td><td>103 </td></tr>
	<tr><td>GSM892740</td><td>7.097211</td><td>73</td><td>Plaques</td><td>686</td><td>FALSE</td><td>1479</td></tr>
</tbody>
</table>




```R
######select the plaues
plaque_info  <- data.table::fread('~/AS_Plaques/Pos_Add_Description.csv')   %>% select('Gene','degree_pos')
```


```R
head(plaque_info)
```


<table class="dataframe">
<caption>A data.table: 6 × 2</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>degree_pos</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>PARVB</td><td>6</td></tr>
	<tr><td>LY6E </td><td>6</td></tr>
	<tr><td>LY6E </td><td>6</td></tr>
	<tr><td>CXCL1</td><td>6</td></tr>
	<tr><td>DAGLA</td><td>6</td></tr>
	<tr><td>PTK2B</td><td>6</td></tr>
</tbody>
</table>




```R
pearson_tibble  %>% 
  filter(pearson_value < 0.05)  %>% 
    arrange(pearson_index) %>% 
     left_join(plaque_info, by = 'Gene') %>% 
      filter(Gene =='UPP1')
```


<table class="dataframe">
<caption>A tibble: 1 × 4</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>pearson_index</th><th scope=col>pearson_value</th><th scope=col>degree_pos</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>UPP1</td><td>-0.1996363</td><td>0.002747554</td><td>5</td></tr>
</tbody>
</table>




```R
pearson_tibble  %>% 
  filter(pearson_value < 0.05)  %>% 
    arrange(pearson_index) %>% 
     left_join(plaque_info, by = 'Gene') %>%   
       na.omit() %>% 
        filter(degree_pos > 2) 
```


<table class="dataframe">
<caption>A tibble: 611 × 4</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>pearson_index</th><th scope=col>pearson_value</th><th scope=col>degree_pos</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>FHOD1    </td><td>-0.3474296</td><td>1.004488e-07</td><td>4</td></tr>
	<tr><td>EAF1     </td><td>-0.3318437</td><td>3.934643e-07</td><td>3</td></tr>
	<tr><td>PGAP2    </td><td>-0.3246960</td><td>7.179530e-07</td><td>3</td></tr>
	<tr><td>TYMP     </td><td>-0.3088560</td><td>2.579496e-06</td><td>3</td></tr>
	<tr><td>SLC35D2  </td><td>-0.2975581</td><td>6.144573e-06</td><td>3</td></tr>
	<tr><td>SLC35D2  </td><td>-0.2975581</td><td>6.144573e-06</td><td>3</td></tr>
	<tr><td>ATP1A1   </td><td>-0.2971546</td><td>6.333787e-06</td><td>3</td></tr>
	<tr><td>PSMD3    </td><td>-0.2957989</td><td>7.010981e-06</td><td>3</td></tr>
	<tr><td>BAK1     </td><td>-0.2953944</td><td>7.225968e-06</td><td>5</td></tr>
	<tr><td>TNFAIP8L1</td><td>-0.2937222</td><td>8.183148e-06</td><td>3</td></tr>
	<tr><td>NRD1     </td><td>-0.2862905</td><td>1.409081e-05</td><td>3</td></tr>
	<tr><td>IFI44    </td><td>-0.2853742</td><td>1.505146e-05</td><td>4</td></tr>
	<tr><td>TAZ      </td><td>-0.2843991</td><td>1.614158e-05</td><td>4</td></tr>
	<tr><td>MAP3K8   </td><td>-0.2807482</td><td>2.092464e-05</td><td>4</td></tr>
	<tr><td>ZBED1    </td><td>-0.2738008</td><td>3.394748e-05</td><td>4</td></tr>
	<tr><td>STARD3   </td><td>-0.2727789</td><td>3.641197e-05</td><td>5</td></tr>
	<tr><td>FAM91A1  </td><td>-0.2706074</td><td>4.221986e-05</td><td>3</td></tr>
	<tr><td>RPN1     </td><td>-0.2666418</td><td>5.514286e-05</td><td>3</td></tr>
	<tr><td>DUSP5    </td><td>-0.2666002</td><td>5.529643e-05</td><td>4</td></tr>
	<tr><td>SRBD1    </td><td>-0.2656481</td><td>5.892046e-05</td><td>3</td></tr>
	<tr><td>DNAJB11  </td><td>-0.2602397</td><td>8.412394e-05</td><td>3</td></tr>
	<tr><td>TARBP1   </td><td>-0.2600434</td><td>8.520597e-05</td><td>4</td></tr>
	<tr><td>APH1B    </td><td>-0.2599494</td><td>8.572867e-05</td><td>4</td></tr>
	<tr><td>ABCF3    </td><td>-0.2571738</td><td>1.025926e-04</td><td>4</td></tr>
	<tr><td>CPT2     </td><td>-0.2558670</td><td>1.115662e-04</td><td>3</td></tr>
	<tr><td>SHMT2    </td><td>-0.2558523</td><td>1.116716e-04</td><td>4</td></tr>
	<tr><td>PPP1R16B </td><td>-0.2537194</td><td>1.279283e-04</td><td>4</td></tr>
	<tr><td>ARRDC3   </td><td>-0.2530899</td><td>1.331343e-04</td><td>3</td></tr>
	<tr><td>TSPAN14  </td><td>-0.2505182</td><td>1.565311e-04</td><td>4</td></tr>
	<tr><td>LRPAP1   </td><td>-0.2504085</td><td>1.576102e-04</td><td>3</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>RILP   </td><td>0.2306793</td><td>5.158914e-04</td><td>5</td></tr>
	<tr><td>SPI1   </td><td>0.2309619</td><td>5.075520e-04</td><td>4</td></tr>
	<tr><td>FURIN  </td><td>0.2331572</td><td>4.468915e-04</td><td>4</td></tr>
	<tr><td>UQCRC1 </td><td>0.2377998</td><td>3.400844e-04</td><td>3</td></tr>
	<tr><td>ITGA5  </td><td>0.2381833</td><td>3.324189e-04</td><td>3</td></tr>
	<tr><td>S100A11</td><td>0.2401140</td><td>2.962032e-04</td><td>5</td></tr>
	<tr><td>CDK4   </td><td>0.2425096</td><td>2.563699e-04</td><td>3</td></tr>
	<tr><td>SNX27  </td><td>0.2429238</td><td>2.500104e-04</td><td>4</td></tr>
	<tr><td>GUSB   </td><td>0.2453292</td><td>2.158923e-04</td><td>4</td></tr>
	<tr><td>PPP2R4 </td><td>0.2499176</td><td>1.625234e-04</td><td>3</td></tr>
	<tr><td>ZNF385A</td><td>0.2567677</td><td>1.053064e-04</td><td>4</td></tr>
	<tr><td>EPHX1  </td><td>0.2590461</td><td>9.090812e-05</td><td>4</td></tr>
	<tr><td>EMP3   </td><td>0.2626576</td><td>7.181112e-05</td><td>4</td></tr>
	<tr><td>G6PC3  </td><td>0.2665117</td><td>5.562416e-05</td><td>4</td></tr>
	<tr><td>COMT   </td><td>0.2816340</td><td>1.965422e-05</td><td>4</td></tr>
	<tr><td>LGALS9 </td><td>0.2877423</td><td>1.268671e-05</td><td>5</td></tr>
	<tr><td>CSF1R  </td><td>0.2908048</td><td>1.014758e-05</td><td>5</td></tr>
	<tr><td>HECTD3 </td><td>0.2912009</td><td>9.856765e-06</td><td>4</td></tr>
	<tr><td>ALG3   </td><td>0.2923370</td><td>9.065945e-06</td><td>3</td></tr>
	<tr><td>PTGES2 </td><td>0.2974479</td><td>6.195704e-06</td><td>3</td></tr>
	<tr><td>PFKL   </td><td>0.3042432</td><td>3.692719e-06</td><td>4</td></tr>
	<tr><td>NRIP1  </td><td>0.3218751</td><td>9.064667e-07</td><td>3</td></tr>
	<tr><td>VMP1   </td><td>0.3258950</td><td>6.497454e-07</td><td>3</td></tr>
	<tr><td>CD74   </td><td>0.3306902</td><td>4.340173e-07</td><td>5</td></tr>
	<tr><td>AP2M1  </td><td>0.3321197</td><td>3.843133e-07</td><td>4</td></tr>
	<tr><td>PRKCSH </td><td>0.3365815</td><td>2.618681e-07</td><td>3</td></tr>
	<tr><td>TOMM40 </td><td>0.3385098</td><td>2.214452e-07</td><td>4</td></tr>
	<tr><td>TRAPPC1</td><td>0.3515733</td><td>6.898916e-08</td><td>4</td></tr>
	<tr><td>CFL1   </td><td>0.3556685</td><td>4.733778e-08</td><td>4</td></tr>
	<tr><td>MOSPD3 </td><td>0.4026712</td><td>4.216687e-10</td><td>3</td></tr>
</tbody>
</table>




```R

plot_pearson  <- function(x){
p <- normalized_gset_mean[x,]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples') %>% 
      as_tibble() %>% 
        select(x,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = x, xlab = 'Time',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
 return(p)
}

```


```R
plot_pearson('FHOD1')
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.8_Survival_Curve_Filter_Gene_files/Step11.8_Survival_Curve_Filter_Gene_32_1.png)
    



```R
normalized_gset_mean['KIAA0753',]  %>% 
  t()  %>% 
   as.data.frame()  %>% 
    rownames_to_column(var = 'Samples')  %>% 
     left_join(pheno_group, by  ='Samples') %>% 
      as_tibble() %>% 
        select(KIAA0753,Type,Time) %>% 
        mutate(Time = as.integer(Time)) %>% 
           ggscatter(x = "Time", y = "KIAA0753", xlab = 'Time',size = 0.5, alpha = 0.5,
                   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE ) + stat_cor(method = "pearson",size=2) +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
```

    [1m[22m`geom_smooth()` using formula = 'y ~ x'



    
![png](Step11.8_Survival_Curve_Filter_Gene_files/Step11.8_Survival_Curve_Filter_Gene_33_1.png)
    

