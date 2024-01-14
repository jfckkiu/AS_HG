```R
library(tidyverse)
library(GSVA)
library(GEOquery) 
library(limma)
library(ggpubr)
```


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
<caption>A data.frame: 23521 × 29</caption>
<thead>
	<tr><th></th><th scope=col>GSM714070</th><th scope=col>GSM714071</th><th scope=col>GSM714072</th><th scope=col>GSM714073</th><th scope=col>GSM714074</th><th scope=col>GSM714075</th><th scope=col>GSM714076</th><th scope=col>GSM714077</th><th scope=col>GSM714078</th><th scope=col>GSM714079</th><th scope=col>⋯</th><th scope=col>GSM714089</th><th scope=col>GSM714090</th><th scope=col>GSM714091</th><th scope=col>GSM714092</th><th scope=col>GSM714093</th><th scope=col>GSM714094</th><th scope=col>GSM714095</th><th scope=col>GSM714096</th><th scope=col>GSM714097</th><th scope=col>GSM714098</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row></th><td>4.204701</td><td>4.178371</td><td>4.143596</td><td>4.166885</td><td>4.108082</td><td>4.151953</td><td>4.280475</td><td>4.134108</td><td>4.185679</td><td>4.244531</td><td>⋯</td><td>4.173463</td><td>4.167156</td><td>4.211431</td><td>4.262509</td><td>4.043401</td><td>4.265211</td><td>4.264777</td><td>4.102540</td><td>4.141335</td><td>4.137400</td></tr>
	<tr><th scope=row>A1BG</th><td>4.106102</td><td>3.947136</td><td>4.113395</td><td>4.180741</td><td>4.110530</td><td>4.254026</td><td>4.369771</td><td>4.076951</td><td>4.063928</td><td>4.508848</td><td>⋯</td><td>3.862815</td><td>3.999267</td><td>4.661946</td><td>4.104036</td><td>4.029722</td><td>3.869728</td><td>4.567856</td><td>4.405135</td><td>3.983564</td><td>4.302158</td></tr>
	<tr><th scope=row>A1BG-AS1</th><td>5.336627</td><td>4.811346</td><td>4.925941</td><td>5.221658</td><td>4.664895</td><td>4.686854</td><td>5.652197</td><td>5.586011</td><td>5.146922</td><td>5.239102</td><td>⋯</td><td>5.072761</td><td>4.845185</td><td>5.524900</td><td>4.978780</td><td>4.729793</td><td>4.595183</td><td>5.176554</td><td>4.593427</td><td>5.242100</td><td>5.307885</td></tr>
	<tr><th scope=row>A1CF</th><td>3.779809</td><td>3.904868</td><td>3.858832</td><td>3.853183</td><td>3.762851</td><td>3.932020</td><td>4.038929</td><td>3.807236</td><td>4.068498</td><td>3.942386</td><td>⋯</td><td>3.875870</td><td>3.891494</td><td>4.014399</td><td>4.123775</td><td>3.817804</td><td>3.979776</td><td>3.770665</td><td>3.866540</td><td>3.960822</td><td>3.764874</td></tr>
	<tr><th scope=row>A2M</th><td>8.659909</td><td>8.651674</td><td>8.674632</td><td>8.606335</td><td>8.533633</td><td>8.564443</td><td>8.463527</td><td>8.866571</td><td>8.680580</td><td>8.668496</td><td>⋯</td><td>8.576654</td><td>8.472089</td><td>8.813100</td><td>8.690898</td><td>8.534455</td><td>8.546497</td><td>8.623689</td><td>8.380991</td><td>8.776827</td><td>8.531878</td></tr>
	<tr><th scope=row>A2M-AS1</th><td>4.935758</td><td>5.013675</td><td>4.942323</td><td>4.726485</td><td>5.785828</td><td>5.953109</td><td>5.437650</td><td>4.760124</td><td>5.081644</td><td>4.970680</td><td>⋯</td><td>5.231338</td><td>5.660150</td><td>5.234488</td><td>6.537811</td><td>5.968495</td><td>5.531858</td><td>5.421051</td><td>5.815954</td><td>4.822342</td><td>5.137411</td></tr>
	<tr><th scope=row>A2ML1</th><td>3.836972</td><td>3.960079</td><td>3.888047</td><td>3.924307</td><td>3.672496</td><td>3.915952</td><td>4.035383</td><td>3.728304</td><td>3.838071</td><td>3.877560</td><td>⋯</td><td>3.930211</td><td>3.795777</td><td>4.114669</td><td>3.926175</td><td>3.878757</td><td>4.004829</td><td>3.745379</td><td>3.869285</td><td>3.773188</td><td>3.972736</td></tr>
	<tr><th scope=row>A2MP1</th><td>3.843107</td><td>3.858819</td><td>3.852384</td><td>3.738754</td><td>3.814954</td><td>4.038877</td><td>3.985710</td><td>3.593841</td><td>3.817349</td><td>3.882286</td><td>⋯</td><td>3.951803</td><td>3.757386</td><td>3.672199</td><td>4.294355</td><td>3.783103</td><td>4.139190</td><td>3.840809</td><td>3.772074</td><td>3.675566</td><td>3.788172</td></tr>
	<tr><th scope=row>A4GALT</th><td>6.176025</td><td>6.188115</td><td>6.328069</td><td>6.585833</td><td>5.934496</td><td>5.937293</td><td>6.867352</td><td>7.005228</td><td>6.591046</td><td>6.371748</td><td>⋯</td><td>6.504925</td><td>5.997439</td><td>6.814168</td><td>6.618655</td><td>6.255978</td><td>5.957907</td><td>6.441533</td><td>6.295238</td><td>6.619972</td><td>6.536432</td></tr>
	<tr><th scope=row>A4GNT</th><td>4.266891</td><td>4.314068</td><td>4.351827</td><td>4.373050</td><td>4.181610</td><td>4.327797</td><td>4.454779</td><td>4.139948</td><td>4.721252</td><td>4.530140</td><td>⋯</td><td>4.329388</td><td>4.297643</td><td>4.301049</td><td>4.355148</td><td>4.310243</td><td>4.738760</td><td>4.222349</td><td>4.063275</td><td>4.374438</td><td>4.246503</td></tr>
	<tr><th scope=row>AA06</th><td>5.829061</td><td>5.182662</td><td>5.314325</td><td>5.361708</td><td>5.080983</td><td>5.202658</td><td>5.285089</td><td>5.327272</td><td>5.468409</td><td>5.437074</td><td>⋯</td><td>5.454110</td><td>5.256406</td><td>5.383308</td><td>5.861998</td><td>5.214543</td><td>4.988837</td><td>5.503563</td><td>5.012132</td><td>5.285089</td><td>5.159861</td></tr>
	<tr><th scope=row>AAAS</th><td>7.375949</td><td>6.868985</td><td>6.987717</td><td>7.498753</td><td>7.033968</td><td>7.318524</td><td>7.677559</td><td>7.466543</td><td>7.131611</td><td>7.343791</td><td>⋯</td><td>7.188650</td><td>6.922810</td><td>7.450386</td><td>7.629798</td><td>6.924842</td><td>7.207743</td><td>7.107549</td><td>6.795420</td><td>7.434865</td><td>7.437234</td></tr>
	<tr><th scope=row>AACS</th><td>6.259346</td><td>6.290238</td><td>6.663290</td><td>7.123100</td><td>6.986371</td><td>6.911871</td><td>6.516723</td><td>6.644207</td><td>6.444059</td><td>6.242975</td><td>⋯</td><td>6.504593</td><td>6.320884</td><td>5.481974</td><td>6.361695</td><td>6.553592</td><td>6.850658</td><td>6.408710</td><td>7.325834</td><td>6.862727</td><td>6.726357</td></tr>
	<tr><th scope=row>AACSP1</th><td>3.366667</td><td>3.264495</td><td>3.304917</td><td>3.244466</td><td>3.315976</td><td>3.532193</td><td>3.526237</td><td>3.372873</td><td>3.315420</td><td>3.366667</td><td>⋯</td><td>3.329105</td><td>3.366667</td><td>3.362210</td><td>3.465694</td><td>3.478380</td><td>3.432777</td><td>3.510820</td><td>3.168870</td><td>3.298755</td><td>3.404514</td></tr>
	<tr><th scope=row>AADAC</th><td>4.633496</td><td>4.518946</td><td>4.387143</td><td>4.462769</td><td>4.142629</td><td>4.298071</td><td>4.387143</td><td>4.324682</td><td>4.387143</td><td>4.604419</td><td>⋯</td><td>4.570373</td><td>4.319873</td><td>4.875610</td><td>4.527193</td><td>4.358048</td><td>4.538196</td><td>4.357443</td><td>4.195952</td><td>4.622276</td><td>4.388126</td></tr>
	<tr><th scope=row>AADACL2</th><td>3.572257</td><td>3.415901</td><td>3.145400</td><td>3.344613</td><td>3.467589</td><td>3.385383</td><td>3.564092</td><td>3.316250</td><td>3.222415</td><td>3.431241</td><td>⋯</td><td>4.219802</td><td>4.049939</td><td>3.400620</td><td>3.258379</td><td>4.045701</td><td>4.500324</td><td>3.267457</td><td>3.346536</td><td>3.361568</td><td>3.365454</td></tr>
	<tr><th scope=row>AADACP1</th><td>3.282612</td><td>3.163265</td><td>3.282571</td><td>3.290759</td><td>3.148048</td><td>2.976001</td><td>3.239670</td><td>3.147820</td><td>3.126730</td><td>3.184131</td><td>⋯</td><td>3.203397</td><td>3.302927</td><td>3.245092</td><td>3.311634</td><td>3.133207</td><td>3.163704</td><td>3.109630</td><td>3.475713</td><td>3.101089</td><td>3.215238</td></tr>
	<tr><th scope=row>AADAT</th><td>4.699800</td><td>4.519791</td><td>4.773737</td><td>4.572053</td><td>4.994756</td><td>4.720551</td><td>4.647409</td><td>4.046638</td><td>4.309393</td><td>4.494054</td><td>⋯</td><td>4.733684</td><td>5.158550</td><td>4.873562</td><td>5.395247</td><td>4.729718</td><td>5.739162</td><td>5.057908</td><td>4.529368</td><td>4.716826</td><td>4.964841</td></tr>
	<tr><th scope=row>AAED1</th><td>7.031979</td><td>7.373937</td><td>7.119238</td><td>6.982464</td><td>7.500704</td><td>7.439378</td><td>6.230202</td><td>6.785997</td><td>6.781792</td><td>6.699520</td><td>⋯</td><td>7.061056</td><td>6.936909</td><td>6.323160</td><td>6.726016</td><td>6.921136</td><td>6.729404</td><td>7.256034</td><td>7.559638</td><td>6.520067</td><td>6.565175</td></tr>
	<tr><th scope=row>AAGAB</th><td>5.272612</td><td>5.041865</td><td>5.294904</td><td>5.139230</td><td>5.377591</td><td>5.169908</td><td>5.127413</td><td>4.912737</td><td>5.039895</td><td>5.168866</td><td>⋯</td><td>5.350644</td><td>5.192675</td><td>5.297305</td><td>5.296730</td><td>5.063447</td><td>5.008165</td><td>5.097654</td><td>5.529641</td><td>5.235556</td><td>5.073676</td></tr>
	<tr><th scope=row>AAK1</th><td>5.524248</td><td>5.504725</td><td>5.550905</td><td>5.510843</td><td>5.378284</td><td>5.448990</td><td>5.503373</td><td>5.423780</td><td>5.535307</td><td>5.580593</td><td>⋯</td><td>5.577838</td><td>5.415347</td><td>5.217837</td><td>5.612288</td><td>5.276701</td><td>5.619645</td><td>5.405663</td><td>5.347350</td><td>5.431114</td><td>5.378499</td></tr>
	<tr><th scope=row>AAMDC</th><td>5.845038</td><td>5.838212</td><td>5.649203</td><td>5.909766</td><td>6.008041</td><td>5.945756</td><td>6.052510</td><td>5.935154</td><td>5.842478</td><td>5.825416</td><td>⋯</td><td>5.888667</td><td>5.953371</td><td>6.176645</td><td>5.758899</td><td>6.059375</td><td>5.746431</td><td>6.256036</td><td>5.905339</td><td>5.933274</td><td>6.143611</td></tr>
	<tr><th scope=row>AAMP</th><td>7.806008</td><td>7.631903</td><td>7.951581</td><td>7.909716</td><td>8.039593</td><td>7.779688</td><td>7.894021</td><td>7.998238</td><td>7.813201</td><td>7.878115</td><td>⋯</td><td>7.860414</td><td>7.826723</td><td>8.024891</td><td>8.173666</td><td>7.915873</td><td>7.471637</td><td>7.735428</td><td>7.950192</td><td>7.988500</td><td>8.225056</td></tr>
	<tr><th scope=row>AANAT</th><td>3.474932</td><td>3.558082</td><td>3.643273</td><td>3.477994</td><td>3.445380</td><td>3.409878</td><td>3.832212</td><td>3.482291</td><td>3.651300</td><td>3.802784</td><td>⋯</td><td>3.470705</td><td>3.417737</td><td>3.917020</td><td>3.684232</td><td>3.551115</td><td>3.496414</td><td>3.405148</td><td>3.674841</td><td>3.290204</td><td>3.613069</td></tr>
	<tr><th scope=row>AAR2</th><td>7.050068</td><td>7.085621</td><td>7.390868</td><td>7.319545</td><td>7.019556</td><td>7.081887</td><td>7.286499</td><td>7.236298</td><td>7.470641</td><td>6.960413</td><td>⋯</td><td>7.076317</td><td>7.184909</td><td>7.416425</td><td>7.578767</td><td>7.242981</td><td>7.236705</td><td>6.888615</td><td>7.367334</td><td>7.678296</td><td>7.394968</td></tr>
	<tr><th scope=row>AARS</th><td>8.923542</td><td>8.781250</td><td>8.892499</td><td>9.205739</td><td>8.603629</td><td>8.756985</td><td>9.086555</td><td>9.084149</td><td>8.983862</td><td>9.130559</td><td>⋯</td><td>8.962659</td><td>9.036867</td><td>9.221237</td><td>9.194308</td><td>8.871197</td><td>8.907193</td><td>8.712182</td><td>8.938274</td><td>9.155971</td><td>9.290644</td></tr>
	<tr><th scope=row>AARS2</th><td>5.947106</td><td>5.938774</td><td>5.893777</td><td>6.063958</td><td>5.841475</td><td>5.796699</td><td>6.106543</td><td>5.927696</td><td>5.886258</td><td>5.873204</td><td>⋯</td><td>5.889057</td><td>5.754701</td><td>5.966599</td><td>6.052349</td><td>5.572122</td><td>5.904109</td><td>5.680229</td><td>5.728520</td><td>5.762313</td><td>5.583508</td></tr>
	<tr><th scope=row>AARSD1 /// PTGES3L /// PTGES3L-AARSD1</th><td>5.739798</td><td>6.140150</td><td>6.099193</td><td>5.963655</td><td>6.205131</td><td>6.219933</td><td>6.327263</td><td>6.068039</td><td>6.019527</td><td>6.110080</td><td>⋯</td><td>6.173645</td><td>6.291676</td><td>6.411047</td><td>6.723465</td><td>6.722855</td><td>6.240752</td><td>6.102623</td><td>6.318823</td><td>6.278466</td><td>6.312054</td></tr>
	<tr><th scope=row>AASDH</th><td>5.450204</td><td>5.915240</td><td>5.991232</td><td>5.862259</td><td>6.063213</td><td>5.910326</td><td>5.360013</td><td>5.528700</td><td>5.544879</td><td>5.516378</td><td>⋯</td><td>5.443694</td><td>5.531546</td><td>5.400134</td><td>5.004096</td><td>5.740375</td><td>5.507284</td><td>5.833667</td><td>6.094287</td><td>5.515274</td><td>5.321592</td></tr>
	<tr><th scope=row>AASDHPPT</th><td>5.455375</td><td>5.748347</td><td>6.127687</td><td>6.105790</td><td>6.069672</td><td>5.739659</td><td>5.555763</td><td>5.548367</td><td>5.318349</td><td>5.694886</td><td>⋯</td><td>5.618564</td><td>5.933225</td><td>5.322185</td><td>4.722247</td><td>6.076346</td><td>5.667254</td><td>5.718957</td><td>6.276515</td><td>5.732382</td><td>5.650236</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>ZSCAN25</th><td>5.576850</td><td>5.501148</td><td>5.566974</td><td>5.631587</td><td>5.493802</td><td>5.445467</td><td>5.331722</td><td>5.817135</td><td>5.773384</td><td>5.464429</td><td>⋯</td><td>5.591010</td><td>5.299933</td><td>4.730151</td><td>5.459925</td><td>5.580999</td><td>5.503551</td><td>5.323730</td><td>5.497026</td><td>5.826958</td><td>5.525523</td></tr>
	<tr><th scope=row>ZSCAN26</th><td>5.917053</td><td>6.282629</td><td>6.512053</td><td>6.175719</td><td>6.454005</td><td>6.766410</td><td>5.863947</td><td>5.888614</td><td>6.151991</td><td>6.063222</td><td>⋯</td><td>6.172473</td><td>6.545453</td><td>5.566559</td><td>7.068375</td><td>7.057580</td><td>6.436746</td><td>6.383092</td><td>6.378708</td><td>6.065440</td><td>5.870974</td></tr>
	<tr><th scope=row>ZSCAN29</th><td>5.372817</td><td>6.066830</td><td>6.110933</td><td>6.056618</td><td>6.244892</td><td>6.260856</td><td>5.818202</td><td>5.882775</td><td>6.491942</td><td>5.874470</td><td>⋯</td><td>6.042757</td><td>6.162003</td><td>6.033561</td><td>5.561968</td><td>5.988882</td><td>6.642746</td><td>6.086095</td><td>6.244934</td><td>6.211146</td><td>6.087452</td></tr>
	<tr><th scope=row>ZSCAN30</th><td>5.167281</td><td>5.292257</td><td>5.246287</td><td>5.157484</td><td>5.117607</td><td>5.218470</td><td>5.293157</td><td>5.158242</td><td>5.319522</td><td>5.236330</td><td>⋯</td><td>5.151904</td><td>5.082524</td><td>5.120998</td><td>5.414138</td><td>5.253158</td><td>5.334716</td><td>5.174419</td><td>5.044210</td><td>5.149220</td><td>5.166467</td></tr>
	<tr><th scope=row>ZSCAN31</th><td>4.861661</td><td>4.712895</td><td>4.854374</td><td>4.775828</td><td>4.798931</td><td>4.842801</td><td>4.446394</td><td>4.580120</td><td>4.566709</td><td>4.806338</td><td>⋯</td><td>4.897953</td><td>5.289462</td><td>4.697163</td><td>5.035382</td><td>5.791903</td><td>5.051937</td><td>5.159487</td><td>4.956967</td><td>4.749387</td><td>4.854211</td></tr>
	<tr><th scope=row>ZSCAN32</th><td>5.754153</td><td>5.700887</td><td>5.575510</td><td>5.436720</td><td>5.613704</td><td>5.812824</td><td>5.542960</td><td>5.686676</td><td>5.606083</td><td>5.842161</td><td>⋯</td><td>5.537563</td><td>5.575797</td><td>5.319181</td><td>5.750625</td><td>6.128203</td><td>5.693581</td><td>5.285557</td><td>5.749901</td><td>5.576617</td><td>6.052531</td></tr>
	<tr><th scope=row>ZSCAN4</th><td>3.261880</td><td>3.371152</td><td>3.387520</td><td>3.218107</td><td>3.252132</td><td>3.114244</td><td>3.415945</td><td>3.185195</td><td>3.223856</td><td>3.155678</td><td>⋯</td><td>3.420478</td><td>3.220099</td><td>3.387573</td><td>3.453639</td><td>3.195021</td><td>3.304791</td><td>3.270123</td><td>3.224576</td><td>3.254093</td><td>3.188966</td></tr>
	<tr><th scope=row>ZSCAN5A</th><td>4.809591</td><td>4.872012</td><td>4.804760</td><td>4.962774</td><td>4.709530</td><td>4.816773</td><td>5.178146</td><td>4.695498</td><td>4.931747</td><td>4.911126</td><td>⋯</td><td>4.652444</td><td>4.625903</td><td>5.327661</td><td>4.746256</td><td>4.776717</td><td>4.944288</td><td>4.661506</td><td>4.834602</td><td>4.685713</td><td>4.789441</td></tr>
	<tr><th scope=row>ZSCAN9</th><td>5.216342</td><td>5.257293</td><td>5.191099</td><td>5.164968</td><td>5.011562</td><td>5.157011</td><td>4.999741</td><td>5.004980</td><td>4.962485</td><td>5.151312</td><td>⋯</td><td>5.047638</td><td>4.996475</td><td>4.779671</td><td>5.375870</td><td>5.251515</td><td>5.037287</td><td>4.727732</td><td>5.023123</td><td>5.102834</td><td>5.152413</td></tr>
	<tr><th scope=row>ZSWIM1</th><td>5.513959</td><td>5.534741</td><td>5.647308</td><td>5.633362</td><td>5.632852</td><td>5.629261</td><td>5.442234</td><td>5.481038</td><td>5.712950</td><td>5.884571</td><td>⋯</td><td>5.578564</td><td>5.435245</td><td>5.513611</td><td>5.712802</td><td>5.309378</td><td>5.686290</td><td>5.785001</td><td>5.361599</td><td>5.526915</td><td>5.399139</td></tr>
	<tr><th scope=row>ZSWIM2</th><td>3.147564</td><td>2.990286</td><td>3.046064</td><td>3.053675</td><td>3.067725</td><td>2.955893</td><td>3.140944</td><td>3.008819</td><td>2.954865</td><td>3.136060</td><td>⋯</td><td>3.053675</td><td>3.079424</td><td>3.137307</td><td>3.315474</td><td>3.188288</td><td>3.154208</td><td>2.930063</td><td>3.053675</td><td>2.930063</td><td>3.165373</td></tr>
	<tr><th scope=row>ZSWIM3</th><td>5.292392</td><td>5.179029</td><td>5.422553</td><td>5.777641</td><td>5.543375</td><td>5.675699</td><td>5.543375</td><td>5.499818</td><td>5.543375</td><td>5.562396</td><td>⋯</td><td>5.575770</td><td>5.192373</td><td>5.684009</td><td>5.911968</td><td>5.661382</td><td>5.516274</td><td>5.108468</td><td>5.543375</td><td>5.479237</td><td>5.375196</td></tr>
	<tr><th scope=row>ZSWIM4</th><td>4.556306</td><td>4.433170</td><td>4.142952</td><td>4.652882</td><td>3.912856</td><td>3.995772</td><td>4.548857</td><td>4.572304</td><td>4.332467</td><td>4.530762</td><td>⋯</td><td>3.981093</td><td>4.258388</td><td>4.865989</td><td>4.092057</td><td>4.161559</td><td>4.300452</td><td>4.150107</td><td>4.146677</td><td>4.336046</td><td>4.311450</td></tr>
	<tr><th scope=row>ZSWIM5</th><td>5.092565</td><td>5.133174</td><td>5.182286</td><td>5.089342</td><td>5.133174</td><td>5.224193</td><td>5.156336</td><td>4.819692</td><td>5.006631</td><td>5.113858</td><td>⋯</td><td>5.040788</td><td>5.108242</td><td>5.233987</td><td>5.945725</td><td>5.297582</td><td>5.149006</td><td>5.230768</td><td>5.218156</td><td>5.183394</td><td>4.648839</td></tr>
	<tr><th scope=row>ZSWIM6</th><td>8.055820</td><td>8.274527</td><td>8.148936</td><td>8.247590</td><td>8.115983</td><td>7.902798</td><td>7.322463</td><td>7.913948</td><td>7.696705</td><td>8.051980</td><td>⋯</td><td>7.800950</td><td>8.028108</td><td>7.170549</td><td>6.861107</td><td>7.579658</td><td>8.021557</td><td>8.397336</td><td>8.152633</td><td>7.999835</td><td>7.613235</td></tr>
	<tr><th scope=row>ZSWIM7</th><td>5.928983</td><td>6.161880</td><td>6.330669</td><td>6.617260</td><td>6.386632</td><td>6.340583</td><td>6.357119</td><td>6.232419</td><td>6.350131</td><td>6.108196</td><td>⋯</td><td>5.924743</td><td>5.909421</td><td>5.747335</td><td>5.195109</td><td>6.297969</td><td>6.123288</td><td>5.894485</td><td>6.401027</td><td>6.163282</td><td>6.176595</td></tr>
	<tr><th scope=row>ZSWIM8</th><td>6.526152</td><td>6.304656</td><td>6.151715</td><td>6.869644</td><td>6.114023</td><td>6.213919</td><td>7.042453</td><td>6.863168</td><td>6.597746</td><td>6.664900</td><td>⋯</td><td>6.468841</td><td>6.181210</td><td>7.138317</td><td>6.612834</td><td>5.795318</td><td>6.261933</td><td>6.998414</td><td>5.774104</td><td>6.910997</td><td>6.848129</td></tr>
	<tr><th scope=row>ZUFSP</th><td>5.045251</td><td>5.524085</td><td>5.400220</td><td>5.240170</td><td>5.728763</td><td>5.420924</td><td>4.729734</td><td>5.057539</td><td>4.861343</td><td>4.821492</td><td>⋯</td><td>5.030220</td><td>5.522211</td><td>4.769830</td><td>4.783476</td><td>6.028794</td><td>5.694988</td><td>5.312554</td><td>5.768123</td><td>5.240170</td><td>4.900788</td></tr>
	<tr><th scope=row>ZW10</th><td>6.532947</td><td>6.844585</td><td>7.090090</td><td>6.638143</td><td>7.069657</td><td>6.987296</td><td>6.232262</td><td>6.453781</td><td>6.601562</td><td>6.478863</td><td>⋯</td><td>6.802123</td><td>6.805328</td><td>6.443610</td><td>6.820394</td><td>7.077828</td><td>6.863640</td><td>6.453746</td><td>7.213092</td><td>6.852016</td><td>6.532196</td></tr>
	<tr><th scope=row>ZWILCH</th><td>5.179557</td><td>5.938250</td><td>6.031706</td><td>6.018961</td><td>6.356358</td><td>6.348315</td><td>4.700116</td><td>5.300679</td><td>5.262361</td><td>5.390895</td><td>⋯</td><td>5.243686</td><td>5.399875</td><td>5.094221</td><td>5.809554</td><td>5.989448</td><td>5.905992</td><td>5.417796</td><td>6.333875</td><td>5.221244</td><td>5.140981</td></tr>
	<tr><th scope=row>ZWINT</th><td>5.187827</td><td>5.694767</td><td>6.301238</td><td>6.068179</td><td>6.721838</td><td>6.564375</td><td>3.915004</td><td>4.926547</td><td>4.941605</td><td>5.568757</td><td>⋯</td><td>4.356972</td><td>4.709569</td><td>4.074976</td><td>4.530009</td><td>4.752366</td><td>5.951338</td><td>4.748276</td><td>4.693989</td><td>4.140109</td><td>3.929971</td></tr>
	<tr><th scope=row>ZXDA</th><td>4.758527</td><td>4.995421</td><td>5.130285</td><td>4.957908</td><td>5.100419</td><td>4.627589</td><td>4.331602</td><td>4.613176</td><td>4.690691</td><td>4.780953</td><td>⋯</td><td>4.940262</td><td>5.054431</td><td>4.477671</td><td>4.679559</td><td>5.644872</td><td>4.798523</td><td>5.083370</td><td>4.871985</td><td>5.009182</td><td>4.711844</td></tr>
	<tr><th scope=row>ZXDA /// ZXDB</th><td>4.152408</td><td>4.099890</td><td>4.153377</td><td>4.007248</td><td>3.970115</td><td>3.881231</td><td>4.065842</td><td>4.090594</td><td>4.106367</td><td>4.197949</td><td>⋯</td><td>4.184791</td><td>4.148682</td><td>4.176648</td><td>4.387018</td><td>4.255884</td><td>3.988697</td><td>4.093996</td><td>3.994381</td><td>4.010759</td><td>4.006372</td></tr>
	<tr><th scope=row>ZXDB</th><td>4.555366</td><td>4.416483</td><td>4.766444</td><td>4.858439</td><td>4.765345</td><td>4.702525</td><td>4.626552</td><td>4.680426</td><td>4.664092</td><td>4.613109</td><td>⋯</td><td>4.707630</td><td>4.671234</td><td>4.781594</td><td>4.347300</td><td>4.914092</td><td>4.721728</td><td>4.828412</td><td>4.440002</td><td>4.970231</td><td>4.985610</td></tr>
	<tr><th scope=row>ZXDC</th><td>4.737652</td><td>4.893519</td><td>4.948034</td><td>4.998556</td><td>4.985985</td><td>4.990143</td><td>5.011917</td><td>4.990650</td><td>4.862903</td><td>4.888992</td><td>⋯</td><td>4.901736</td><td>4.956096</td><td>4.916287</td><td>5.026349</td><td>4.778690</td><td>5.399621</td><td>5.198812</td><td>4.914514</td><td>4.999590</td><td>5.016674</td></tr>
	<tr><th scope=row>ZYG11A</th><td>4.579707</td><td>4.228833</td><td>4.274655</td><td>4.274655</td><td>4.218237</td><td>4.127405</td><td>4.588900</td><td>4.372920</td><td>4.277907</td><td>4.241648</td><td>⋯</td><td>4.051505</td><td>4.147688</td><td>4.330925</td><td>4.245348</td><td>4.194728</td><td>4.108309</td><td>4.226417</td><td>4.058446</td><td>4.263681</td><td>4.359055</td></tr>
	<tr><th scope=row>ZYG11B</th><td>7.098187</td><td>7.400314</td><td>7.467173</td><td>7.231039</td><td>7.403864</td><td>7.314083</td><td>6.959457</td><td>7.147285</td><td>6.850487</td><td>6.987738</td><td>⋯</td><td>7.303608</td><td>7.565854</td><td>6.966647</td><td>6.625811</td><td>7.479652</td><td>7.466676</td><td>7.270319</td><td>7.614041</td><td>7.399101</td><td>7.327149</td></tr>
	<tr><th scope=row>ZYX</th><td>9.300109</td><td>9.028097</td><td>9.102436</td><td>9.101108</td><td>8.716208</td><td>8.582236</td><td>8.557131</td><td>9.404469</td><td>9.065257</td><td>9.228102</td><td>⋯</td><td>9.544422</td><td>8.834466</td><td>9.658579</td><td>9.739418</td><td>8.440841</td><td>9.031375</td><td>8.849947</td><td>8.172197</td><td>8.998627</td><td>9.727029</td></tr>
	<tr><th scope=row>ZZEF1</th><td>5.154425</td><td>5.354146</td><td>5.373565</td><td>5.332031</td><td>5.197398</td><td>5.411819</td><td>5.278294</td><td>5.237773</td><td>5.381028</td><td>5.656935</td><td>⋯</td><td>5.361512</td><td>5.367635</td><td>5.601214</td><td>5.745473</td><td>5.141358</td><td>5.495377</td><td>5.378354</td><td>5.249917</td><td>5.274691</td><td>5.490840</td></tr>
	<tr><th scope=row>ZZZ3</th><td>6.535700</td><td>6.738932</td><td>6.745918</td><td>6.924736</td><td>7.115012</td><td>7.075149</td><td>6.480572</td><td>6.445403</td><td>6.707779</td><td>6.449606</td><td>⋯</td><td>6.565413</td><td>6.796559</td><td>6.280222</td><td>6.359303</td><td>7.354997</td><td>6.949927</td><td>6.773739</td><td>7.111721</td><td>6.998961</td><td>6.931274</td></tr>
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
<caption>A data.frame: 29 × 4</caption>
<thead>
	<tr><th scope=col>sample</th><th scope=col>Yellow_Module</th><th scope=col>Brown_Module</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM714070</td><td>2.063423</td><td>2.729803</td><td>ADV</td></tr>
	<tr><td>GSM714071</td><td>2.080312</td><td>2.757582</td><td>ADV</td></tr>
	<tr><td>GSM714072</td><td>2.171248</td><td>2.765137</td><td>ADV</td></tr>
	<tr><td>GSM714073</td><td>2.161142</td><td>2.747168</td><td>ADV</td></tr>
	<tr><td>GSM714074</td><td>2.052517</td><td>2.594697</td><td>ADV</td></tr>
	<tr><td>GSM714075</td><td>2.001150</td><td>2.548091</td><td>ADV</td></tr>
	<tr><td>GSM714076</td><td>2.015772</td><td>2.398323</td><td>ADV</td></tr>
	<tr><td>GSM714077</td><td>2.071035</td><td>2.709875</td><td>ADV</td></tr>
	<tr><td>GSM714078</td><td>2.009346</td><td>2.647153</td><td>ADV</td></tr>
	<tr><td>GSM714079</td><td>2.123672</td><td>2.594856</td><td>ADV</td></tr>
	<tr><td>GSM714080</td><td>2.093609</td><td>2.607381</td><td>ADV</td></tr>
	<tr><td>GSM714081</td><td>2.011439</td><td>2.465251</td><td>ADV</td></tr>
	<tr><td>GSM714082</td><td>2.006742</td><td>2.423636</td><td>ADV</td></tr>
	<tr><td>GSM714083</td><td>2.045299</td><td>2.524459</td><td>ADV</td></tr>
	<tr><td>GSM714084</td><td>2.235076</td><td>2.854767</td><td>ADV</td></tr>
	<tr><td>GSM714085</td><td>2.208613</td><td>2.824144</td><td>ADV</td></tr>
	<tr><td>GSM714086</td><td>2.138156</td><td>2.556621</td><td>EAR</td></tr>
	<tr><td>GSM714087</td><td>2.109113</td><td>2.556298</td><td>EAR</td></tr>
	<tr><td>GSM714088</td><td>1.991937</td><td>2.477699</td><td>EAR</td></tr>
	<tr><td>GSM714089</td><td>2.021021</td><td>2.454346</td><td>EAR</td></tr>
	<tr><td>GSM714090</td><td>2.028091</td><td>2.477310</td><td>EAR</td></tr>
	<tr><td>GSM714091</td><td>1.904850</td><td>2.133473</td><td>EAR</td></tr>
	<tr><td>GSM714092</td><td>1.854767</td><td>2.176780</td><td>EAR</td></tr>
	<tr><td>GSM714093</td><td>1.911551</td><td>2.427933</td><td>EAR</td></tr>
	<tr><td>GSM714094</td><td>2.017781</td><td>2.500247</td><td>EAR</td></tr>
	<tr><td>GSM714095</td><td>2.113632</td><td>2.442093</td><td>EAR</td></tr>
	<tr><td>GSM714096</td><td>2.043185</td><td>2.654559</td><td>EAR</td></tr>
	<tr><td>GSM714097</td><td>2.027081</td><td>2.441265</td><td>EAR</td></tr>
	<tr><td>GSM714098</td><td>2.014746</td><td>2.409001</td><td>EAR</td></tr>
</tbody>
</table>




```R
p_yellow  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      mutate(group =group_list)  %>% 
       ggboxplot(x = 'group', y = 'Yellow_Module', ylab = 'GSVA Score',title = 'Yellow Module',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('EAR','ADV')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_yellow
```


    
![png](Step5.1_Two_Group_Compare_GSE28829_files/Step5.1_Two_Group_Compare_GSE28829_28_0.png)
    



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
         stat_compare_means(comparisons = list(c('EAR','ADV')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_brown
```


    
![png](Step5.1_Two_Group_Compare_GSE28829_files/Step5.1_Two_Group_Compare_GSE28829_30_0.png)
    



```R
library(patchwork)
```


```R
p0 <- p_yellow + p_brown + plot_layout(guides='collect') +  plot_annotation(title = 'GSE28829',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p0
ggsave(p0, file = './GSE28829_Two_Module.pdf', height = 6, width = 9, units = 'cm')
```


    
![png](Step5.1_Two_Group_Compare_GSE28829_files/Step5.1_Two_Group_Compare_GSE28829_32_0.png)
    

