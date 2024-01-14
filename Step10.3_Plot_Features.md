```R
rm(list=ls())
options(stringsAsFactors = F)
```


```R
library(Seurat)
library(tidyverse)
library(patchwork)
```

    Attaching SeuratObject
    
    Attaching sp
    
    Warning message in system("timedatectl", intern = TRUE):
    “running command 'timedatectl' had status 1”
    ── [1mAttaching packages[22m ─────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
    
    [32m✔[39m [34mggplot2[39m 3.4.0     [32m✔[39m [34mpurrr  [39m 0.3.4
    [32m✔[39m [34mtibble [39m 3.1.7     [32m✔[39m [34mdplyr  [39m 1.0.7
    [32m✔[39m [34mtidyr  [39m 1.1.4     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 2.1.1     [32m✔[39m [34mforcats[39m 0.5.1
    
    ── [1mConflicts[22m ────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    



```R
immune.combined <- readRDS('~/AS/AS_Mouse/AS_Mouse2/Final_result/Figure8_Human_Plaques/Step1_AddDoublet_Seurat.Rds')
```


```R
DimPlot(immune.combined, group.by ='celltype') 
```


    
![png](Step10.3_Plot_Features_files/Step10.3_Plot_Features_3_0.png)
    



```R
source('~/Resources/Functions/Plot_colorPaletters.R')
```


```R
Palettes
```


<dl>
	<dt>$group_pal</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#D51F26'</li><li>'#272E6A'</li><li>'#208A42'</li><li>'#89288F'</li><li>'#6387C5'</li></ol>
</dd>
	<dt>$group_pal2</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#E64B35'</li><li>'#4DBBD5'</li><li>'#00A087'</li><li>'#3C5488'</li><li>'#F39B7F'</li></ol>
</dd>
	<dt>$seurat_palletes_group</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#DA4C35'</li><li>'#4EB1C9'</li><li>'#169982'</li><li>'#3B5282'</li><li>'#EB977D'</li><li>'#818DAF'</li><li>'#7C6048'</li><li>'#8ECABC'</li><li>'#AF9B85'</li><li>'#4456A5'</li><li>'#C43E32'</li><li>'#729757'</li><li>'#E9E085'</li><li>'#466780'</li><li>'#B46139'</li><li>'#7B2063'</li><li>'#73BE69'</li><li>'#7B171C'</li><li>'#D5D5CD'</li><li>'#F6B04A'</li><li>'#AAAE7C'</li><li>'#5A8AA2'</li><li>'#CE9163'</li><li>'#AC726D'</li><li>'#898B79'</li><li>'#715562'</li><li>'#ED6E1C'</li><li>'#BB1B21'</li><li>'#823F8F'</li><li>'#589195'</li><li>'#F090A2'</li><li>'#000000'</li></ol>
</dd>
	<dt>$seurat_palletes</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#DA4C35'</dd><dt>2</dt><dd>'#4EB1C9'</dd><dt>3</dt><dd>'#169982'</dd><dt>4</dt><dd>'#3B5282'</dd><dt>5</dt><dd>'#EB977D'</dd><dt>6</dt><dd>'#818DAF'</dd><dt>7</dt><dd>'#7C6048'</dd><dt>8</dt><dd>'#8ECABC'</dd><dt>9</dt><dd>'#AF9B85'</dd><dt>10</dt><dd>'#4456A5'</dd><dt>11</dt><dd>'#C43E32'</dd><dt>12</dt><dd>'#729757'</dd><dt>13</dt><dd>'#E9E085'</dd><dt>14</dt><dd>'#466780'</dd><dt>15</dt><dd>'#B46139'</dd><dt>16</dt><dd>'#7B2063'</dd><dt>17</dt><dd>'#73BE69'</dd><dt>18</dt><dd>'#7B171C'</dd><dt>19</dt><dd>'#D5D5CD'</dd><dt>20</dt><dd>'#F6B04A'</dd><dt>21</dt><dd>'#AAAE7C'</dd><dt>22</dt><dd>'#5A8AA2'</dd><dt>23</dt><dd>'#CE9163'</dd><dt>24</dt><dd>'#AC726D'</dd><dt>25</dt><dd>'#898B79'</dd><dt>26</dt><dd>'#715562'</dd><dt>27</dt><dd>'#ED6E1C'</dd><dt>28</dt><dd>'#BB1B21'</dd><dt>29</dt><dd>'#823F8F'</dd><dt>30</dt><dd>'#589195'</dd><dt>31</dt><dd>'#F090A2'</dd><dt>32</dt><dd>'#000000'</dd></dl>
</dd>
	<dt>$group_stallion</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#F13A13'</dd><dt>2</dt><dd>'#F47D2B'</dd><dt>3</dt><dd>'#FEE500'</dd><dt>4</dt><dd>'#83AD00'</dd><dt>5</dt><dd>'#208A42'</dd><dt>6</dt><dd>'#00C0BA'</dd><dt>7</dt><dd>'#7B6FD0'</dd><dt>8</dt><dd>'#C06CAB'</dd><dt>19</dt><dd>'#F37B7D'</dd><dt>10</dt><dd>'#e0598b'</dd></dl>
</dd>
	<dt>$stallion</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#EF7F48'</dd><dt>2</dt><dd>'#D69100'</dd><dt>3</dt><dd>'#C69900'</dd><dt>4</dt><dd>''</dd><dt>5</dt><dd>'#00BE6B'</dd><dt>6</dt><dd>'#00C0BA'</dd><dt>7</dt><dd>'#00B9E1'</dd><dt>8</dt><dd>'#00B3F1'</dd><dt>19</dt><dd>'#E6C2DC'</dd><dt>10</dt><dd>'#8794FF'</dd><dt>11</dt><dd>'#B086FF'</dd><dt>12</dt><dd>'#E46DF6'</dd><dt>13</dt><dd>'#F365E6'</dd><dt>14</dt><dd>'#D24B27'</dd><dt>15</dt><dd>'#3BBCA8'</dd><dt>16</dt><dd>'#6E4B9E'</dd><dt>17</dt><dd>'#0C727C'</dd><dt>18</dt><dd>'#7E1416'</dd><dt>9</dt><dd>'#00ABFD'</dd><dt>20</dt><dd>'#3D3D3D'</dd></dl>
</dd>
	<dt>$stallion2</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#EF7F48'</dd><dt>2</dt><dd>'#D69100'</dd><dt>3</dt><dd>'#C69900'</dd><dt>4</dt><dd>'#83AD00'</dd><dt>5</dt><dd>'#00BE6B'</dd><dt>6</dt><dd>'#FEE500'</dd><dt>7</dt><dd>'#8A9FD1'</dd><dt>8</dt><dd>'#C06CAB'</dd><dt>19</dt><dd>'#E6C2DC'</dd><dt>10</dt><dd>'#90D5E4'</dd><dt>11</dt><dd>'#89C75F'</dd><dt>12</dt><dd>'#F37B7D'</dd><dt>13</dt><dd>'#9983BD'</dd><dt>14</dt><dd>'#D24B27'</dd><dt>15</dt><dd>'#3BBCA8'</dd><dt>16</dt><dd>'#6E4B9E'</dd><dt>17</dt><dd>'#0C727C'</dd><dt>18</dt><dd>'#7E1416'</dd><dt>9</dt><dd>'#D8A767'</dd></dl>
</dd>
	<dt>$calm</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#7DD06F'</dd><dt>2</dt><dd>'#844081'</dd><dt>3</dt><dd>'#688EC1'</dd><dt>4</dt><dd>'#C17E73'</dd><dt>5</dt><dd>'#484125'</dd><dt>6</dt><dd>'#6CD3A7'</dd><dt>7</dt><dd>'#597873'</dd><dt>8</dt><dd>'#7B6FD0'</dd><dt>9</dt><dd>'#CF4A31'</dd><dt>10</dt><dd>'#D0CD47'</dd><dt>11</dt><dd>'#722A2D'</dd><dt>12</dt><dd>'#CBC594'</dd><dt>13</dt><dd>'#D19EC4'</dd><dt>14</dt><dd>'#5A7E36'</dd><dt>15</dt><dd>'#D4477D'</dd><dt>16</dt><dd>'#403552'</dd><dt>17</dt><dd>'#76D73C'</dd><dt>18</dt><dd>'#96CED5'</dd><dt>19</dt><dd>'#CE54D1'</dd><dt>20</dt><dd>'#C48736'</dd></dl>
</dd>
	<dt>$kelly</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#FFB300'</dd><dt>2</dt><dd>'#803E75'</dd><dt>3</dt><dd>'#FF6800'</dd><dt>4</dt><dd>'#A6BDD7'</dd><dt>5</dt><dd>'#C10020'</dd><dt>6</dt><dd>'#CEA262'</dd><dt>7</dt><dd>'#817066'</dd><dt>8</dt><dd>'#007D34'</dd><dt>9</dt><dd>'#F6768E'</dd><dt>10</dt><dd>'#00538A'</dd><dt>11</dt><dd>'#FF7A5C'</dd><dt>12</dt><dd>'#53377A'</dd><dt>13</dt><dd>'#FF8E00'</dd><dt>14</dt><dd>'#B32851'</dd><dt>15</dt><dd>'#F4C800'</dd><dt>16</dt><dd>'#7F180D'</dd><dt>17</dt><dd>'#93AA00'</dd><dt>18</dt><dd>'#593315'</dd><dt>19</dt><dd>'#F13A13'</dd><dt>20</dt><dd>'#232C16'</dd></dl>
</dd>
	<dt>$bear</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#faa818'</dd><dt>2</dt><dd>'#41a30d'</dd><dt>3</dt><dd>'#fbdf72'</dd><dt>4</dt><dd>'#367d7d'</dd><dt>5</dt><dd>'#d33502'</dd><dt>6</dt><dd>'#6ebcbc'</dd><dt>7</dt><dd>'#37526d'</dd><dt>8</dt><dd>'#916848'</dd><dt>9</dt><dd>'#f5b390'</dd><dt>10</dt><dd>'#342739'</dd><dt>11</dt><dd>'#bed678'</dd><dt>12</dt><dd>'#a6d9ee'</dd><dt>13</dt><dd>'#0d74b6'</dd><dt>14</dt><dd>'#60824f'</dd><dt>15</dt><dd>'#725ca5'</dd><dt>16</dt><dd>'#e0598b'</dd></dl>
</dd>
	<dt>$ironMan</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>9</dt><dd>'#371377'</dd><dt>3</dt><dd>'#7700FF'</dd><dt>2</dt><dd>'#9E0142'</dd><dt>10</dt><dd>'#FF0080'</dd><dt>14</dt><dd>'#DC494C'</dd><dt>12</dt><dd>'#F88D51'</dd><dt>1</dt><dd>'#FAD510'</dd><dt>8</dt><dd>'#FFFF5F'</dd><dt>4</dt><dd>'#88CFA4'</dd><dt>13</dt><dd>'#238B45'</dd><dt>5</dt><dd>'#02401B'</dd><dt>7</dt><dd>'#0AD7D3'</dd><dt>11</dt><dd>'#046C9A'</dd><dt>6</dt><dd>'#A2A475'</dd><dt>15</dt><dd>'grey35'</dd></dl>
</dd>
	<dt>$circus</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#D52126'</dd><dt>2</dt><dd>'#88CCEE'</dd><dt>3</dt><dd>'#FEE52C'</dd><dt>4</dt><dd>'#117733'</dd><dt>5</dt><dd>'#CC61B0'</dd><dt>6</dt><dd>'#99C945'</dd><dt>7</dt><dd>'#2F8AC4'</dd><dt>8</dt><dd>'#332288'</dd><dt>9</dt><dd>'#E68316'</dd><dt>10</dt><dd>'#661101'</dd><dt>11</dt><dd>'#F97B72'</dd><dt>12</dt><dd>'#DDCC77'</dd><dt>13</dt><dd>'#11A579'</dd><dt>14</dt><dd>'#89288F'</dd><dt>15</dt><dd>'#E73F74'</dd></dl>
</dd>
	<dt>$paired</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>9</dt><dd>'#A6CDE2'</dd><dt>1</dt><dd>'#1E78B4'</dd><dt>3</dt><dd>'#74C476'</dd><dt>12</dt><dd>'#34A047'</dd><dt>11</dt><dd>'#F59899'</dd><dt>2</dt><dd>'#E11E26'</dd><dt>10</dt><dd>'#FCBF6E'</dd><dt>4</dt><dd>'#F47E1F'</dd><dt>5</dt><dd>'#CAB2D6'</dd><dt>8</dt><dd>'#6A3E98'</dd><dt>6</dt><dd>'#FAF39B'</dd><dt>7</dt><dd>'#B15928'</dd></dl>
</dd>
	<dt>$grove</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>11</dt><dd>'#1a1334'</dd><dt>9</dt><dd>'#01545a'</dd><dt>1</dt><dd>'#017351'</dd><dt>6</dt><dd>'#03c383'</dd><dt>8</dt><dd>'#aad962'</dd><dt>2</dt><dd>'#fbbf45'</dd><dt>10</dt><dd>'#ef6a32'</dd><dt>3</dt><dd>'#ed0345'</dd><dt>7</dt><dd>'#a12a5e'</dd><dt>5</dt><dd>'#710162'</dd><dt>4</dt><dd>'#3B9AB2'</dd></dl>
</dd>
	<dt>$summerNight</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#2a7185'</dd><dt>2</dt><dd>'#a64027'</dd><dt>3</dt><dd>'#fbdf72'</dd><dt>4</dt><dd>'#60824f'</dd><dt>5</dt><dd>'#9cdff0'</dd><dt>6</dt><dd>'#022336'</dd><dt>7</dt><dd>'#725ca5'</dd></dl>
</dd>
	<dt>$zissou</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#3B9AB2'</dd><dt>4</dt><dd>'#78B7C5'</dd><dt>3</dt><dd>'#EBCC2A'</dd><dt>5</dt><dd>'#E1AF00'</dd><dt>2</dt><dd>'#F21A00'</dd></dl>
</dd>
	<dt>$darjeeling</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#FF0000'</dd><dt>2</dt><dd>'#00A08A'</dd><dt>3</dt><dd>'#F2AD00'</dd><dt>4</dt><dd>'#F98400'</dd><dt>5</dt><dd>'#5BBCD6'</dd></dl>
</dd>
	<dt>$rushmore</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#E1BD6D'</dd><dt>5</dt><dd>'#EABE94'</dd><dt>2</dt><dd>'#0B775E'</dd><dt>4</dt><dd>'#35274A'</dd><dt>3</dt><dd>'#F2300F'</dd></dl>
</dd>
	<dt>$captain</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'grey'</dd><dt>2</dt><dd>'#A1CDE1'</dd><dt>3</dt><dd>'#12477C'</dd><dt>4</dt><dd>'#EC9274'</dd><dt>5</dt><dd>'#67001E'</dd></dl>
</dd>
	<dt>$horizon</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#000075'</dd><dt>4</dt><dd>'#2E00FF'</dd><dt>6</dt><dd>'#9408F7'</dd><dt>10</dt><dd>'#C729D6'</dd><dt>8</dt><dd>'#FA4AB5'</dd><dt>3</dt><dd>'#FF6A95'</dd><dt>7</dt><dd>'#FF8B74'</dd><dt>5</dt><dd>'#FFAC53'</dd><dt>9</dt><dd>'#FFCD32'</dd><dt>2</dt><dd>'#FFFF60'</dd></dl>
</dd>
	<dt>$horizonExtra</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#000436'</dd><dt>4</dt><dd>'#021EA9'</dd><dt>6</dt><dd>'#1632FB'</dd><dt>8</dt><dd>'#6E34FC'</dd><dt>3</dt><dd>'#C732D5'</dd><dt>9</dt><dd>'#FD619D'</dd><dt>7</dt><dd>'#FF9965'</dd><dt>5</dt><dd>'#FFD32B'</dd><dt>2</dt><dd>'#FFFC5A'</dd></dl>
</dd>
	<dt>$blueYellow</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#352A86'</dd><dt>2</dt><dd>'#343DAE'</dd><dt>3</dt><dd>'#0262E0'</dd><dt>4</dt><dd>'#1389D2'</dd><dt>5</dt><dd>'#2DB7A3'</dd><dt>6</dt><dd>'#A5BE6A'</dd><dt>7</dt><dd>'#F8BA43'</dd><dt>8</dt><dd>'#F6DA23'</dd><dt>9</dt><dd>'#F8FA0D'</dd></dl>
</dd>
	<dt>$sambaNight</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>6</dt><dd>'#1873CC'</dd><dt>2</dt><dd>'#1798E5'</dd><dt>8</dt><dd>'#00BFFF'</dd><dt>5</dt><dd>'#4AC596'</dd><dt>1</dt><dd>'#00CC00'</dd><dt>4</dt><dd>'#A2E700'</dd><dt>9</dt><dd>'#FFFF00'</dd><dt>7</dt><dd>'#FFD200'</dd><dt>3</dt><dd>'#FFA500'</dd></dl>
</dd>
	<dt>$solarExtra</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>5</dt><dd>'#3361A5'</dd><dt>7</dt><dd>'#248AF3'</dd><dt>1</dt><dd>'#14B3FF'</dd><dt>8</dt><dd>'#88CEEF'</dd><dt>9</dt><dd>'#C1D5DC'</dd><dt>4</dt><dd>'#EAD397'</dd><dt>3</dt><dd>'#FDB31A'</dd><dt>2</dt><dd>'#E42A2A'</dd><dt>6</dt><dd>'#A31D1D'</dd></dl>
</dd>
	<dt>$whitePurple</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>9</dt><dd>'#f7fcfd'</dd><dt>6</dt><dd>'#e0ecf4'</dd><dt>8</dt><dd>'#bfd3e6'</dd><dt>5</dt><dd>'#9ebcda'</dd><dt>2</dt><dd>'#8c96c6'</dd><dt>4</dt><dd>'#8c6bb1'</dd><dt>7</dt><dd>'#88419d'</dd><dt>3</dt><dd>'#810f7c'</dd><dt>1</dt><dd>'#4d004b'</dd></dl>
</dd>
	<dt>$whiteBlue</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>9</dt><dd>'#fff7fb'</dd><dt>6</dt><dd>'#ece7f2'</dd><dt>8</dt><dd>'#d0d1e6'</dd><dt>5</dt><dd>'#a6bddb'</dd><dt>2</dt><dd>'#74a9cf'</dd><dt>4</dt><dd>'#3690c0'</dd><dt>7</dt><dd>'#0570b0'</dd><dt>3</dt><dd>'#045a8d'</dd><dt>1</dt><dd>'#023858'</dd></dl>
</dd>
	<dt>$whiteRed</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'white'</dd><dt>2</dt><dd>'red'</dd></dl>
</dd>
	<dt>$comet</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#E6E7E8'</dd><dt>2</dt><dd>'#3A97FF'</dd><dt>3</dt><dd>'#8816A7'</dd><dt>4</dt><dd>'black'</dd></dl>
</dd>
	<dt>$greenBlue</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>4</dt><dd>'#e0f3db'</dd><dt>7</dt><dd>'#ccebc5'</dd><dt>2</dt><dd>'#a8ddb5'</dd><dt>5</dt><dd>'#4eb3d3'</dd><dt>3</dt><dd>'#2b8cbe'</dd><dt>6</dt><dd>'#0868ac'</dd><dt>1</dt><dd>'#084081'</dd></dl>
</dd>
	<dt>$beach</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>4</dt><dd>'#87D2DB'</dd><dt>1</dt><dd>'#5BB1CB'</dd><dt>6</dt><dd>'#4F66AF'</dd><dt>3</dt><dd>'#F15F30'</dd><dt>5</dt><dd>'#F7962E'</dd><dt>2</dt><dd>'#FCEE2B'</dd></dl>
</dd>
	<dt>$coolwarm</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#4858A7'</dd><dt>4</dt><dd>'#788FC8'</dd><dt>5</dt><dd>'#D6DAE1'</dd><dt>3</dt><dd>'#F49B7C'</dd><dt>2</dt><dd>'#B51F29'</dd></dl>
</dd>
	<dt>$fireworks</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>5</dt><dd>'white'</dd><dt>2</dt><dd>'#2488F0'</dd><dt>4</dt><dd>'#7F3F98'</dd><dt>3</dt><dd>'#E22929'</dd><dt>1</dt><dd>'#FCB31A'</dd></dl>
</dd>
	<dt>$greyMagma</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>2</dt><dd>'grey'</dd><dt>4</dt><dd>'#FB8861FF'</dd><dt>5</dt><dd>'#B63679FF'</dd><dt>3</dt><dd>'#51127CFF'</dd><dt>1</dt><dd>'#000004FF'</dd></dl>
</dd>
	<dt>$fireworks2</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>5</dt><dd>'black'</dd><dt>2</dt><dd>'#2488F0'</dd><dt>4</dt><dd>'#7F3F98'</dd><dt>3</dt><dd>'#E22929'</dd><dt>1</dt><dd>'#FCB31A'</dd></dl>
</dd>
	<dt>$purpleOrange</dt>
		<dd><style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>5</dt><dd>'#581845'</dd><dt>2</dt><dd>'#900C3F'</dd><dt>4</dt><dd>'#C70039'</dd><dt>3</dt><dd>'#FF5744'</dd><dt>1</dt><dd>'#FFC30F'</dd></dl>
</dd>
</dl>




```R
immune.combined
```


    An object of class Seurat 
    25621 features across 48624 samples within 2 assays 
    Active assay: RNA (22621 features, 0 variable features)
     1 other assay present: integrated
     2 dimensional reductions calculated: pca, umap



```R
 DimPlot(immune.combined, group.by ='sample',cols =Palettes[['calm']]  %>% unname() )  + theme_void() + 
NoLegend()  + labs(title = '')
```


    
![png](Step10.3_Plot_Features_files/Step10.3_Plot_Features_7_0.png)
    



```R
pdf('~/AS_HG/Final_Results/Figure4_New/Seurat_Sample.pdf')
DimPlot(immune.combined, label = T, cols =Palettes[['calm']]  %>% unname()  %>% .[c(1:12)], pt.size = 0.1,group.by = 'sample') + labs (title = '48624 cells') +  theme_void() + theme(plot.title=element_text(hjust=0.5)) + NoLegend()
DimPlot(immune.combined, label = F, cols =Palettes[['calm']]  %>% unname()  %>% .[c(1:12)],, pt.size = 0.1,group.by = 'sample')  + labs (title = '48624 cells') +  theme_void() + theme(plot.title=element_text(hjust=0.5)) + NoLegend()
DimPlot(immune.combined, label = F, cols =Palettes[['calm']]  %>% unname()  %>% .[c(1:12)],, pt.size = 0.1,group.by = 'sample')  + labs (title = '48624 cells') +  theme_void() + theme(plot.title=element_text(hjust=0.5)) 
dev.off()
```


<strong>png:</strong> 2



```R
tiff('~/AS_HG/Final_Results/Figure4_New/Seurat_Sample.tiff')
DimPlot(immune.combined, label = F, cols =Palettes[['calm']]  %>% unname()  %>% .[c(1:12)],, pt.size = 0.1,group.by = 'sample')  + labs (title = '') +  theme_void() + theme(plot.title=element_text(hjust=0.5)) + NoLegend()
dev.off()
```


<strong>png:</strong> 2



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

```


```R
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(yellow_top100),
                 ctrl = length(yellow_top100),
                 name = 'Yellow_Module',
                 search = FALSE,
                 assay = "integrated")


```

    Warning message:
    “The following features are not present in the object: STXBP5, SLC38A2, ARHGEF6, CCDC93, ABCA7, FNBP4, TNFRSF26, HEXIM1, SLC5A3, SLC6A12, AMPD2, SLC7A6, PROKR1, KLHL21, ZNRF1, BCL6, NUMB, HIPK1, CELF1, TTLL4, LDLRAD3, RAF1, YY1, DPP7, ARRB1, RNASEH2C, C77080, ARMT1, ZFP628, SUMO3, ZDHHC14, LENG8, PPP1R2, GM10521, IRAK1, CEP78, PBXIP1, FAM193B, MUS_MUSCULUS_NEWGENE_530, MACROD1, ATG16L2, PHOSPHO2, DNMT1, MAST2, P2RY2, ZFP503, MUS_MUSCULUS_NEWGENE_707, GTF2IRD2, TRIM3, MSH3, NKTR, NFAT5, PGS1, ZFC3H1, GSTP1, ADCY7, KANSL2, GMIP, HERC1, ANKRD13B, PRPF39, GALNT6, HAUS4, PLCB3, PLD2, UBAP1, TGFBR1, DCLRE1B, SETD1A, MKRN2, MBP, TMEM63A, CCNL2, KLF1, ZFP335, RNASEH2A, WDR26, SNX27, WDFY1, MUS_MUSCULUS_NEWGENE_1249, LUC7L3, 4833420G17RIK, SLC26A11, LARP1B, RNF220, TRIM28, not searching for symbol synonyms”



    Error in .local(x, na.rm, dims, ...): invalid object passed to as_cholmod_sparse
    Traceback:


    1. AddModuleScore(object = immune.combined, features = list(yellow_top100), 
     .     ctrl = length(yellow_top100), name = "Yellow_Module", search = FALSE, 
     .     assay = "integrated")

    2. Matrix::rowMeans(x = assay.data[pool, ])

    3. Matrix::rowMeans(x = assay.data[pool, ])

    4. .local(x, na.rm, dims, ...)



```R
immune.combined <- AddModuleScore(
                 object = immune.combined,
                 features = list(brown_top100),
                 ctrl = length(brown_top100),
                 name = 'Brown_Module',
                 search = FALSE,
                 assay = "integrated")
```
