```R
library(tidyverse)
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
    



```R
###analyze aldoa
```


```R
data.frame(L_C = c(1.254,0.946,0.843), L_O = c(1.403,1.245,1.937), H_C = c(2.360,2.716,1.638), H_O = c(4.322,3.623,2.454 ))   
```


<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th scope=col>L_C</th><th scope=col>L_O</th><th scope=col>H_C</th><th scope=col>H_O</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1.254</td><td>1.403</td><td>2.360</td><td>4.322</td></tr>
	<tr><td>0.946</td><td>1.245</td><td>2.716</td><td>3.623</td></tr>
	<tr><td>0.843</td><td>1.937</td><td>1.638</td><td>2.454</td></tr>
</tbody>
</table>




```R
(1.254+0.946+0.843+1.312+0.372+0.653)/6
```


0.896666666666667



```R
p1 <- tibble(value = c(1.254,0.946,0.843,1.312,0.372,0.653,1.403,1.245,1.937,1.447,1.037,1.328,
                2.36,2.716,1.638,2.137,1.389,1.383, 4.322,3.623,2.654,2.617,2.021,2.263), group = rep(c('L_C', 'L_O', 'H_C', 'H_O'),c(6,6,6,6))) %>% 
  mutate(group = factor(group, levels = c('L_C','L_O','H_C','H_O'))) %>% 
  mutate(value = value/((1.254+0.946+0.843+1.312+0.372+0.653)/6))%>% 
   ggboxplot(x = 'group', y = 'value', title = 'ALDOA',add = 'point', fill = 'group',palette = 'lancet',legend = 'right', ylab = 'Relative Expression of ALDOA', xlab = '')+
    stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('H_C', 'H_O'),c('L_O', 'H_O')),label = 'p.signif',method = "t.test") + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), axis.title.y = element_text(size=8), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p1
```


    
![png](Step11.7_Four_Group_Compare_files/Step11.7_Four_Group_Compare_4_0.png)
    



```R
(0.999+0.609+0.676+1.079+1.02+0.909)/6
```


0.882



```R
##########CREG1
p2 <- tibble(value = c(0.999,0.609,0.676,1.079,1.02,0.909,1.24,1.138,1.579,1.157,0.901,0.983,
                 0.811,0.947,0.956,0.608,0.744,0.599,1.631,1.344,1.73,1.113,2.505,2.207), group = rep(c('L_C', 'L_O', 'H_C', 'H_O'),c(6,6,6,6))) %>% 
  mutate(group = factor(group, levels = c('L_C','L_O','H_C','H_O'))) %>% 
  mutate(value = value/((0.999+0.609+0.676+1.079+1.02+0.909)/6))%>% 
   ggboxplot(x = 'group', y = 'value', title = 'CREG1', add = 'point', fill = 'group',palette = 'lancet',legend = 'right', ylab = 'Relative Expression of CREG1', xlab = '')+
    stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('H_C', 'H_O'),c('L_O', 'H_O')),label = 'p.signif',method = "t.test") + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), axis.title.y = element_text(size=8), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p2 
```


    
![png](Step11.7_Four_Group_Compare_files/Step11.7_Four_Group_Compare_6_0.png)
    



```R
##########LGMN
p3 <- tibble(value = c(1.052,1.2,0.792,1.244,0.985,0.816,2.07,2.166,1.255,1.333,2.138,1.140,
                      1.25,1.581,1.957,1.267,1.202,1.041,2.748,2.382,1.871,1.606,1.752,1.907), group = rep(c('L_C', 'L_O', 'H_C', 'H_O'),c(6,6,6,6))) %>% 
  mutate(group = factor(group, levels = c('L_C','L_O','H_C','H_O'))) %>% 
  mutate(value = value/((1.052+1.2+0.792+1.244+0.985+0.816)/6))%>% 
   ggboxplot(x = 'group', y = 'value',  title = 'LGMN', add = 'point', fill = 'group',palette = 'lancet',legend = 'right', ylab = 'Relative Expression of LGMN', xlab = '')+
    stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('H_C', 'H_O'),c('L_O', 'H_O')),label = 'p.signif',method = "t.test") + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), axis.title.y = element_text(size=8), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p3 
```


    
![png](Step11.7_Four_Group_Compare_files/Step11.7_Four_Group_Compare_7_0.png)
    



```R
##########PKM
p4 <- tibble(value = c(1.182115,1.232007,0.585878,1.47121,0.555293,1.273495,1.473712,2.162904,2.215157,1.152817,2.89905,2.937274,
                       1.286489,1.834675,1.383181,3.503636,4.771254,3.429995,2.958919,2.447312,3.163532,5.647751,4.006297,5.20181), group = rep(c('L_C', 'L_O', 'H_C', 'H_O'),c(6,6,6,6))) %>% 
  mutate(group = factor(group, levels = c('L_C','L_O','H_C','H_O'))) %>% 
  mutate(value = value/(((1.182115+1.232007+0.585878+1.47121+0.555293+1.273495)/6)))%>% 
   ggboxplot(x = 'group', y = 'value',  title = 'PKM', add = 'point', fill = 'group',palette = 'lancet',legend = 'right', ylab = 'Relative Expression of PKM', xlab = '')+
    stat_compare_means(comparisons = list(c('L_C','H_C'),c('L_C','L_O'),c('H_C', 'H_O'),c('L_O', 'H_O')),label = 'p.signif',method = "t.test") + 
           theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), axis.title.y = element_text(size=8), strip.background = element_blank(),strip.text = element_text(size=6),plot.title=element_text(hjust=0.5, size = 8)) 
p4 
```


    
![png](Step11.7_Four_Group_Compare_files/Step11.7_Four_Group_Compare_8_0.png)
    



```R
library(patchwork)
```


```R
p <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2,guides = 'collect')
p
ggsave(p, file = './4_Gene_QPCR.pdf', height = 9, width = 11,units = 'cm')
```


    
![png](Step11.7_Four_Group_Compare_files/Step11.7_Four_Group_Compare_10_0.png)
    



```R

```
