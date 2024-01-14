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
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*6)
```


```R
gset  <-  GEOquery::getGEO('GSE100927',getGPL = F)
```

    Found 1 file(s)
    
    GSE100927_series_matrix.txt.gz
    
    [1mRows: [22m[34m47264[39m [1mColumns: [22m[34m105[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m "\t"
    [31mchr[39m   (1): ID_REF
    [32mdbl[39m (104): GSM2696609, GSM2696610, GSM2696611, GSM2696612, GSM2696613, GSM26...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
normalized_gset  <- gset$GSE100927_series_matrix.txt.gz@assayData$exprs
```


```R
range(normalized_gset)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>2.28</li><li>14.1</li></ol>




```R
normalized_gset
```


<table class="dataframe">
<caption>A matrix: 47264 × 104 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>GSM2696609</th><th scope=col>GSM2696610</th><th scope=col>GSM2696611</th><th scope=col>GSM2696612</th><th scope=col>GSM2696613</th><th scope=col>GSM2696614</th><th scope=col>GSM2696615</th><th scope=col>GSM2696616</th><th scope=col>GSM2696617</th><th scope=col>GSM2696618</th><th scope=col>⋯</th><th scope=col>GSM2696703</th><th scope=col>GSM2696704</th><th scope=col>GSM2696705</th><th scope=col>GSM2696706</th><th scope=col>GSM2696707</th><th scope=col>GSM2696708</th><th scope=col>GSM2696709</th><th scope=col>GSM2696710</th><th scope=col>GSM2696711</th><th scope=col>GSM2696712</th></tr>
</thead>
<tbody>
	<tr><th scope=row>A_19_P00315452</th><td>4.45</td><td>4.47</td><td>5.16</td><td>4.96</td><td>4.76</td><td>4.93</td><td>5.03</td><td>5.54</td><td>4.50</td><td>4.98</td><td>⋯</td><td>5.01</td><td>4.81</td><td>4.82</td><td>4.70</td><td>4.55</td><td>6.17</td><td>4.57</td><td>4.65</td><td>4.62</td><td>4.69</td></tr>
	<tr><th scope=row>A_19_P00315459</th><td>4.99</td><td>5.12</td><td>4.56</td><td>5.02</td><td>4.88</td><td>5.13</td><td>5.03</td><td>4.92</td><td>4.47</td><td>5.24</td><td>⋯</td><td>5.23</td><td>4.75</td><td>4.79</td><td>5.14</td><td>4.86</td><td>4.88</td><td>5.30</td><td>5.00</td><td>4.68</td><td>4.86</td></tr>
	<tr><th scope=row>A_19_P00315482</th><td>4.00</td><td>4.06</td><td>3.92</td><td>3.93</td><td>3.96</td><td>3.95</td><td>3.88</td><td>3.87</td><td>3.83</td><td>3.87</td><td>⋯</td><td>4.09</td><td>4.05</td><td>3.82</td><td>4.03</td><td>3.90</td><td>4.08</td><td>3.70</td><td>3.89</td><td>3.92</td><td>3.77</td></tr>
	<tr><th scope=row>A_19_P00315492</th><td>4.13</td><td>4.41</td><td>4.15</td><td>3.83</td><td>4.11</td><td>3.99</td><td>3.91</td><td>4.09</td><td>3.91</td><td>3.74</td><td>⋯</td><td>4.08</td><td>4.27</td><td>3.83</td><td>3.94</td><td>3.91</td><td>3.76</td><td>3.96</td><td>4.11</td><td>3.99</td><td>3.83</td></tr>
	<tr><th scope=row>A_19_P00315493</th><td>6.30</td><td>5.21</td><td>5.10</td><td>5.35</td><td>6.08</td><td>5.82</td><td>6.50</td><td>4.89</td><td>5.22</td><td>5.26</td><td>⋯</td><td>6.11</td><td>5.57</td><td>5.39</td><td>5.63</td><td>5.21</td><td>5.36</td><td>5.38</td><td>5.44</td><td>6.40</td><td>5.57</td></tr>
	<tr><th scope=row>A_19_P00315502</th><td>3.91</td><td>3.61</td><td>4.24</td><td>3.76</td><td>3.84</td><td>3.89</td><td>3.77</td><td>3.79</td><td>3.82</td><td>3.70</td><td>⋯</td><td>3.75</td><td>3.87</td><td>3.94</td><td>3.68</td><td>3.73</td><td>3.72</td><td>3.92</td><td>3.82</td><td>3.82</td><td>3.83</td></tr>
	<tr><th scope=row>A_19_P00315506</th><td>6.38</td><td>6.33</td><td>6.70</td><td>6.26</td><td>6.52</td><td>6.58</td><td>6.69</td><td>6.26</td><td>7.22</td><td>6.81</td><td>⋯</td><td>6.25</td><td>6.56</td><td>6.82</td><td>6.70</td><td>7.04</td><td>6.32</td><td>6.66</td><td>7.12</td><td>6.33</td><td>6.81</td></tr>
	<tr><th scope=row>A_19_P00315518</th><td>3.93</td><td>3.48</td><td>3.75</td><td>3.67</td><td>3.84</td><td>3.82</td><td>3.82</td><td>3.88</td><td>3.65</td><td>3.78</td><td>⋯</td><td>4.02</td><td>3.74</td><td>3.66</td><td>3.87</td><td>3.86</td><td>3.84</td><td>3.75</td><td>3.85</td><td>3.81</td><td>3.81</td></tr>
	<tr><th scope=row>A_19_P00315519</th><td>3.40</td><td>3.57</td><td>3.66</td><td>3.58</td><td>3.71</td><td>3.78</td><td>3.78</td><td>3.84</td><td>3.59</td><td>3.69</td><td>⋯</td><td>3.74</td><td>3.63</td><td>3.38</td><td>3.68</td><td>3.67</td><td>3.51</td><td>3.67</td><td>3.82</td><td>3.74</td><td>3.60</td></tr>
	<tr><th scope=row>A_19_P00315524</th><td>5.25</td><td>5.71</td><td>4.60</td><td>5.15</td><td>5.78</td><td>6.43</td><td>6.78</td><td>4.88</td><td>4.43</td><td>5.35</td><td>⋯</td><td>6.58</td><td>5.36</td><td>5.68</td><td>6.06</td><td>6.61</td><td>6.53</td><td>4.41</td><td>5.32</td><td>6.58</td><td>6.37</td></tr>
	<tr><th scope=row>A_19_P00315528</th><td>3.94</td><td>4.41</td><td>4.09</td><td>4.03</td><td>4.08</td><td>4.03</td><td>4.33</td><td>3.98</td><td>3.93</td><td>4.36</td><td>⋯</td><td>4.26</td><td>3.83</td><td>4.02</td><td>4.07</td><td>4.15</td><td>4.39</td><td>3.91</td><td>4.03</td><td>4.21</td><td>4.19</td></tr>
	<tr><th scope=row>A_19_P00315529</th><td>3.86</td><td>4.44</td><td>3.85</td><td>4.10</td><td>4.19</td><td>4.05</td><td>4.41</td><td>4.30</td><td>4.01</td><td>4.20</td><td>⋯</td><td>4.08</td><td>4.06</td><td>3.91</td><td>4.22</td><td>4.10</td><td>4.36</td><td>4.31</td><td>4.38</td><td>4.17</td><td>4.06</td></tr>
	<tr><th scope=row>A_19_P00315538</th><td>3.69</td><td>4.18</td><td>3.76</td><td>3.51</td><td>3.83</td><td>4.46</td><td>3.87</td><td>3.83</td><td>3.91</td><td>3.72</td><td>⋯</td><td>4.12</td><td>3.82</td><td>3.93</td><td>3.61</td><td>3.72</td><td>3.72</td><td>3.81</td><td>3.61</td><td>3.88</td><td>3.76</td></tr>
	<tr><th scope=row>A_19_P00315541</th><td>3.73</td><td>3.77</td><td>3.80</td><td>3.88</td><td>3.53</td><td>3.75</td><td>3.61</td><td>3.58</td><td>3.79</td><td>3.83</td><td>⋯</td><td>3.91</td><td>4.03</td><td>3.80</td><td>3.78</td><td>3.71</td><td>3.69</td><td>3.55</td><td>3.69</td><td>3.76</td><td>3.83</td></tr>
	<tr><th scope=row>A_19_P00315543</th><td>4.78</td><td>4.42</td><td>5.36</td><td>4.76</td><td>4.88</td><td>4.29</td><td>4.20</td><td>5.43</td><td>4.57</td><td>4.43</td><td>⋯</td><td>4.99</td><td>4.49</td><td>4.69</td><td>4.91</td><td>4.85</td><td>4.65</td><td>4.62</td><td>4.73</td><td>4.79</td><td>4.80</td></tr>
	<tr><th scope=row>A_19_P00315550</th><td>6.81</td><td>6.45</td><td>6.44</td><td>5.41</td><td>6.27</td><td>5.79</td><td>5.59</td><td>6.74</td><td>5.25</td><td>5.72</td><td>⋯</td><td>6.40</td><td>6.53</td><td>6.40</td><td>6.82</td><td>6.57</td><td>5.58</td><td>5.06</td><td>6.50</td><td>5.75</td><td>6.44</td></tr>
	<tr><th scope=row>A_19_P00315551</th><td>6.88</td><td>6.50</td><td>6.08</td><td>5.93</td><td>6.38</td><td>6.08</td><td>5.80</td><td>6.40</td><td>5.58</td><td>5.86</td><td>⋯</td><td>6.10</td><td>6.41</td><td>6.40</td><td>6.76</td><td>6.34</td><td>5.99</td><td>5.31</td><td>6.58</td><td>6.03</td><td>6.42</td></tr>
	<tr><th scope=row>A_19_P00315554</th><td>4.06</td><td>3.94</td><td>3.77</td><td>3.68</td><td>3.92</td><td>5.27</td><td>4.14</td><td>3.77</td><td>3.93</td><td>3.85</td><td>⋯</td><td>4.15</td><td>3.86</td><td>4.08</td><td>3.99</td><td>4.33</td><td>4.29</td><td>3.85</td><td>3.85</td><td>4.11</td><td>4.20</td></tr>
	<tr><th scope=row>A_19_P00315581</th><td>7.84</td><td>7.93</td><td>7.87</td><td>7.92</td><td>7.87</td><td>7.89</td><td>8.01</td><td>7.19</td><td>7.58</td><td>7.71</td><td>⋯</td><td>7.78</td><td>8.19</td><td>7.95</td><td>7.99</td><td>8.02</td><td>7.84</td><td>7.51</td><td>7.93</td><td>7.96</td><td>7.85</td></tr>
	<tr><th scope=row>A_19_P00315583</th><td>4.77</td><td>4.94</td><td>4.05</td><td>4.80</td><td>4.88</td><td>5.15</td><td>4.58</td><td>4.64</td><td>4.21</td><td>4.55</td><td>⋯</td><td>4.81</td><td>4.87</td><td>4.94</td><td>4.87</td><td>4.73</td><td>4.70</td><td>4.84</td><td>4.59</td><td>4.78</td><td>4.87</td></tr>
	<tr><th scope=row>A_19_P00315584</th><td>5.44</td><td>5.50</td><td>4.68</td><td>5.33</td><td>5.48</td><td>5.82</td><td>5.28</td><td>5.29</td><td>4.55</td><td>5.25</td><td>⋯</td><td>5.43</td><td>5.58</td><td>5.63</td><td>5.50</td><td>5.31</td><td>5.32</td><td>5.28</td><td>5.07</td><td>5.53</td><td>5.60</td></tr>
	<tr><th scope=row>A_19_P00315587</th><td>3.70</td><td>3.71</td><td>3.91</td><td>3.85</td><td>4.13</td><td>3.90</td><td>3.68</td><td>3.78</td><td>3.83</td><td>3.86</td><td>⋯</td><td>3.77</td><td>3.88</td><td>3.42</td><td>3.82</td><td>3.77</td><td>3.79</td><td>3.67</td><td>3.75</td><td>3.86</td><td>4.05</td></tr>
	<tr><th scope=row>A_19_P00315593</th><td>4.19</td><td>4.38</td><td>4.10</td><td>4.08</td><td>4.07</td><td>4.35</td><td>4.04</td><td>3.80</td><td>4.12</td><td>4.21</td><td>⋯</td><td>3.91</td><td>3.98</td><td>3.90</td><td>4.08</td><td>4.12</td><td>3.91</td><td>3.90</td><td>4.00</td><td>3.77</td><td>4.17</td></tr>
	<tr><th scope=row>A_19_P00315601</th><td>8.43</td><td>8.58</td><td>7.09</td><td>8.83</td><td>9.05</td><td>7.90</td><td>8.24</td><td>7.73</td><td>7.99</td><td>8.23</td><td>⋯</td><td>8.88</td><td>8.64</td><td>8.61</td><td>8.92</td><td>8.84</td><td>8.55</td><td>8.71</td><td>8.43</td><td>8.79</td><td>8.73</td></tr>
	<tr><th scope=row>A_19_P00315603</th><td>7.99</td><td>8.11</td><td>6.49</td><td>8.22</td><td>8.48</td><td>7.66</td><td>7.89</td><td>6.85</td><td>7.65</td><td>7.84</td><td>⋯</td><td>8.39</td><td>8.08</td><td>8.21</td><td>8.50</td><td>8.39</td><td>8.02</td><td>8.36</td><td>8.00</td><td>8.39</td><td>8.25</td></tr>
	<tr><th scope=row>A_19_P00315625</th><td>3.88</td><td>3.62</td><td>3.74</td><td>3.91</td><td>4.07</td><td>3.84</td><td>3.73</td><td>4.05</td><td>3.89</td><td>4.17</td><td>⋯</td><td>3.95</td><td>3.80</td><td>4.05</td><td>3.80</td><td>3.88</td><td>3.93</td><td>3.89</td><td>3.95</td><td>3.67</td><td>3.75</td></tr>
	<tr><th scope=row>A_19_P00315627</th><td>4.47</td><td>4.28</td><td>4.40</td><td>4.41</td><td>4.81</td><td>4.04</td><td>4.28</td><td>4.36</td><td>4.40</td><td>5.27</td><td>⋯</td><td>4.34</td><td>4.68</td><td>4.28</td><td>4.40</td><td>4.17</td><td>4.09</td><td>4.57</td><td>5.03</td><td>4.39</td><td>4.31</td></tr>
	<tr><th scope=row>A_19_P00315631</th><td>4.61</td><td>4.65</td><td>5.40</td><td>3.76</td><td>4.21</td><td>4.03</td><td>3.87</td><td>4.14</td><td>3.80</td><td>3.75</td><td>⋯</td><td>4.14</td><td>4.62</td><td>5.00</td><td>3.89</td><td>4.04</td><td>3.87</td><td>3.63</td><td>4.57</td><td>4.03</td><td>3.94</td></tr>
	<tr><th scope=row>A_19_P00315633</th><td>5.01</td><td>4.74</td><td>5.89</td><td>3.45</td><td>4.51</td><td>4.00</td><td>3.77</td><td>4.23</td><td>3.84</td><td>3.72</td><td>⋯</td><td>4.20</td><td>5.16</td><td>5.52</td><td>4.07</td><td>4.05</td><td>3.84</td><td>3.80</td><td>4.97</td><td>4.07</td><td>3.84</td></tr>
	<tr><th scope=row>A_19_P00315641</th><td>4.20</td><td>4.52</td><td>5.12</td><td>4.99</td><td>4.49</td><td>3.95</td><td>4.92</td><td>5.05</td><td>5.11</td><td>4.67</td><td>⋯</td><td>4.84</td><td>4.92</td><td>4.67</td><td>4.81</td><td>4.66</td><td>4.28</td><td>5.34</td><td>4.80</td><td>4.53</td><td>4.66</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>A_33_P3880302</th><td>4.64</td><td>4.86</td><td>4.90</td><td>4.39</td><td>4.26</td><td>5.87</td><td>4.81</td><td>4.55</td><td>4.50</td><td>4.83</td><td>⋯</td><td>4.25</td><td>4.62</td><td>4.83</td><td>4.38</td><td>4.73</td><td>5.30</td><td>4.19</td><td>4.52</td><td>4.69</td><td>4.69</td></tr>
	<tr><th scope=row>A_33_P3881056</th><td>6.92</td><td>6.82</td><td>5.52</td><td>7.03</td><td>7.28</td><td>6.51</td><td>6.69</td><td>5.48</td><td>6.19</td><td>6.01</td><td>⋯</td><td>7.00</td><td>7.39</td><td>6.74</td><td>6.69</td><td>6.67</td><td>6.31</td><td>6.56</td><td>6.85</td><td>7.42</td><td>6.65</td></tr>
	<tr><th scope=row>A_33_P3881262</th><td>7.51</td><td>6.56</td><td>5.50</td><td>6.76</td><td>5.84</td><td>8.74</td><td>7.86</td><td>5.29</td><td>5.57</td><td>6.34</td><td>⋯</td><td>7.30</td><td>6.35</td><td>7.79</td><td>7.03</td><td>7.26</td><td>7.51</td><td>5.52</td><td>6.24</td><td>6.78</td><td>6.48</td></tr>
	<tr><th scope=row>A_33_P3881812</th><td>3.88</td><td>4.29</td><td>3.90</td><td>3.86</td><td>4.09</td><td>4.07</td><td>4.17</td><td>4.39</td><td>3.84</td><td>3.93</td><td>⋯</td><td>3.89</td><td>4.08</td><td>3.91</td><td>3.84</td><td>4.10</td><td>4.05</td><td>3.90</td><td>4.04</td><td>4.06</td><td>4.04</td></tr>
	<tr><th scope=row>A_33_P3882264</th><td>5.04</td><td>4.54</td><td>4.91</td><td>4.60</td><td>4.61</td><td>5.54</td><td>4.58</td><td>4.78</td><td>5.19</td><td>4.86</td><td>⋯</td><td>4.86</td><td>4.88</td><td>4.65</td><td>5.03</td><td>4.89</td><td>4.75</td><td>4.45</td><td>4.83</td><td>4.51</td><td>4.78</td></tr>
	<tr><th scope=row>A_33_P3882659</th><td>8.20</td><td>7.92</td><td>8.34</td><td>7.89</td><td>7.43</td><td>8.21</td><td>8.01</td><td>7.80</td><td>7.87</td><td>7.87</td><td>⋯</td><td>7.37</td><td>7.90</td><td>7.91</td><td>7.92</td><td>7.77</td><td>7.96</td><td>7.90</td><td>7.84</td><td>7.53</td><td>7.68</td></tr>
	<tr><th scope=row>A_33_P3883116</th><td>3.96</td><td>4.41</td><td>3.81</td><td>3.95</td><td>3.71</td><td>3.96</td><td>3.88</td><td>3.93</td><td>4.62</td><td>3.85</td><td>⋯</td><td>4.25</td><td>3.97</td><td>3.72</td><td>3.95</td><td>3.66</td><td>3.83</td><td>4.12</td><td>3.58</td><td>3.94</td><td>4.03</td></tr>
	<tr><th scope=row>A_33_P3883912</th><td>6.68</td><td>6.23</td><td>6.99</td><td>6.74</td><td>6.30</td><td>5.92</td><td>6.07</td><td>6.46</td><td>6.78</td><td>6.55</td><td>⋯</td><td>6.90</td><td>6.15</td><td>6.31</td><td>6.31</td><td>6.16</td><td>6.28</td><td>6.38</td><td>6.53</td><td>6.14</td><td>6.07</td></tr>
	<tr><th scope=row>A_33_P3883985</th><td>6.78</td><td>6.03</td><td>5.85</td><td>7.33</td><td>7.10</td><td>6.33</td><td>6.96</td><td>7.35</td><td>6.58</td><td>7.55</td><td>⋯</td><td>6.18</td><td>6.75</td><td>6.42</td><td>6.64</td><td>5.67</td><td>6.63</td><td>7.09</td><td>7.07</td><td>7.00</td><td>5.84</td></tr>
	<tr><th scope=row>A_33_P3884005</th><td>3.94</td><td>3.75</td><td>3.92</td><td>4.30</td><td>4.12</td><td>3.98</td><td>3.95</td><td>4.01</td><td>4.17</td><td>4.11</td><td>⋯</td><td>4.21</td><td>3.91</td><td>3.92</td><td>3.90</td><td>4.16</td><td>4.03</td><td>4.20</td><td>4.25</td><td>3.65</td><td>4.09</td></tr>
	<tr><th scope=row>A_33_P3884179</th><td>6.23</td><td>6.34</td><td>4.59</td><td>6.46</td><td>6.27</td><td>6.10</td><td>6.01</td><td>5.40</td><td>5.73</td><td>6.17</td><td>⋯</td><td>6.14</td><td>6.28</td><td>6.03</td><td>6.52</td><td>6.12</td><td>6.54</td><td>6.08</td><td>6.11</td><td>6.23</td><td>6.49</td></tr>
	<tr><th scope=row>A_33_P3884230</th><td>6.69</td><td>6.27</td><td>6.30</td><td>6.95</td><td>6.58</td><td>5.44</td><td>6.69</td><td>7.07</td><td>6.89</td><td>6.70</td><td>⋯</td><td>6.43</td><td>6.77</td><td>6.62</td><td>6.44</td><td>6.75</td><td>5.88</td><td>6.73</td><td>6.84</td><td>6.65</td><td>6.78</td></tr>
	<tr><th scope=row>A_33_P3884610</th><td>3.55</td><td>3.97</td><td>3.53</td><td>3.78</td><td>3.63</td><td>3.57</td><td>3.70</td><td>3.84</td><td>3.78</td><td>3.76</td><td>⋯</td><td>3.65</td><td>3.67</td><td>3.88</td><td>3.85</td><td>3.61</td><td>3.78</td><td>3.71</td><td>3.52</td><td>3.60</td><td>3.39</td></tr>
	<tr><th scope=row>A_33_P3885084</th><td>3.87</td><td>3.66</td><td>3.62</td><td>4.33</td><td>3.84</td><td>3.83</td><td>3.69</td><td>3.47</td><td>3.72</td><td>3.80</td><td>⋯</td><td>3.77</td><td>3.87</td><td>3.81</td><td>3.85</td><td>3.81</td><td>3.99</td><td>3.77</td><td>3.82</td><td>3.82</td><td>3.87</td></tr>
	<tr><th scope=row>A_33_P3886707</th><td>7.72</td><td>7.97</td><td>7.96</td><td>7.94</td><td>7.90</td><td>7.79</td><td>7.99</td><td>8.11</td><td>8.27</td><td>7.94</td><td>⋯</td><td>7.96</td><td>7.82</td><td>7.92</td><td>7.84</td><td>7.87</td><td>7.92</td><td>7.98</td><td>7.69</td><td>8.00</td><td>7.81</td></tr>
	<tr><th scope=row>A_33_P3886737</th><td>3.99</td><td>4.19</td><td>4.05</td><td>3.80</td><td>4.08</td><td>4.08</td><td>3.90</td><td>3.94</td><td>3.89</td><td>3.89</td><td>⋯</td><td>4.06</td><td>4.41</td><td>4.09</td><td>3.72</td><td>3.74</td><td>4.07</td><td>4.06</td><td>3.77</td><td>3.91</td><td>3.73</td></tr>
	<tr><th scope=row>A_33_P3886938</th><td>3.88</td><td>4.07</td><td>4.09</td><td>3.61</td><td>3.51</td><td>3.69</td><td>3.65</td><td>3.57</td><td>3.72</td><td>3.49</td><td>⋯</td><td>3.94</td><td>3.47</td><td>3.90</td><td>3.49</td><td>3.75</td><td>3.75</td><td>3.74</td><td>3.89</td><td>3.86</td><td>3.73</td></tr>
	<tr><th scope=row>A_33_P3887081</th><td>4.79</td><td>4.25</td><td>4.29</td><td>4.14</td><td>5.26</td><td>4.21</td><td>4.21</td><td>4.41</td><td>4.21</td><td>4.50</td><td>⋯</td><td>4.62</td><td>4.53</td><td>4.51</td><td>4.64</td><td>4.12</td><td>4.26</td><td>4.27</td><td>4.51</td><td>4.81</td><td>4.02</td></tr>
	<tr><th scope=row>A_33_P3887888</th><td>3.93</td><td>4.16</td><td>4.00</td><td>4.26</td><td>4.64</td><td>4.34</td><td>4.39</td><td>4.61</td><td>5.38</td><td>4.71</td><td>⋯</td><td>4.40</td><td>4.11</td><td>4.22</td><td>4.42</td><td>4.20</td><td>4.79</td><td>4.83</td><td>4.53</td><td>4.28</td><td>4.65</td></tr>
	<tr><th scope=row>A_33_P3888365</th><td>6.44</td><td>6.39</td><td>6.34</td><td>6.67</td><td>6.42</td><td>5.91</td><td>6.54</td><td>6.78</td><td>6.44</td><td>6.09</td><td>⋯</td><td>6.38</td><td>6.40</td><td>6.36</td><td>6.33</td><td>6.22</td><td>6.17</td><td>6.58</td><td>6.30</td><td>6.41</td><td>6.37</td></tr>
	<tr><th scope=row>A_33_P3888568</th><td>4.00</td><td>4.05</td><td>3.96</td><td>4.02</td><td>3.98</td><td>3.87</td><td>4.03</td><td>3.96</td><td>4.09</td><td>4.00</td><td>⋯</td><td>4.22</td><td>3.90</td><td>3.88</td><td>3.94</td><td>4.02</td><td>3.98</td><td>4.12</td><td>4.02</td><td>3.97</td><td>4.01</td></tr>
	<tr><th scope=row>A_33_P3888629</th><td>5.25</td><td>5.75</td><td>5.48</td><td>6.98</td><td>5.74</td><td>5.09</td><td>6.19</td><td>6.07</td><td>6.88</td><td>6.63</td><td>⋯</td><td>5.17</td><td>5.69</td><td>5.78</td><td>5.29</td><td>5.31</td><td>5.66</td><td>7.38</td><td>5.76</td><td>6.32</td><td>5.16</td></tr>
	<tr><th scope=row>A_33_P3889179</th><td>3.96</td><td>3.70</td><td>4.04</td><td>3.72</td><td>3.77</td><td>4.31</td><td>3.99</td><td>3.83</td><td>3.67</td><td>3.88</td><td>⋯</td><td>3.92</td><td>4.03</td><td>3.85</td><td>3.78</td><td>3.98</td><td>3.96</td><td>3.67</td><td>3.71</td><td>3.73</td><td>3.89</td></tr>
	<tr><th scope=row>A_33_P3889634</th><td>3.64</td><td>3.81</td><td>3.83</td><td>3.67</td><td>3.76</td><td>3.94</td><td>3.74</td><td>3.58</td><td>3.59</td><td>3.52</td><td>⋯</td><td>3.52</td><td>3.83</td><td>3.78</td><td>3.65</td><td>3.63</td><td>3.68</td><td>3.71</td><td>3.76</td><td>3.69</td><td>3.62</td></tr>
	<tr><th scope=row>A_33_P3891747</th><td>3.68</td><td>3.65</td><td>3.81</td><td>4.05</td><td>4.15</td><td>3.95</td><td>4.01</td><td>3.95</td><td>3.88</td><td>3.88</td><td>⋯</td><td>4.03</td><td>3.96</td><td>3.75</td><td>4.00</td><td>3.88</td><td>3.93</td><td>3.66</td><td>3.91</td><td>3.90</td><td>4.07</td></tr>
	<tr><th scope=row>A_33_P3892608</th><td>5.07</td><td>5.20</td><td>5.00</td><td>4.75</td><td>5.47</td><td>4.72</td><td>5.40</td><td>5.22</td><td>4.32</td><td>4.90</td><td>⋯</td><td>5.64</td><td>5.19</td><td>5.30</td><td>5.54</td><td>5.58</td><td>5.01</td><td>4.45</td><td>5.33</td><td>5.22</td><td>5.50</td></tr>
	<tr><th scope=row>A_33_P3892710</th><td>3.95</td><td>4.34</td><td>3.96</td><td>4.30</td><td>4.24</td><td>4.16</td><td>4.02</td><td>3.91</td><td>3.88</td><td>4.08</td><td>⋯</td><td>3.81</td><td>3.82</td><td>3.92</td><td>4.20</td><td>3.76</td><td>4.24</td><td>4.06</td><td>4.04</td><td>3.92</td><td>3.98</td></tr>
	<tr><th scope=row>A_33_P3893191</th><td>6.90</td><td>6.57</td><td>5.33</td><td>6.30</td><td>6.70</td><td>6.50</td><td>6.53</td><td>6.13</td><td>6.33</td><td>6.26</td><td>⋯</td><td>5.93</td><td>6.06</td><td>6.19</td><td>6.80</td><td>6.04</td><td>6.59</td><td>6.31</td><td>6.62</td><td>6.64</td><td>6.52</td></tr>
	<tr><th scope=row>A_33_P3897734</th><td>3.80</td><td>3.87</td><td>3.81</td><td>3.81</td><td>3.75</td><td>3.59</td><td>3.97</td><td>3.57</td><td>3.72</td><td>3.98</td><td>⋯</td><td>3.82</td><td>3.72</td><td>3.69</td><td>3.75</td><td>3.75</td><td>3.59</td><td>3.66</td><td>3.72</td><td>3.93</td><td>3.90</td></tr>
	<tr><th scope=row>A_33_P3901921</th><td>3.71</td><td>3.92</td><td>4.15</td><td>3.90</td><td>3.87</td><td>3.97</td><td>4.03</td><td>3.82</td><td>4.07</td><td>3.87</td><td>⋯</td><td>4.06</td><td>4.10</td><td>3.98</td><td>3.73</td><td>3.96</td><td>4.01</td><td>4.05</td><td>4.09</td><td>3.94</td><td>4.15</td></tr>
</tbody>
</table>




```R
ids <- AnnoProbe::idmap('GPL17077',type = 'soft')
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
<caption>A data.frame: 47264 × 106</caption>
<thead>
	<tr><th scope=col>Probe</th><th scope=col>GSM2696609</th><th scope=col>GSM2696610</th><th scope=col>GSM2696611</th><th scope=col>GSM2696612</th><th scope=col>GSM2696613</th><th scope=col>GSM2696614</th><th scope=col>GSM2696615</th><th scope=col>GSM2696616</th><th scope=col>GSM2696617</th><th scope=col>⋯</th><th scope=col>GSM2696704</th><th scope=col>GSM2696705</th><th scope=col>GSM2696706</th><th scope=col>GSM2696707</th><th scope=col>GSM2696708</th><th scope=col>GSM2696709</th><th scope=col>GSM2696710</th><th scope=col>GSM2696711</th><th scope=col>GSM2696712</th><th scope=col>symbol</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>A_19_P00315452</td><td>4.45</td><td>4.47</td><td>5.16</td><td>4.96</td><td>4.76</td><td>4.93</td><td>5.03</td><td>5.54</td><td>4.50</td><td>⋯</td><td>4.81</td><td>4.82</td><td>4.70</td><td>4.55</td><td>6.17</td><td>4.57</td><td>4.65</td><td>4.62</td><td>4.69</td><td>LOC100130938  </td></tr>
	<tr><td>A_19_P00315459</td><td>4.99</td><td>5.12</td><td>4.56</td><td>5.02</td><td>4.88</td><td>5.13</td><td>5.03</td><td>4.92</td><td>4.47</td><td>⋯</td><td>4.75</td><td>4.79</td><td>5.14</td><td>4.86</td><td>4.88</td><td>5.30</td><td>5.00</td><td>4.68</td><td>4.86</td><td>              </td></tr>
	<tr><td>A_19_P00315482</td><td>4.00</td><td>4.06</td><td>3.92</td><td>3.93</td><td>3.96</td><td>3.95</td><td>3.88</td><td>3.87</td><td>3.83</td><td>⋯</td><td>4.05</td><td>3.82</td><td>4.03</td><td>3.90</td><td>4.08</td><td>3.70</td><td>3.89</td><td>3.92</td><td>3.77</td><td>LOC729177     </td></tr>
	<tr><td>A_19_P00315492</td><td>4.13</td><td>4.41</td><td>4.15</td><td>3.83</td><td>4.11</td><td>3.99</td><td>3.91</td><td>4.09</td><td>3.91</td><td>⋯</td><td>4.27</td><td>3.83</td><td>3.94</td><td>3.91</td><td>3.76</td><td>3.96</td><td>4.11</td><td>3.99</td><td>3.83</td><td>Q73P46        </td></tr>
	<tr><td>A_19_P00315493</td><td>6.30</td><td>5.21</td><td>5.10</td><td>5.35</td><td>6.08</td><td>5.82</td><td>6.50</td><td>4.89</td><td>5.22</td><td>⋯</td><td>5.57</td><td>5.39</td><td>5.63</td><td>5.21</td><td>5.36</td><td>5.38</td><td>5.44</td><td>6.40</td><td>5.57</td><td>LOC145474     </td></tr>
	<tr><td>A_19_P00315502</td><td>3.91</td><td>3.61</td><td>4.24</td><td>3.76</td><td>3.84</td><td>3.89</td><td>3.77</td><td>3.79</td><td>3.82</td><td>⋯</td><td>3.87</td><td>3.94</td><td>3.68</td><td>3.73</td><td>3.72</td><td>3.92</td><td>3.82</td><td>3.82</td><td>3.83</td><td>P01115        </td></tr>
	<tr><td>A_19_P00315506</td><td>6.38</td><td>6.33</td><td>6.70</td><td>6.26</td><td>6.52</td><td>6.58</td><td>6.69</td><td>6.26</td><td>7.22</td><td>⋯</td><td>6.56</td><td>6.82</td><td>6.70</td><td>7.04</td><td>6.32</td><td>6.66</td><td>7.12</td><td>6.33</td><td>6.81</td><td>KIAA0040      </td></tr>
	<tr><td>A_19_P00315518</td><td>3.93</td><td>3.48</td><td>3.75</td><td>3.67</td><td>3.84</td><td>3.82</td><td>3.82</td><td>3.88</td><td>3.65</td><td>⋯</td><td>3.74</td><td>3.66</td><td>3.87</td><td>3.86</td><td>3.84</td><td>3.75</td><td>3.85</td><td>3.81</td><td>3.81</td><td>XLOC_005810   </td></tr>
	<tr><td>A_19_P00315519</td><td>3.40</td><td>3.57</td><td>3.66</td><td>3.58</td><td>3.71</td><td>3.78</td><td>3.78</td><td>3.84</td><td>3.59</td><td>⋯</td><td>3.63</td><td>3.38</td><td>3.68</td><td>3.67</td><td>3.51</td><td>3.67</td><td>3.82</td><td>3.74</td><td>3.60</td><td>XLOC_004914   </td></tr>
	<tr><td>A_19_P00315524</td><td>5.25</td><td>5.71</td><td>4.60</td><td>5.15</td><td>5.78</td><td>6.43</td><td>6.78</td><td>4.88</td><td>4.43</td><td>⋯</td><td>5.36</td><td>5.68</td><td>6.06</td><td>6.61</td><td>6.53</td><td>4.41</td><td>5.32</td><td>6.58</td><td>6.37</td><td>MIAT          </td></tr>
	<tr><td>A_19_P00315528</td><td>3.94</td><td>4.41</td><td>4.09</td><td>4.03</td><td>4.08</td><td>4.03</td><td>4.33</td><td>3.98</td><td>3.93</td><td>⋯</td><td>3.83</td><td>4.02</td><td>4.07</td><td>4.15</td><td>4.39</td><td>3.91</td><td>4.03</td><td>4.21</td><td>4.19</td><td>XLOC_008370   </td></tr>
	<tr><td>A_19_P00315529</td><td>3.86</td><td>4.44</td><td>3.85</td><td>4.10</td><td>4.19</td><td>4.05</td><td>4.41</td><td>4.30</td><td>4.01</td><td>⋯</td><td>4.06</td><td>3.91</td><td>4.22</td><td>4.10</td><td>4.36</td><td>4.31</td><td>4.38</td><td>4.17</td><td>4.06</td><td>XLOC_008370   </td></tr>
	<tr><td>A_19_P00315538</td><td>3.69</td><td>4.18</td><td>3.76</td><td>3.51</td><td>3.83</td><td>4.46</td><td>3.87</td><td>3.83</td><td>3.91</td><td>⋯</td><td>3.82</td><td>3.93</td><td>3.61</td><td>3.72</td><td>3.72</td><td>3.81</td><td>3.61</td><td>3.88</td><td>3.76</td><td>              </td></tr>
	<tr><td>A_19_P00315541</td><td>3.73</td><td>3.77</td><td>3.80</td><td>3.88</td><td>3.53</td><td>3.75</td><td>3.61</td><td>3.58</td><td>3.79</td><td>⋯</td><td>4.03</td><td>3.80</td><td>3.78</td><td>3.71</td><td>3.69</td><td>3.55</td><td>3.69</td><td>3.76</td><td>3.83</td><td>LOC100507428  </td></tr>
	<tr><td>A_19_P00315543</td><td>4.78</td><td>4.42</td><td>5.36</td><td>4.76</td><td>4.88</td><td>4.29</td><td>4.20</td><td>5.43</td><td>4.57</td><td>⋯</td><td>4.49</td><td>4.69</td><td>4.91</td><td>4.85</td><td>4.65</td><td>4.62</td><td>4.73</td><td>4.79</td><td>4.80</td><td>LOC100132354  </td></tr>
	<tr><td>A_19_P00315550</td><td>6.81</td><td>6.45</td><td>6.44</td><td>5.41</td><td>6.27</td><td>5.79</td><td>5.59</td><td>6.74</td><td>5.25</td><td>⋯</td><td>6.53</td><td>6.40</td><td>6.82</td><td>6.57</td><td>5.58</td><td>5.06</td><td>6.50</td><td>5.75</td><td>6.44</td><td>LOC400043     </td></tr>
	<tr><td>A_19_P00315551</td><td>6.88</td><td>6.50</td><td>6.08</td><td>5.93</td><td>6.38</td><td>6.08</td><td>5.80</td><td>6.40</td><td>5.58</td><td>⋯</td><td>6.41</td><td>6.40</td><td>6.76</td><td>6.34</td><td>5.99</td><td>5.31</td><td>6.58</td><td>6.03</td><td>6.42</td><td>LOC400043     </td></tr>
	<tr><td>A_19_P00315554</td><td>4.06</td><td>3.94</td><td>3.77</td><td>3.68</td><td>3.92</td><td>5.27</td><td>4.14</td><td>3.77</td><td>3.93</td><td>⋯</td><td>3.86</td><td>4.08</td><td>3.99</td><td>4.33</td><td>4.29</td><td>3.85</td><td>3.85</td><td>4.11</td><td>4.20</td><td>XLOC_006756   </td></tr>
	<tr><td>A_19_P00315581</td><td>7.84</td><td>7.93</td><td>7.87</td><td>7.92</td><td>7.87</td><td>7.89</td><td>8.01</td><td>7.19</td><td>7.58</td><td>⋯</td><td>8.19</td><td>7.95</td><td>7.99</td><td>8.02</td><td>7.84</td><td>7.51</td><td>7.93</td><td>7.96</td><td>7.85</td><td>FLJ30838      </td></tr>
	<tr><td>A_19_P00315583</td><td>4.77</td><td>4.94</td><td>4.05</td><td>4.80</td><td>4.88</td><td>5.15</td><td>4.58</td><td>4.64</td><td>4.21</td><td>⋯</td><td>4.87</td><td>4.94</td><td>4.87</td><td>4.73</td><td>4.70</td><td>4.84</td><td>4.59</td><td>4.78</td><td>4.87</td><td>LOC100130691  </td></tr>
	<tr><td>A_19_P00315584</td><td>5.44</td><td>5.50</td><td>4.68</td><td>5.33</td><td>5.48</td><td>5.82</td><td>5.28</td><td>5.29</td><td>4.55</td><td>⋯</td><td>5.58</td><td>5.63</td><td>5.50</td><td>5.31</td><td>5.32</td><td>5.28</td><td>5.07</td><td>5.53</td><td>5.60</td><td>LOC100130691  </td></tr>
	<tr><td>A_19_P00315587</td><td>3.70</td><td>3.71</td><td>3.91</td><td>3.85</td><td>4.13</td><td>3.90</td><td>3.68</td><td>3.78</td><td>3.83</td><td>⋯</td><td>3.88</td><td>3.42</td><td>3.82</td><td>3.77</td><td>3.79</td><td>3.67</td><td>3.75</td><td>3.86</td><td>4.05</td><td>LOC100506924  </td></tr>
	<tr><td>A_19_P00315593</td><td>4.19</td><td>4.38</td><td>4.10</td><td>4.08</td><td>4.07</td><td>4.35</td><td>4.04</td><td>3.80</td><td>4.12</td><td>⋯</td><td>3.98</td><td>3.90</td><td>4.08</td><td>4.12</td><td>3.91</td><td>3.90</td><td>4.00</td><td>3.77</td><td>4.17</td><td>XLOC_004643   </td></tr>
	<tr><td>A_19_P00315601</td><td>8.43</td><td>8.58</td><td>7.09</td><td>8.83</td><td>9.05</td><td>7.90</td><td>8.24</td><td>7.73</td><td>7.99</td><td>⋯</td><td>8.64</td><td>8.61</td><td>8.92</td><td>8.84</td><td>8.55</td><td>8.71</td><td>8.43</td><td>8.79</td><td>8.73</td><td>XLOC_l2_015760</td></tr>
	<tr><td>A_19_P00315603</td><td>7.99</td><td>8.11</td><td>6.49</td><td>8.22</td><td>8.48</td><td>7.66</td><td>7.89</td><td>6.85</td><td>7.65</td><td>⋯</td><td>8.08</td><td>8.21</td><td>8.50</td><td>8.39</td><td>8.02</td><td>8.36</td><td>8.00</td><td>8.39</td><td>8.25</td><td>XLOC_l2_015760</td></tr>
	<tr><td>A_19_P00315625</td><td>3.88</td><td>3.62</td><td>3.74</td><td>3.91</td><td>4.07</td><td>3.84</td><td>3.73</td><td>4.05</td><td>3.89</td><td>⋯</td><td>3.80</td><td>4.05</td><td>3.80</td><td>3.88</td><td>3.93</td><td>3.89</td><td>3.95</td><td>3.67</td><td>3.75</td><td>XLOC_005441   </td></tr>
	<tr><td>A_19_P00315627</td><td>4.47</td><td>4.28</td><td>4.40</td><td>4.41</td><td>4.81</td><td>4.04</td><td>4.28</td><td>4.36</td><td>4.40</td><td>⋯</td><td>4.68</td><td>4.28</td><td>4.40</td><td>4.17</td><td>4.09</td><td>4.57</td><td>5.03</td><td>4.39</td><td>4.31</td><td>XLOC_l2_015098</td></tr>
	<tr><td>A_19_P00315631</td><td>4.61</td><td>4.65</td><td>5.40</td><td>3.76</td><td>4.21</td><td>4.03</td><td>3.87</td><td>4.14</td><td>3.80</td><td>⋯</td><td>4.62</td><td>5.00</td><td>3.89</td><td>4.04</td><td>3.87</td><td>3.63</td><td>4.57</td><td>4.03</td><td>3.94</td><td>XLOC_l2_013734</td></tr>
	<tr><td>A_19_P00315633</td><td>5.01</td><td>4.74</td><td>5.89</td><td>3.45</td><td>4.51</td><td>4.00</td><td>3.77</td><td>4.23</td><td>3.84</td><td>⋯</td><td>5.16</td><td>5.52</td><td>4.07</td><td>4.05</td><td>3.84</td><td>3.80</td><td>4.97</td><td>4.07</td><td>3.84</td><td>              </td></tr>
	<tr><td>A_19_P00315641</td><td>4.20</td><td>4.52</td><td>5.12</td><td>4.99</td><td>4.49</td><td>3.95</td><td>4.92</td><td>5.05</td><td>5.11</td><td>⋯</td><td>4.92</td><td>4.67</td><td>4.81</td><td>4.66</td><td>4.28</td><td>5.34</td><td>4.80</td><td>4.53</td><td>4.66</td><td>XLOC_008079   </td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>A_33_P3880302</td><td>4.64</td><td>4.86</td><td>4.90</td><td>4.39</td><td>4.26</td><td>5.87</td><td>4.81</td><td>4.55</td><td>4.50</td><td>⋯</td><td>4.62</td><td>4.83</td><td>4.38</td><td>4.73</td><td>5.30</td><td>4.19</td><td>4.52</td><td>4.69</td><td>4.69</td><td>EPHB2         </td></tr>
	<tr><td>A_33_P3881056</td><td>6.92</td><td>6.82</td><td>5.52</td><td>7.03</td><td>7.28</td><td>6.51</td><td>6.69</td><td>5.48</td><td>6.19</td><td>⋯</td><td>7.39</td><td>6.74</td><td>6.69</td><td>6.67</td><td>6.31</td><td>6.56</td><td>6.85</td><td>7.42</td><td>6.65</td><td>              </td></tr>
	<tr><td>A_33_P3881262</td><td>7.51</td><td>6.56</td><td>5.50</td><td>6.76</td><td>5.84</td><td>8.74</td><td>7.86</td><td>5.29</td><td>5.57</td><td>⋯</td><td>6.35</td><td>7.79</td><td>7.03</td><td>7.26</td><td>7.51</td><td>5.52</td><td>6.24</td><td>6.78</td><td>6.48</td><td>CSF3R         </td></tr>
	<tr><td>A_33_P3881812</td><td>3.88</td><td>4.29</td><td>3.90</td><td>3.86</td><td>4.09</td><td>4.07</td><td>4.17</td><td>4.39</td><td>3.84</td><td>⋯</td><td>4.08</td><td>3.91</td><td>3.84</td><td>4.10</td><td>4.05</td><td>3.90</td><td>4.04</td><td>4.06</td><td>4.04</td><td>IRGM          </td></tr>
	<tr><td>A_33_P3882264</td><td>5.04</td><td>4.54</td><td>4.91</td><td>4.60</td><td>4.61</td><td>5.54</td><td>4.58</td><td>4.78</td><td>5.19</td><td>⋯</td><td>4.88</td><td>4.65</td><td>5.03</td><td>4.89</td><td>4.75</td><td>4.45</td><td>4.83</td><td>4.51</td><td>4.78</td><td>GRASPOS       </td></tr>
	<tr><td>A_33_P3882659</td><td>8.20</td><td>7.92</td><td>8.34</td><td>7.89</td><td>7.43</td><td>8.21</td><td>8.01</td><td>7.80</td><td>7.87</td><td>⋯</td><td>7.90</td><td>7.91</td><td>7.92</td><td>7.77</td><td>7.96</td><td>7.90</td><td>7.84</td><td>7.53</td><td>7.68</td><td>HSP90AB5P     </td></tr>
	<tr><td>A_33_P3883116</td><td>3.96</td><td>4.41</td><td>3.81</td><td>3.95</td><td>3.71</td><td>3.96</td><td>3.88</td><td>3.93</td><td>4.62</td><td>⋯</td><td>3.97</td><td>3.72</td><td>3.95</td><td>3.66</td><td>3.83</td><td>4.12</td><td>3.58</td><td>3.94</td><td>4.03</td><td>LOC253962     </td></tr>
	<tr><td>A_33_P3883912</td><td>6.68</td><td>6.23</td><td>6.99</td><td>6.74</td><td>6.30</td><td>5.92</td><td>6.07</td><td>6.46</td><td>6.78</td><td>⋯</td><td>6.15</td><td>6.31</td><td>6.31</td><td>6.16</td><td>6.28</td><td>6.38</td><td>6.53</td><td>6.14</td><td>6.07</td><td>ZCCHC10       </td></tr>
	<tr><td>A_33_P3883985</td><td>6.78</td><td>6.03</td><td>5.85</td><td>7.33</td><td>7.10</td><td>6.33</td><td>6.96</td><td>7.35</td><td>6.58</td><td>⋯</td><td>6.75</td><td>6.42</td><td>6.64</td><td>5.67</td><td>6.63</td><td>7.09</td><td>7.07</td><td>7.00</td><td>5.84</td><td>LMF1          </td></tr>
	<tr><td>A_33_P3884005</td><td>3.94</td><td>3.75</td><td>3.92</td><td>4.30</td><td>4.12</td><td>3.98</td><td>3.95</td><td>4.01</td><td>4.17</td><td>⋯</td><td>3.91</td><td>3.92</td><td>3.90</td><td>4.16</td><td>4.03</td><td>4.20</td><td>4.25</td><td>3.65</td><td>4.09</td><td>LOC93444      </td></tr>
	<tr><td>A_33_P3884179</td><td>6.23</td><td>6.34</td><td>4.59</td><td>6.46</td><td>6.27</td><td>6.10</td><td>6.01</td><td>5.40</td><td>5.73</td><td>⋯</td><td>6.28</td><td>6.03</td><td>6.52</td><td>6.12</td><td>6.54</td><td>6.08</td><td>6.11</td><td>6.23</td><td>6.49</td><td>LOC100506123  </td></tr>
	<tr><td>A_33_P3884230</td><td>6.69</td><td>6.27</td><td>6.30</td><td>6.95</td><td>6.58</td><td>5.44</td><td>6.69</td><td>7.07</td><td>6.89</td><td>⋯</td><td>6.77</td><td>6.62</td><td>6.44</td><td>6.75</td><td>5.88</td><td>6.73</td><td>6.84</td><td>6.65</td><td>6.78</td><td>NFIX          </td></tr>
	<tr><td>A_33_P3884610</td><td>3.55</td><td>3.97</td><td>3.53</td><td>3.78</td><td>3.63</td><td>3.57</td><td>3.70</td><td>3.84</td><td>3.78</td><td>⋯</td><td>3.67</td><td>3.88</td><td>3.85</td><td>3.61</td><td>3.78</td><td>3.71</td><td>3.52</td><td>3.60</td><td>3.39</td><td>PIP5K1P1      </td></tr>
	<tr><td>A_33_P3885084</td><td>3.87</td><td>3.66</td><td>3.62</td><td>4.33</td><td>3.84</td><td>3.83</td><td>3.69</td><td>3.47</td><td>3.72</td><td>⋯</td><td>3.87</td><td>3.81</td><td>3.85</td><td>3.81</td><td>3.99</td><td>3.77</td><td>3.82</td><td>3.82</td><td>3.87</td><td>LOC147004     </td></tr>
	<tr><td>A_33_P3886707</td><td>7.72</td><td>7.97</td><td>7.96</td><td>7.94</td><td>7.90</td><td>7.79</td><td>7.99</td><td>8.11</td><td>8.27</td><td>⋯</td><td>7.82</td><td>7.92</td><td>7.84</td><td>7.87</td><td>7.92</td><td>7.98</td><td>7.69</td><td>8.00</td><td>7.81</td><td>PPA2          </td></tr>
	<tr><td>A_33_P3886737</td><td>3.99</td><td>4.19</td><td>4.05</td><td>3.80</td><td>4.08</td><td>4.08</td><td>3.90</td><td>3.94</td><td>3.89</td><td>⋯</td><td>4.41</td><td>4.09</td><td>3.72</td><td>3.74</td><td>4.07</td><td>4.06</td><td>3.77</td><td>3.91</td><td>3.73</td><td>DKFZp686L13185</td></tr>
	<tr><td>A_33_P3886938</td><td>3.88</td><td>4.07</td><td>4.09</td><td>3.61</td><td>3.51</td><td>3.69</td><td>3.65</td><td>3.57</td><td>3.72</td><td>⋯</td><td>3.47</td><td>3.90</td><td>3.49</td><td>3.75</td><td>3.75</td><td>3.74</td><td>3.89</td><td>3.86</td><td>3.73</td><td>              </td></tr>
	<tr><td>A_33_P3887081</td><td>4.79</td><td>4.25</td><td>4.29</td><td>4.14</td><td>5.26</td><td>4.21</td><td>4.21</td><td>4.41</td><td>4.21</td><td>⋯</td><td>4.53</td><td>4.51</td><td>4.64</td><td>4.12</td><td>4.26</td><td>4.27</td><td>4.51</td><td>4.81</td><td>4.02</td><td>LOC338817     </td></tr>
	<tr><td>A_33_P3887888</td><td>3.93</td><td>4.16</td><td>4.00</td><td>4.26</td><td>4.64</td><td>4.34</td><td>4.39</td><td>4.61</td><td>5.38</td><td>⋯</td><td>4.11</td><td>4.22</td><td>4.42</td><td>4.20</td><td>4.79</td><td>4.83</td><td>4.53</td><td>4.28</td><td>4.65</td><td>FLJ44511      </td></tr>
	<tr><td>A_33_P3888365</td><td>6.44</td><td>6.39</td><td>6.34</td><td>6.67</td><td>6.42</td><td>5.91</td><td>6.54</td><td>6.78</td><td>6.44</td><td>⋯</td><td>6.40</td><td>6.36</td><td>6.33</td><td>6.22</td><td>6.17</td><td>6.58</td><td>6.30</td><td>6.41</td><td>6.37</td><td>RSBN1         </td></tr>
	<tr><td>A_33_P3888568</td><td>4.00</td><td>4.05</td><td>3.96</td><td>4.02</td><td>3.98</td><td>3.87</td><td>4.03</td><td>3.96</td><td>4.09</td><td>⋯</td><td>3.90</td><td>3.88</td><td>3.94</td><td>4.02</td><td>3.98</td><td>4.12</td><td>4.02</td><td>3.97</td><td>4.01</td><td>SNORA70B      </td></tr>
	<tr><td>A_33_P3888629</td><td>5.25</td><td>5.75</td><td>5.48</td><td>6.98</td><td>5.74</td><td>5.09</td><td>6.19</td><td>6.07</td><td>6.88</td><td>⋯</td><td>5.69</td><td>5.78</td><td>5.29</td><td>5.31</td><td>5.66</td><td>7.38</td><td>5.76</td><td>6.32</td><td>5.16</td><td>MECOM         </td></tr>
	<tr><td>A_33_P3889179</td><td>3.96</td><td>3.70</td><td>4.04</td><td>3.72</td><td>3.77</td><td>4.31</td><td>3.99</td><td>3.83</td><td>3.67</td><td>⋯</td><td>4.03</td><td>3.85</td><td>3.78</td><td>3.98</td><td>3.96</td><td>3.67</td><td>3.71</td><td>3.73</td><td>3.89</td><td>LOC285401     </td></tr>
	<tr><td>A_33_P3889634</td><td>3.64</td><td>3.81</td><td>3.83</td><td>3.67</td><td>3.76</td><td>3.94</td><td>3.74</td><td>3.58</td><td>3.59</td><td>⋯</td><td>3.83</td><td>3.78</td><td>3.65</td><td>3.63</td><td>3.68</td><td>3.71</td><td>3.76</td><td>3.69</td><td>3.62</td><td>LOC284798     </td></tr>
	<tr><td>A_33_P3891747</td><td>3.68</td><td>3.65</td><td>3.81</td><td>4.05</td><td>4.15</td><td>3.95</td><td>4.01</td><td>3.95</td><td>3.88</td><td>⋯</td><td>3.96</td><td>3.75</td><td>4.00</td><td>3.88</td><td>3.93</td><td>3.66</td><td>3.91</td><td>3.90</td><td>4.07</td><td>LOC284072     </td></tr>
	<tr><td>A_33_P3892608</td><td>5.07</td><td>5.20</td><td>5.00</td><td>4.75</td><td>5.47</td><td>4.72</td><td>5.40</td><td>5.22</td><td>4.32</td><td>⋯</td><td>5.19</td><td>5.30</td><td>5.54</td><td>5.58</td><td>5.01</td><td>4.45</td><td>5.33</td><td>5.22</td><td>5.50</td><td>KCNT2         </td></tr>
	<tr><td>A_33_P3892710</td><td>3.95</td><td>4.34</td><td>3.96</td><td>4.30</td><td>4.24</td><td>4.16</td><td>4.02</td><td>3.91</td><td>3.88</td><td>⋯</td><td>3.82</td><td>3.92</td><td>4.20</td><td>3.76</td><td>4.24</td><td>4.06</td><td>4.04</td><td>3.92</td><td>3.98</td><td>SNORA22       </td></tr>
	<tr><td>A_33_P3893191</td><td>6.90</td><td>6.57</td><td>5.33</td><td>6.30</td><td>6.70</td><td>6.50</td><td>6.53</td><td>6.13</td><td>6.33</td><td>⋯</td><td>6.06</td><td>6.19</td><td>6.80</td><td>6.04</td><td>6.59</td><td>6.31</td><td>6.62</td><td>6.64</td><td>6.52</td><td>DBIL5P2       </td></tr>
	<tr><td>A_33_P3897734</td><td>3.80</td><td>3.87</td><td>3.81</td><td>3.81</td><td>3.75</td><td>3.59</td><td>3.97</td><td>3.57</td><td>3.72</td><td>⋯</td><td>3.72</td><td>3.69</td><td>3.75</td><td>3.75</td><td>3.59</td><td>3.66</td><td>3.72</td><td>3.93</td><td>3.90</td><td>              </td></tr>
	<tr><td>A_33_P3901921</td><td>3.71</td><td>3.92</td><td>4.15</td><td>3.90</td><td>3.87</td><td>3.97</td><td>4.03</td><td>3.82</td><td>4.07</td><td>⋯</td><td>4.10</td><td>3.98</td><td>3.73</td><td>3.96</td><td>4.01</td><td>4.05</td><td>4.09</td><td>3.94</td><td>4.15</td><td>C12orf48      </td></tr>
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
<ol class=list-inline><li>'GSM2696609'</li><li>'GSM2696610'</li><li>'GSM2696611'</li><li>'GSM2696612'</li><li>'GSM2696613'</li><li>'GSM2696614'</li><li>'GSM2696615'</li><li>'GSM2696616'</li><li>'GSM2696617'</li><li>'GSM2696618'</li><li>'GSM2696619'</li><li>'GSM2696620'</li><li>'GSM2696621'</li><li>'GSM2696622'</li><li>'GSM2696623'</li><li>'GSM2696624'</li><li>'GSM2696625'</li><li>'GSM2696626'</li><li>'GSM2696627'</li><li>'GSM2696628'</li><li>'GSM2696629'</li><li>'GSM2696630'</li><li>'GSM2696631'</li><li>'GSM2696632'</li><li>'GSM2696633'</li><li>'GSM2696634'</li><li>'GSM2696635'</li><li>'GSM2696636'</li><li>'GSM2696637'</li><li>'GSM2696638'</li><li>'GSM2696639'</li><li>'GSM2696640'</li><li>'GSM2696641'</li><li>'GSM2696642'</li><li>'GSM2696643'</li><li>'GSM2696644'</li><li>'GSM2696645'</li><li>'GSM2696646'</li><li>'GSM2696647'</li><li>'GSM2696648'</li><li>'GSM2696649'</li><li>'GSM2696650'</li><li>'GSM2696651'</li><li>'GSM2696652'</li><li>'GSM2696653'</li><li>'GSM2696654'</li><li>'GSM2696655'</li><li>'GSM2696656'</li><li>'GSM2696657'</li><li>'GSM2696658'</li><li>'GSM2696659'</li><li>'GSM2696660'</li><li>'GSM2696661'</li><li>'GSM2696662'</li><li>'GSM2696663'</li><li>'GSM2696664'</li><li>'GSM2696665'</li><li>'GSM2696666'</li><li>'GSM2696667'</li><li>'GSM2696668'</li><li>'GSM2696669'</li><li>'GSM2696670'</li><li>'GSM2696671'</li><li>'GSM2696672'</li><li>'GSM2696673'</li><li>'GSM2696674'</li><li>'GSM2696675'</li><li>'GSM2696676'</li><li>'GSM2696677'</li><li>'GSM2696678'</li><li>'GSM2696679'</li><li>'GSM2696680'</li><li>'GSM2696681'</li><li>'GSM2696682'</li><li>'GSM2696683'</li><li>'GSM2696684'</li><li>'GSM2696685'</li><li>'GSM2696686'</li><li>'GSM2696687'</li><li>'GSM2696688'</li><li>'GSM2696689'</li><li>'GSM2696690'</li><li>'GSM2696691'</li><li>'GSM2696692'</li><li>'GSM2696693'</li><li>'GSM2696694'</li><li>'GSM2696695'</li><li>'GSM2696696'</li><li>'GSM2696697'</li><li>'GSM2696698'</li><li>'GSM2696699'</li><li>'GSM2696700'</li><li>'GSM2696701'</li><li>'GSM2696702'</li><li>'GSM2696703'</li><li>'GSM2696704'</li><li>'GSM2696705'</li><li>'GSM2696706'</li><li>'GSM2696707'</li><li>'GSM2696708'</li><li>'GSM2696709'</li><li>'GSM2696710'</li><li>'GSM2696711'</li><li>'GSM2696712'</li></ol>




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


'GSM2696609=mean(GSM2696609),GSM2696610=mean(GSM2696610),GSM2696611=mean(GSM2696611),GSM2696612=mean(GSM2696612),GSM2696613=mean(GSM2696613),GSM2696614=mean(GSM2696614),GSM2696615=mean(GSM2696615),GSM2696616=mean(GSM2696616),GSM2696617=mean(GSM2696617),GSM2696618=mean(GSM2696618),GSM2696619=mean(GSM2696619),GSM2696620=mean(GSM2696620),GSM2696621=mean(GSM2696621),GSM2696622=mean(GSM2696622),GSM2696623=mean(GSM2696623),GSM2696624=mean(GSM2696624),GSM2696625=mean(GSM2696625),GSM2696626=mean(GSM2696626),GSM2696627=mean(GSM2696627),GSM2696628=mean(GSM2696628),GSM2696629=mean(GSM2696629),GSM2696630=mean(GSM2696630),GSM2696631=mean(GSM2696631),GSM2696632=mean(GSM2696632),GSM2696633=mean(GSM2696633),GSM2696634=mean(GSM2696634),GSM2696635=mean(GSM2696635),GSM2696636=mean(GSM2696636),GSM2696637=mean(GSM2696637),GSM2696638=mean(GSM2696638),GSM2696639=mean(GSM2696639),GSM2696640=mean(GSM2696640),GSM2696641=mean(GSM2696641),GSM2696642=mean(GSM2696642),GSM2696643=mean(GSM2696643),GSM2696644=mean(GSM2696644),GSM2696645=mean(GSM2696645),GSM2696646=mean(GSM2696646),GSM2696647=mean(GSM2696647),GSM2696648=mean(GSM2696648),GSM2696649=mean(GSM2696649),GSM2696650=mean(GSM2696650),GSM2696651=mean(GSM2696651),GSM2696652=mean(GSM2696652),GSM2696653=mean(GSM2696653),GSM2696654=mean(GSM2696654),GSM2696655=mean(GSM2696655),GSM2696656=mean(GSM2696656),GSM2696657=mean(GSM2696657),GSM2696658=mean(GSM2696658),GSM2696659=mean(GSM2696659),GSM2696660=mean(GSM2696660),GSM2696661=mean(GSM2696661),GSM2696662=mean(GSM2696662),GSM2696663=mean(GSM2696663),GSM2696664=mean(GSM2696664),GSM2696665=mean(GSM2696665),GSM2696666=mean(GSM2696666),GSM2696667=mean(GSM2696667),GSM2696668=mean(GSM2696668),GSM2696669=mean(GSM2696669),GSM2696670=mean(GSM2696670),GSM2696671=mean(GSM2696671),GSM2696672=mean(GSM2696672),GSM2696673=mean(GSM2696673),GSM2696674=mean(GSM2696674),GSM2696675=mean(GSM2696675),GSM2696676=mean(GSM2696676),GSM2696677=mean(GSM2696677),GSM2696678=mean(GSM2696678),GSM2696679=mean(GSM2696679),GSM2696680=mean(GSM2696680),GSM2696681=mean(GSM2696681),GSM2696682=mean(GSM2696682),GSM2696683=mean(GSM2696683),GSM2696684=mean(GSM2696684),GSM2696685=mean(GSM2696685),GSM2696686=mean(GSM2696686),GSM2696687=mean(GSM2696687),GSM2696688=mean(GSM2696688),GSM2696689=mean(GSM2696689),GSM2696690=mean(GSM2696690),GSM2696691=mean(GSM2696691),GSM2696692=mean(GSM2696692),GSM2696693=mean(GSM2696693),GSM2696694=mean(GSM2696694),GSM2696695=mean(GSM2696695),GSM2696696=mean(GSM2696696),GSM2696697=mean(GSM2696697),GSM2696698=mean(GSM2696698),GSM2696699=mean(GSM2696699),GSM2696700=mean(GSM2696700),GSM2696701=mean(GSM2696701),GSM2696702=mean(GSM2696702),GSM2696703=mean(GSM2696703),GSM2696704=mean(GSM2696704),GSM2696705=mean(GSM2696705),GSM2696706=mean(GSM2696706),GSM2696707=mean(GSM2696707),GSM2696708=mean(GSM2696708),GSM2696709=mean(GSM2696709),GSM2696710=mean(GSM2696710),GSM2696711=mean(GSM2696711),GSM2696712=mean(GSM2696712),'



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
<caption>A data.frame: 6 × 105</caption>
<thead>
	<tr><th></th><th scope=col>GSM2696609</th><th scope=col>GSM2696610</th><th scope=col>GSM2696611</th><th scope=col>GSM2696612</th><th scope=col>GSM2696613</th><th scope=col>GSM2696614</th><th scope=col>GSM2696615</th><th scope=col>GSM2696616</th><th scope=col>GSM2696617</th><th scope=col>GSM2696618</th><th scope=col>⋯</th><th scope=col>GSM2696704</th><th scope=col>GSM2696705</th><th scope=col>GSM2696706</th><th scope=col>GSM2696707</th><th scope=col>GSM2696708</th><th scope=col>GSM2696709</th><th scope=col>GSM2696710</th><th scope=col>GSM2696711</th><th scope=col>GSM2696712</th><th scope=col>symbol</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>4.45</td><td>4.47</td><td>5.16</td><td>4.96</td><td>4.76</td><td>4.93</td><td>5.03</td><td>5.54</td><td>4.50</td><td>4.98</td><td>⋯</td><td>4.81</td><td>4.82</td><td>4.70</td><td>4.55</td><td>6.17</td><td>4.57</td><td>4.65</td><td>4.62</td><td>4.69</td><td>LOC100130938</td></tr>
	<tr><th scope=row>2</th><td>4.99</td><td>5.12</td><td>4.56</td><td>5.02</td><td>4.88</td><td>5.13</td><td>5.03</td><td>4.92</td><td>4.47</td><td>5.24</td><td>⋯</td><td>4.75</td><td>4.79</td><td>5.14</td><td>4.86</td><td>4.88</td><td>5.30</td><td>5.00</td><td>4.68</td><td>4.86</td><td>            </td></tr>
	<tr><th scope=row>3</th><td>4.00</td><td>4.06</td><td>3.92</td><td>3.93</td><td>3.96</td><td>3.95</td><td>3.88</td><td>3.87</td><td>3.83</td><td>3.87</td><td>⋯</td><td>4.05</td><td>3.82</td><td>4.03</td><td>3.90</td><td>4.08</td><td>3.70</td><td>3.89</td><td>3.92</td><td>3.77</td><td>LOC729177   </td></tr>
	<tr><th scope=row>4</th><td>4.13</td><td>4.41</td><td>4.15</td><td>3.83</td><td>4.11</td><td>3.99</td><td>3.91</td><td>4.09</td><td>3.91</td><td>3.74</td><td>⋯</td><td>4.27</td><td>3.83</td><td>3.94</td><td>3.91</td><td>3.76</td><td>3.96</td><td>4.11</td><td>3.99</td><td>3.83</td><td>Q73P46      </td></tr>
	<tr><th scope=row>5</th><td>6.30</td><td>5.21</td><td>5.10</td><td>5.35</td><td>6.08</td><td>5.82</td><td>6.50</td><td>4.89</td><td>5.22</td><td>5.26</td><td>⋯</td><td>5.57</td><td>5.39</td><td>5.63</td><td>5.21</td><td>5.36</td><td>5.38</td><td>5.44</td><td>6.40</td><td>5.57</td><td>LOC145474   </td></tr>
	<tr><th scope=row>6</th><td>3.91</td><td>3.61</td><td>4.24</td><td>3.76</td><td>3.84</td><td>3.89</td><td>3.77</td><td>3.79</td><td>3.82</td><td>3.70</td><td>⋯</td><td>3.87</td><td>3.94</td><td>3.68</td><td>3.73</td><td>3.72</td><td>3.92</td><td>3.82</td><td>3.82</td><td>3.83</td><td>P01115      </td></tr>
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
              summarise(GSM2696609=mean(GSM2696609),GSM2696610=mean(GSM2696610),GSM2696611=mean(GSM2696611),GSM2696612=mean(GSM2696612),GSM2696613=mean(GSM2696613),GSM2696614=mean(GSM2696614),GSM2696615=mean(GSM2696615),GSM2696616=mean(GSM2696616),GSM2696617=mean(GSM2696617),GSM2696618=mean(GSM2696618),GSM2696619=mean(GSM2696619),GSM2696620=mean(GSM2696620),GSM2696621=mean(GSM2696621),GSM2696622=mean(GSM2696622),GSM2696623=mean(GSM2696623),GSM2696624=mean(GSM2696624),GSM2696625=mean(GSM2696625),GSM2696626=mean(GSM2696626),GSM2696627=mean(GSM2696627),GSM2696628=mean(GSM2696628),GSM2696629=mean(GSM2696629),GSM2696630=mean(GSM2696630),GSM2696631=mean(GSM2696631),GSM2696632=mean(GSM2696632),GSM2696633=mean(GSM2696633),GSM2696634=mean(GSM2696634),GSM2696635=mean(GSM2696635),GSM2696636=mean(GSM2696636),GSM2696637=mean(GSM2696637),GSM2696638=mean(GSM2696638),GSM2696639=mean(GSM2696639),GSM2696640=mean(GSM2696640),GSM2696641=mean(GSM2696641),GSM2696642=mean(GSM2696642),GSM2696643=mean(GSM2696643),GSM2696644=mean(GSM2696644),GSM2696645=mean(GSM2696645),GSM2696646=mean(GSM2696646),GSM2696647=mean(GSM2696647),GSM2696648=mean(GSM2696648),GSM2696649=mean(GSM2696649),GSM2696650=mean(GSM2696650),GSM2696651=mean(GSM2696651),GSM2696652=mean(GSM2696652),GSM2696653=mean(GSM2696653),GSM2696654=mean(GSM2696654),GSM2696655=mean(GSM2696655),GSM2696656=mean(GSM2696656),GSM2696657=mean(GSM2696657),GSM2696658=mean(GSM2696658),GSM2696659=mean(GSM2696659),GSM2696660=mean(GSM2696660),GSM2696661=mean(GSM2696661),GSM2696662=mean(GSM2696662),GSM2696663=mean(GSM2696663),GSM2696664=mean(GSM2696664),GSM2696665=mean(GSM2696665),GSM2696666=mean(GSM2696666),GSM2696667=mean(GSM2696667),GSM2696668=mean(GSM2696668),GSM2696669=mean(GSM2696669),GSM2696670=mean(GSM2696670),GSM2696671=mean(GSM2696671),GSM2696672=mean(GSM2696672),GSM2696673=mean(GSM2696673),GSM2696674=mean(GSM2696674),GSM2696675=mean(GSM2696675),GSM2696676=mean(GSM2696676),GSM2696677=mean(GSM2696677),GSM2696678=mean(GSM2696678),GSM2696679=mean(GSM2696679),GSM2696680=mean(GSM2696680),GSM2696681=mean(GSM2696681),GSM2696682=mean(GSM2696682),GSM2696683=mean(GSM2696683),GSM2696684=mean(GSM2696684),GSM2696685=mean(GSM2696685),GSM2696686=mean(GSM2696686),GSM2696687=mean(GSM2696687),GSM2696688=mean(GSM2696688),GSM2696689=mean(GSM2696689),GSM2696690=mean(GSM2696690),GSM2696691=mean(GSM2696691),GSM2696692=mean(GSM2696692),GSM2696693=mean(GSM2696693),GSM2696694=mean(GSM2696694),GSM2696695=mean(GSM2696695),GSM2696696=mean(GSM2696696),GSM2696697=mean(GSM2696697),GSM2696698=mean(GSM2696698),GSM2696699=mean(GSM2696699),GSM2696700=mean(GSM2696700),GSM2696701=mean(GSM2696701),GSM2696702=mean(GSM2696702),GSM2696703=mean(GSM2696703),GSM2696704=mean(GSM2696704),GSM2696705=mean(GSM2696705),GSM2696706=mean(GSM2696706),GSM2696707=mean(GSM2696707),GSM2696708=mean(GSM2696708),GSM2696709=mean(GSM2696709),GSM2696710=mean(GSM2696710),GSM2696711=mean(GSM2696711),GSM2696712=mean(GSM2696712))  %>% 
         column_to_rownames(var ='symbol')
```


```R
gset$GSE100927_series_matrix.txt.gz@phenoData@data  %>%  head()
```


<table class="dataframe">
<caption>A data.frame: 6 × 32</caption>
<thead>
	<tr><th></th><th scope=col>title</th><th scope=col>geo_accession</th><th scope=col>status</th><th scope=col>submission_date</th><th scope=col>last_update_date</th><th scope=col>type</th><th scope=col>channel_count</th><th scope=col>source_name_ch1</th><th scope=col>organism_ch1</th><th scope=col>characteristics_ch1</th><th scope=col>⋯</th><th scope=col>contact_email</th><th scope=col>contact_institute</th><th scope=col>contact_address</th><th scope=col>contact_city</th><th scope=col>contact_zip/postal_code</th><th scope=col>contact_country</th><th scope=col>supplementary_file</th><th scope=col>data_row_count</th><th scope=col>gender:ch1</th><th scope=col>tissue:ch1</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM2696609</th><td>H1_87: Femoral artery_3                    </td><td>GSM2696609</td><td>Public on Mar 07 2018</td><td>Jul 07 2017</td><td>Mar 07 2018</td><td>RNA</td><td>1</td><td>Atherosclerotic femoral artery</td><td>Homo sapiens</td><td>tissue: Atherosclerotic femoral artery</td><td>⋯</td><td>marja.steenman@inserm.fr</td><td>L'institut du thorax</td><td>8 quai Moncousu</td><td>Nantes</td><td>44007</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2696nnn/GSM2696609/suppl/GSM2696609_US82400123_253949440798_S01_GE1_107_Sep09_1_1.txt.gz</td><td>47264</td><td>male</td><td>Atherosclerotic femoral artery</td></tr>
	<tr><th scope=row>GSM2696610</th><td>H10_109: Femoral artery_7                  </td><td>GSM2696610</td><td>Public on Mar 07 2018</td><td>Jul 07 2017</td><td>Mar 07 2018</td><td>RNA</td><td>1</td><td>Atherosclerotic femoral artery</td><td>Homo sapiens</td><td>tissue: Atherosclerotic femoral artery</td><td>⋯</td><td>marja.steenman@inserm.fr</td><td>L'institut du thorax</td><td>8 quai Moncousu</td><td>Nantes</td><td>44007</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2696nnn/GSM2696610/suppl/GSM2696610_US82400123_253949440799_S01_GE1_107_Sep09_1_2.txt.gz</td><td>47264</td><td>male</td><td>Atherosclerotic femoral artery</td></tr>
	<tr><th scope=row>GSM2696611</th><td>H102_8: Femoral artery_Control 26          </td><td>GSM2696611</td><td>Public on Mar 07 2018</td><td>Jul 07 2017</td><td>Mar 07 2018</td><td>RNA</td><td>1</td><td>Control femoral artery        </td><td>Homo sapiens</td><td>tissue: Control femoral artery        </td><td>⋯</td><td>marja.steenman@inserm.fr</td><td>L'institut du thorax</td><td>8 quai Moncousu</td><td>Nantes</td><td>44007</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2696nnn/GSM2696611/suppl/GSM2696611_US82400123_253949440824_S01_GE1_107_Sep09_2_2.txt.gz</td><td>47264</td><td>male</td><td>Control femoral artery        </td></tr>
	<tr><th scope=row>GSM2696612</th><td>H103_90: Carotid artery_Control 14         </td><td>GSM2696612</td><td>Public on Mar 07 2018</td><td>Jul 07 2017</td><td>Mar 07 2018</td><td>RNA</td><td>1</td><td>Control carotid artery        </td><td>Homo sapiens</td><td>tissue: Control carotid artery        </td><td>⋯</td><td>marja.steenman@inserm.fr</td><td>L'institut du thorax</td><td>8 quai Moncousu</td><td>Nantes</td><td>44007</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2696nnn/GSM2696612/suppl/GSM2696612_US82400123_253949440824_S01_GE1_107_Sep09_2_3.txt.gz</td><td>47264</td><td>male</td><td>Control carotid artery        </td></tr>
	<tr><th scope=row>GSM2696613</th><td>H104_154: Infra-popliteal artery_Control 01</td><td>GSM2696613</td><td>Public on Mar 07 2018</td><td>Jul 07 2017</td><td>Mar 07 2018</td><td>RNA</td><td>1</td><td>Control infra-popliteal artery</td><td>Homo sapiens</td><td>tissue: Control infra-popliteal artery</td><td>⋯</td><td>marja.steenman@inserm.fr</td><td>L'institut du thorax</td><td>8 quai Moncousu</td><td>Nantes</td><td>44007</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2696nnn/GSM2696613/suppl/GSM2696613_US82400123_253949440824_S01_GE1_107_Sep09_2_4.txt.gz</td><td>47264</td><td>male</td><td>Control infra-popliteal artery</td></tr>
	<tr><th scope=row>GSM2696614</th><td>H105_25: Carotid artery_91                 </td><td>GSM2696614</td><td>Public on Mar 07 2018</td><td>Jul 07 2017</td><td>Mar 07 2018</td><td>RNA</td><td>1</td><td>Atherosclerotic carotid artery</td><td>Homo sapiens</td><td>tissue: Atherosclerotic carotid artery</td><td>⋯</td><td>marja.steenman@inserm.fr</td><td>L'institut du thorax</td><td>8 quai Moncousu</td><td>Nantes</td><td>44007</td><td>France</td><td>ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2696nnn/GSM2696614/suppl/GSM2696614_US82400123_253949440823_S01_GE1_107_Sep09_1_1.txt.gz</td><td>47264</td><td>male</td><td>Atherosclerotic carotid artery</td></tr>
</tbody>
</table>




```R
sample_type  <- factor(gset$GSE100927_series_matrix.txt.gz@phenoData@data[,c('tissue:ch1')]) 
```


```R
levels(sample_type)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Atherosclerotic carotid artery'</li><li>'Atherosclerotic femoral artery'</li><li>'Atherosclerotic infra-popliteal artery'</li><li>'Control carotid artery'</li><li>'Control femoral artery'</li><li>'Control infra-popliteal artery'</li></ol>




```R
group_levels  <- tibble(sample_type) %>% 
 mutate(sites = fct_recode(sample_type,'Carotid' = 'Atherosclerotic carotid artery',
                           'Carotid' = 'Control carotid artery',
                           'Popliteal' = 'Atherosclerotic infra-popliteal artery',
                           'Popliteal' = 'Control infra-popliteal artery',
                           'Femoral' = 'Atherosclerotic femoral artery',
                           'Femoral' = 'Control femoral artery')) %>% 
  mutate(sites = factor(sites, levels = c('Carotid', 'Popliteal', 'Femoral'))) %>% 
   mutate(group = fct_recode(sample_type,
                           'ATH' = 'Atherosclerotic carotid artery',
                           'CON' = 'Control carotid artery',
                           'ATH' = 'Atherosclerotic infra-popliteal artery',
                           'CON' = 'Control infra-popliteal artery',
                           'ATH' = 'Atherosclerotic femoral artery',
                           'CON' = 'Control femoral artery')) %>%
  mutate(group = factor(group, levels = c('CON', 'ATH')))
```


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
gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
     bind_cols(group_levels)
```


<table class="dataframe">
<caption>A data.frame: 104 × 6</caption>
<thead>
	<tr><th scope=col>sample</th><th scope=col>Yellow_Module</th><th scope=col>Brown_Module</th><th scope=col>sample_type</th><th scope=col>sites</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>GSM2696609</td><td>2.922172</td><td>3.241956</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696610</td><td>2.947696</td><td>3.418891</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696611</td><td>2.795125</td><td>3.352487</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696612</td><td>2.641102</td><td>2.970704</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696613</td><td>2.915399</td><td>3.118116</td><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>GSM2696614</td><td>3.101599</td><td>3.590428</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696615</td><td>2.982895</td><td>3.563155</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696616</td><td>2.738367</td><td>3.011996</td><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>GSM2696617</td><td>2.767654</td><td>3.147533</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696618</td><td>2.854903</td><td>3.394276</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696619</td><td>2.958601</td><td>3.481471</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696620</td><td>3.066650</td><td>3.593381</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696621</td><td>3.042409</td><td>3.472217</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696622</td><td>3.096925</td><td>3.638442</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696623</td><td>3.101370</td><td>3.586684</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696624</td><td>2.933332</td><td>3.197251</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696625</td><td>2.993798</td><td>3.308700</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696626</td><td>2.964801</td><td>3.369674</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696627</td><td>2.893417</td><td>3.520325</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696628</td><td>3.026598</td><td>3.586363</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696629</td><td>2.881737</td><td>3.238043</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696630</td><td>3.065198</td><td>3.579133</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696631</td><td>2.979947</td><td>3.438506</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696632</td><td>2.980167</td><td>3.464823</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696633</td><td>3.071148</td><td>3.632596</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696634</td><td>3.013482</td><td>3.417732</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696635</td><td>2.939110</td><td>3.180019</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696636</td><td>2.902881</td><td>3.104513</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696637</td><td>3.099116</td><td>3.575938</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696638</td><td>2.999599</td><td>3.386269</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>GSM2696683</td><td>2.916666</td><td>3.207415</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696684</td><td>2.887342</td><td>3.193931</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696685</td><td>2.995799</td><td>3.376843</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696686</td><td>3.109810</td><td>3.578816</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696687</td><td>2.923261</td><td>3.295921</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696688</td><td>2.868956</td><td>3.314209</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696689</td><td>3.057154</td><td>3.502628</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696690</td><td>2.981822</td><td>3.367701</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696691</td><td>2.744221</td><td>3.189460</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696692</td><td>3.000323</td><td>3.474908</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696693</td><td>2.941347</td><td>3.282064</td><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>GSM2696694</td><td>2.796333</td><td>3.008676</td><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>GSM2696695</td><td>2.752885</td><td>3.039661</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696696</td><td>3.077242</td><td>3.552464</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696697</td><td>2.944958</td><td>3.518554</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696698</td><td>2.967306</td><td>3.369459</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696699</td><td>2.812704</td><td>3.042287</td><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>GSM2696700</td><td>2.984998</td><td>3.550577</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696701</td><td>3.062510</td><td>3.517986</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696702</td><td>2.964885</td><td>3.521599</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696703</td><td>2.971498</td><td>3.302083</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696704</td><td>2.879398</td><td>3.283639</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696705</td><td>2.950479</td><td>3.435543</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696706</td><td>2.939623</td><td>3.306607</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>GSM2696707</td><td>2.995583</td><td>3.436781</td><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>GSM2696708</td><td>3.012688</td><td>3.443319</td><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>GSM2696709</td><td>2.762053</td><td>3.027811</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696710</td><td>2.841682</td><td>3.250898</td><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>GSM2696711</td><td>3.017484</td><td>3.493204</td><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>GSM2696712</td><td>2.951231</td><td>3.382648</td><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
</tbody>
</table>




```R
group_levels$group %>% table()
group_levels
```


    .
    CON ATH 
     35  69 



<table class="dataframe">
<caption>A tibble: 104 × 3</caption>
<thead>
	<tr><th scope=col>sample_type</th><th scope=col>sites</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Control femoral artery                </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery        </td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic carotid artery        </td><td>Carotid  </td><td>ATH</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Control infra-popliteal artery        </td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>Control carotid artery                </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Atherosclerotic infra-popliteal artery</td><td>Popliteal</td><td>ATH</td></tr>
</tbody>
</table>




```R
head(group_levels)
group_levels$sites %>% unique()
```


<table class="dataframe">
<caption>A tibble: 6 × 3</caption>
<thead>
	<tr><th scope=col>sample_type</th><th scope=col>sites</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Atherosclerotic femoral artery</td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Atherosclerotic femoral artery</td><td>Femoral  </td><td>ATH</td></tr>
	<tr><td>Control femoral artery        </td><td>Femoral  </td><td>CON</td></tr>
	<tr><td>Control carotid artery        </td><td>Carotid  </td><td>CON</td></tr>
	<tr><td>Control infra-popliteal artery</td><td>Popliteal</td><td>CON</td></tr>
	<tr><td>Atherosclerotic carotid artery</td><td>Carotid  </td><td>ATH</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>Femoral</li><li>Carotid</li><li>Popliteal</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Carotid'</li><li>'Popliteal'</li><li>'Femoral'</li></ol>
</details>



```R
group_levels$group  %>% table()
```


    .
    CON ATH 
     35  69 



```R
p_yellow  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
       ggboxplot(x = 'group', y = 'Yellow_Module', ylab = 'GSVA Score',xlab = '',title = 'All sites',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_yellow
```


    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_30_0.png)
    



```R
p_yellow_carotid  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
      filter(sites == 'Carotid') %>% 
       ggboxplot(x = 'group', y = 'Yellow_Module', ylab = 'GSVA Score',xlab = '',title = 'Carotid artery',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_yellow_carotid
p_yellow_Popliteal  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
      filter(sites == 'Popliteal') %>% 
       ggboxplot(x = 'group', y = 'Yellow_Module', ylab = 'GSVA Score',xlab = '',title = 'Popliteal artery',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_yellow_Popliteal
p_yellow_Femoral  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
      filter(sites == 'Femoral') %>% 
       ggboxplot(x = 'group', y = 'Yellow_Module', ylab = 'GSVA Score',xlab = '',title = 'Femoral artery',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_yellow_Femoral
```


    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_31_0.png)
    



    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_31_1.png)
    



    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_31_2.png)
    



```R
library(patchwork)
```


```R
p_brown  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
       ggboxplot(x = 'group', y = 'Brown_Module', ylab = 'GSVA Score',xlab = '',title = 'All sites',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_brown
```


    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_33_0.png)
    



```R
p_brown_carotid  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
      filter(sites == 'Carotid') %>% 
       ggboxplot(x = 'group', y = 'Brown_Module', ylab = 'GSVA Score',xlab = '',title = 'Carotid artery',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_brown_carotid
p_brown_Popliteal  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
      filter(sites == 'Popliteal') %>% 
       ggboxplot(x = 'group', y = 'Brown_Module', ylab = 'GSVA Score',xlab = '',title = 'Popliteal artery',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_brown_Popliteal
p_brown_Femoral  <- gsva_matrix  %>% 
 t() %>% 
     as.data.frame()  %>% 
     rownames_to_column(var = 'sample') %>% 
      bind_cols(group_levels) %>% 
      filter(sites == 'Femoral') %>% 
       ggboxplot(x = 'group', y = 'Brown_Module', ylab = 'GSVA Score',xlab = '',title = 'Femoral artery',add = 'point',fill = 'group',legend = 'right',palette = 'lancet')  + 
         stat_compare_means(comparisons = list(c('CON','ATH')),method = 't.test', paired = F,label = 'p.signif') + theme(plot.title = element_text(hjust = 0.5))  +
         theme(plot.title=element_text(hjust=0.5, size = 8),
            axis.text.x = element_text(hjust=1,size=6),
             axis.text.y = element_text(size=6),
              axis.title.y = element_text(size=8),
                legend.title =  element_text(size = 6),
                 legend.text = element_text(size = 6))
p_brown_Femoral
```


    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_34_0.png)
    



    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_34_1.png)
    



    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_34_2.png)
    



```R
p <- p_yellow + p_yellow_carotid + p_yellow_Popliteal + p_yellow_Femoral + 
p_brown + p_brown_carotid + p_brown_Popliteal + p_brown_Femoral + plot_layout(ncol = 4,guides='collect')  +  
plot_annotation(title = 'GSE100927',theme = theme(plot.title = element_text(size = 8, hjust = 0.5)))
p
ggsave(p, file = './GSE100927_Two_Module_1.pdf', height = 10, width = 16, units = 'cm')
```


    
![png](Step5.4_Two_Group_Compare_GSE100927_files/Step5.4_Two_Group_Compare_GSE100927_35_0.png)
    

