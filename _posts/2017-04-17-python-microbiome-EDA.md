---
layout: post
title: "Python functions for microbiome exploratory data analysis."
date: 2017-04-17
use_math: true
comments: true
---




```python
import pandas as pd
import numpy as np
import imp

import pprint
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import spearmanr,kruskal,mannwhitneyu

%matplotlib inline
sns.set_style('whitegrid')
```


```python
imp.load_source('biom_library','/Users/firasmidani/Downloads/microbiome_routines/biom_library.py');
from biom_library import *
```

## List of functions

* `otuTaxaDict()`
* `subsetTableByMetadata()`
* `subsetTableyBySampleIDs()`
* `summarizeTable()`
* `tss_norm()`
* `phyloSummaryOtuTableSingle()`
* `read_pcoa_file()`
* `pcoa_figure()`


## Import BIOM file
Break down the BIOM file into two items:
1. otu_table:  OTUs by samples
2. otu_taxa_map: mapping of OTU ID to taxonomy


```python
biom = pd.read_csv('./tables/otus_table.filtered.from_biom.txt',sep='\t',header=0,index_col=0,skiprows=1)
otu_taxa_map = pd.DataFrame(biom.iloc[:,-1])
otu_table = biom.iloc[:,:-1]
```


```python
print otu_taxa_map.shape
otu_taxa_map.head()
```

    (2226, 1)





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>taxonomy</th>
    </tr>
    <tr>
      <th>#OTU ID</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1027904</th>
      <td>k__Bacteria; p__Firmicutes; c__Bacilli; o__Bac...</td>
    </tr>
    <tr>
      <th>1050608</th>
      <td>k__Bacteria; p__Actinobacteria; c__Actinobacte...</td>
    </tr>
    <tr>
      <th>127870</th>
      <td>k__Bacteria; p__Proteobacteria; c__Alphaproteo...</td>
    </tr>
    <tr>
      <th>816470</th>
      <td>k__Bacteria; p__Firmicutes; c__Bacilli; o__Bac...</td>
    </tr>
    <tr>
      <th>177792</th>
      <td>k__Bacteria; p__Firmicutes; c__Clostridia; o__...</td>
    </tr>
  </tbody>
</table>
</div>




```python
print otu_table.shape
otu_table.head()
```

    (2226, 138)





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>113</th>
      <th>30</th>
      <th>15</th>
      <th>145</th>
      <th>116</th>
      <th>72</th>
      <th>37</th>
      <th>49</th>
      <th>64</th>
      <th>73</th>
      <th>...</th>
      <th>130</th>
      <th>87</th>
      <th>88</th>
      <th>115</th>
      <th>4</th>
      <th>31</th>
      <th>117</th>
      <th>123</th>
      <th>133</th>
      <th>99</th>
    </tr>
    <tr>
      <th>#OTU ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1027904</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1050608</th>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>127870</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>816470</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>9.0</td>
      <td>4.0</td>
      <td>3.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>2.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>177792</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 138 columns</p>
</div>



Most of the functions in this library are designed to accept an OTU sample with samples as rows and OTUs as columns. So, let's transpose the OTU table. 


```python
otu_table = otu_table.T
```

## Delineate taxonomic levels for each OTU 


```python
otu_taxa_dict = otuTaxaDict(otu_taxa_map)
otu_taxa_dict.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>p</th>
      <th>c</th>
      <th>o</th>
      <th>f</th>
      <th>g</th>
      <th>s</th>
    </tr>
    <tr>
      <th>#OTU ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1027904</th>
      <td>Firmicutes</td>
      <td>Bacilli</td>
      <td>Bacillales</td>
      <td>Staphylococcaceae</td>
      <td>Macrococcus</td>
      <td></td>
    </tr>
    <tr>
      <th>1050608</th>
      <td>Actinobacteria</td>
      <td>Actinobacteria</td>
      <td>Actinomycetales</td>
      <td>Corynebacteriaceae</td>
      <td>Corynebacterium</td>
      <td></td>
    </tr>
    <tr>
      <th>127870</th>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Sphingomonadales</td>
      <td>Sphingomonadaceae</td>
      <td>Novosphingobium</td>
      <td></td>
    </tr>
    <tr>
      <th>816470</th>
      <td>Firmicutes</td>
      <td>Bacilli</td>
      <td>Bacillales</td>
      <td>Bacillaceae</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>177792</th>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
</div>



## Subsettting an OTU table by sample IDs


```python
subsetTableBySampleIDs(otu_table,[90,93,145])
```




    'ERROR : Some samples are missing from the table.'




```python
subsetTableBySampleIDs(otu_table,['90','93','145'])
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>#OTU ID</th>
      <th>816470</th>
      <th>360717</th>
      <th>3571769</th>
      <th>192773</th>
      <th>1111582</th>
      <th>309696</th>
      <th>524318</th>
      <th>196893</th>
      <th>290468</th>
      <th>4369988</th>
      <th>...</th>
      <th>518743</th>
      <th>196664</th>
      <th>531675</th>
      <th>198555</th>
      <th>583117</th>
      <th>315846</th>
      <th>237444</th>
      <th>560336</th>
      <th>580008</th>
      <th>351231</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>90</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>9.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>53.0</td>
    </tr>
    <tr>
      <th>93</th>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>4.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>3.0</td>
    </tr>
    <tr>
      <th>145</th>
      <td>9.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>7.0</td>
      <td>7.0</td>
      <td>1.0</td>
      <td>5.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>48.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>33.0</td>
      <td>69.0</td>
    </tr>
  </tbody>
</table>
<p>3 rows × 201 columns</p>
</div>



## Subsettting an OTU table by metada

First, we need to import metadata 


```python
mapping_df = pd.read_csv('./2017_03_15_mapping_processed.txt',sep='\t',header=0,index_col=0)
print mapping_df.shape
mapping_df.head()
```

    (138, 6)





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>BarcodeSequence</th>
      <th>Plate</th>
      <th>Row</th>
      <th>Column</th>
      <th>Time/Comment</th>
      <th>NumReads</th>
    </tr>
    <tr>
      <th>#SampleID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>113</th>
      <td>ACAATGTCACAG</td>
      <td>2.0</td>
      <td>E</td>
      <td>7.0</td>
      <td>0</td>
      <td>1199.0</td>
    </tr>
    <tr>
      <th>30</th>
      <td>GAAAGGTGAGAA</td>
      <td>4.0</td>
      <td>B</td>
      <td>8.0</td>
      <td>0</td>
      <td>4394.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>TCAATGACCGCA</td>
      <td>2.0</td>
      <td>C</td>
      <td>7.0</td>
      <td>0</td>
      <td>1333.0</td>
    </tr>
    <tr>
      <th>145</th>
      <td>GAATCCTCACCG</td>
      <td>4.0</td>
      <td>C</td>
      <td>7.0</td>
      <td>0</td>
      <td>1697.0</td>
    </tr>
    <tr>
      <th>116</th>
      <td>GTCCCTATTATC</td>
      <td>1.0</td>
      <td>H</td>
      <td>3.0</td>
      <td>0</td>
      <td>2957.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
parameters_dict = {'Plate':[1,2,3],'Row':['C'],'Time/Comment':['24']}
new_table = subsetTableByMetadata(otu_table, mapping_df,parameters_dict)

pprint.pprint(parameters_dict,width=50)
new_table
```

    {'Plate': [1, 2, 3],
     'Row': ['C'],
     'Time/Comment': ['24']}





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>#OTU ID</th>
      <th>310380</th>
      <th>360717</th>
      <th>221784</th>
      <th>4242681</th>
      <th>309696</th>
      <th>290468</th>
      <th>4369988</th>
      <th>182116</th>
      <th>4366525</th>
      <th>365628</th>
      <th>...</th>
      <th>696563</th>
      <th>523542</th>
      <th>193769</th>
      <th>346639</th>
      <th>301270</th>
      <th>196660</th>
      <th>198555</th>
      <th>237444</th>
      <th>580008</th>
      <th>339512</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>46</th>
      <td>3.0</td>
      <td>1.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>128.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>2996.0</td>
      <td>6.0</td>
      <td>20.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>125</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>109.0</td>
      <td>3.0</td>
      <td>13.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>230.0</td>
      <td>0.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>18.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>129</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>78.0</td>
      <td>0.0</td>
      <td>9.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>...</td>
      <td>2.0</td>
      <td>33.0</td>
      <td>0.0</td>
      <td>2.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>12.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>102</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>103.0</td>
      <td>0.0</td>
      <td>11.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>14.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>39</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>66.0</td>
      <td>0.0</td>
      <td>9.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>7.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>10.0</td>
      <td>1.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 130 columns</p>
</div>




```python
summarizeTable(new_table,otu_taxa_dict).head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>mean</th>
      <th>median</th>
      <th>count</th>
      <th>min</th>
      <th>max</th>
      <th>p</th>
      <th>c</th>
      <th>o</th>
      <th>f</th>
      <th>g</th>
      <th>s</th>
    </tr>
    <tr>
      <th>#OTU ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>310380</th>
      <td>0.6</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Lachnospiraceae</td>
      <td>Dorea</td>
      <td></td>
    </tr>
    <tr>
      <th>360717</th>
      <td>0.2</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Lachnospiraceae</td>
      <td>Blautia</td>
      <td></td>
    </tr>
    <tr>
      <th>221784</th>
      <td>0.4</td>
      <td>2.0</td>
      <td>1.0</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Lachnospiraceae</td>
      <td>Blautia</td>
      <td></td>
    </tr>
    <tr>
      <th>4242681</th>
      <td>0.2</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Lachnospiraceae</td>
      <td>Dorea</td>
      <td></td>
    </tr>
    <tr>
      <th>309696</th>
      <td>71.2</td>
      <td>90.5</td>
      <td>4.0</td>
      <td>66.0</td>
      <td>109.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Clostridiaceae</td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
</div>




```python
summarizeTable(new_table,otu_taxa_dict).sort_values(['mean'],ascending=False).head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>mean</th>
      <th>median</th>
      <th>count</th>
      <th>min</th>
      <th>max</th>
      <th>p</th>
      <th>c</th>
      <th>o</th>
      <th>f</th>
      <th>g</th>
      <th>s</th>
    </tr>
    <tr>
      <th>#OTU ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>349024</th>
      <td>3645.8</td>
      <td>4748.0</td>
      <td>5.0</td>
      <td>1.0</td>
      <td>5398.0</td>
      <td>Firmicutes</td>
      <td>Bacilli</td>
      <td>Lactobacillales</td>
      <td>Streptococcaceae</td>
      <td>Streptococcus</td>
      <td></td>
    </tr>
    <tr>
      <th>523542</th>
      <td>653.2</td>
      <td>131.5</td>
      <td>4.0</td>
      <td>7.0</td>
      <td>2996.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Lachnospiraceae</td>
      <td>Dorea</td>
      <td></td>
    </tr>
    <tr>
      <th>309696</th>
      <td>71.2</td>
      <td>90.5</td>
      <td>4.0</td>
      <td>66.0</td>
      <td>109.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Clostridiaceae</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>523782</th>
      <td>54.8</td>
      <td>3.5</td>
      <td>4.0</td>
      <td>2.0</td>
      <td>265.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Lachnospiraceae</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>363321</th>
      <td>43.8</td>
      <td>18.0</td>
      <td>3.0</td>
      <td>5.0</td>
      <td>196.0</td>
      <td>Firmicutes</td>
      <td>Clostridia</td>
      <td>Clostridiales</td>
      <td>Lachnospiraceae</td>
      <td>Dorea</td>
      <td></td>
    </tr>
  </tbody>
</table>
</div>




```python
new_table.head().sum(1)
```




    46     4516.0
    125    5899.0
    129    4983.0
    102    5032.0
    39     3408.0
    dtype: float64



## Convert absolute abundances to relativea abudnances with total sum scaling


```python
norm_otu_table = tss_norm(otu_table)
print norm_otu_table.sum().head()
norm_otu_table.head()
```

    113    1.0
    30     1.0
    15     1.0
    145    1.0
    116    1.0
    dtype: float64





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>113</th>
      <th>30</th>
      <th>15</th>
      <th>145</th>
      <th>116</th>
      <th>72</th>
      <th>37</th>
      <th>49</th>
      <th>64</th>
      <th>73</th>
      <th>...</th>
      <th>130</th>
      <th>87</th>
      <th>88</th>
      <th>115</th>
      <th>4</th>
      <th>31</th>
      <th>117</th>
      <th>123</th>
      <th>133</th>
      <th>99</th>
    </tr>
    <tr>
      <th>#OTU ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1027904</th>
      <td>0.000834</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1050608</th>
      <td>0.000000</td>
      <td>0.000228</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>127870</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00075</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>816470</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.005303</td>
      <td>0.001353</td>
      <td>0.000779</td>
      <td>0.000621</td>
      <td>0.000508</td>
      <td>0.000671</td>
      <td>0.000978</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>177792</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 138 columns</p>
</div>




```python
norm_otu_table_0 = tss_norm(otu_table,axis=0)
print norm_otu_table_0.sum().head()
norm_otu_table_0.head()
```

    #OTU ID
    1027904    1.0
    1050608    1.0
    127870     1.0
    816470     1.0
    177792     1.0
    dtype: float64





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>#OTU ID</th>
      <th>1027904</th>
      <th>1050608</th>
      <th>127870</th>
      <th>816470</th>
      <th>177792</th>
      <th>337735</th>
      <th>193591</th>
      <th>182674</th>
      <th>310380</th>
      <th>366068</th>
      <th>...</th>
      <th>194925</th>
      <th>194924</th>
      <th>351231</th>
      <th>197214</th>
      <th>525698</th>
      <th>577377</th>
      <th>4417325</th>
      <th>503354</th>
      <th>259955</th>
      <th>361398</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>113</th>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.043478</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.096815</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>30</th>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.003822</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.033333</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.014013</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>145</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.333333</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.087898</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>116</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.148148</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.021656</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.005556</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 2226 columns</p>
</div>



Let us generate another one and join these two tables.

## Standard statistical and taxonomic summary of a table


```python
imp.load_source('biom_library','/Users/firasmidani/Downloads/microbiome_routines/biom_library.py');
from biom_library import *
```


```python
phyloSummaryOtuTableSingle(tss_norm(new_table),otu_taxa_dict,['p','c'])
```




    OrderedDict([('p',                Num. OTUs Rel. Abundance
                  Firmicutes           118       0.986714
                  Bacteroidetes          9      0.0104074
                  Actinobacteria         3     0.00287865),
                 ('c',                 Num. OTUs Rel. Abundance
                  Clostridia             95       0.986492
                  Bacteroidia             9      0.0104074
                  Actinobacteria          2     0.00243578
                  Coriobacteriia          1     0.00044287
                  Bacilli                22    0.000221435
                  Erysipelotrichi         1              0)])




```python
new_table_summary = phyloSummaryOtuTableSingle(tss_norm(new_table),otu_taxa_dict)
new_table_summary['g']
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Num. OTUs</th>
      <th>Rel. Abundance</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Dorea</th>
      <td>29</td>
      <td>0.795837</td>
    </tr>
    <tr>
      <th></th>
      <td>42</td>
      <td>0.170283</td>
    </tr>
    <tr>
      <th>[Ruminococcus]</th>
      <td>6</td>
      <td>0.0110717</td>
    </tr>
    <tr>
      <th>Bacteroides</th>
      <td>8</td>
      <td>0.0104074</td>
    </tr>
    <tr>
      <th>Blautia</th>
      <td>11</td>
      <td>0.00531444</td>
    </tr>
    <tr>
      <th>Coprococcus</th>
      <td>4</td>
      <td>0.00376439</td>
    </tr>
    <tr>
      <th>Bifidobacterium</th>
      <td>2</td>
      <td>0.00243578</td>
    </tr>
    <tr>
      <th>Collinsella</th>
      <td>1</td>
      <td>0.00044287</td>
    </tr>
    <tr>
      <th>Streptococcus</th>
      <td>21</td>
      <td>0.000221435</td>
    </tr>
    <tr>
      <th>Roseburia</th>
      <td>2</td>
      <td>0.000221435</td>
    </tr>
    <tr>
      <th>Phascolarctobacterium</th>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Parabacteroides</th>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Lactobacillus</th>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Veillonella</th>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>



## Principal coordinate analysis figure using QIIME-formatted principal components


```python
pcoa_dict = read_pcoa_file('./pc/pcoa_weighted_unifrac_otus_table.filtered.txt')
print pcoa_dict.keys()
```

    ['eigen_vectors', 'eigen_values', 'var_explained']



```python
fig,ax = plt.subplots(figsize=[4,4])

classes = {'24_1':[{'Time/Comment':['24'],'Column':[1]},(0.0,0.0,1.0,0.60),'^','24_1'],
           '24_2':[{'Time/Comment':['24'],'Column':[2]},(0.5,0.0,0.5,0.60),'o','24_2'],
           '24_3':[{'Time/Comment':['24'],'Column':[3]},(1.0,0.0,0.0,0.60),'*','24_3']}

pcoa_figure(pcoa_dict,ax,['0','1'],'Example PCoA',mapping_df,classes)
```


![png](output_31_0.png)



```python
pcoa_dict = {};
for distance_metric in ['weighted_unifrac','unweighted_unifrac','bray_curtis','binary_jaccard']: 
    pcoa_dict[distance_metric] = read_pcoa_file('./pc/pcoa_%s_otus_table.filtered.txt' % distance_metric)
```


```python
#plt.rc('font',family='sans-serif',serif='Arial')

import matplotlib.gridspec as gridspec

fig   = plt.figure(figsize=(13,12))
gs    = gridspec.GridSpec(2,2,height_ratios=[1,1],hspace=0.4,wspace=0.4)

ax_={};
ax_[0] =  plt.subplot(gs[0,0])
ax_[1] =  plt.subplot(gs[0,1])
ax_[2] =  plt.subplot(gs[1,0])
ax_[3] =  plt.subplot(gs[1,1])

classes = {'24_1':[{'Time/Comment':['24'],'Column':[1]},(0.0,0.0,1.0,0.60),'^','24_1'],
           '24_2':[{'Time/Comment':['24'],'Column':[2]},(0.5,0.0,0.5,0.60),'o','24_2'],
           '24_3':[{'Time/Comment':['24'],'Column':[3]},(1.0,0.0,0.0,0.60),'*','24_3']}


pcoa_figure(pcoa_dict['weighted_unifrac'],    ax_[3], ['0','1'], 'Weighted Unifrac Distance', mapping_df, classes);
pcoa_figure(pcoa_dict['unweighted_unifrac'],  ax_[2], ['0','1'], 'Unweighted Unifrac Distance', mapping_df, classes);
pcoa_figure(pcoa_dict['bray_curtis'],         ax_[1], ['0','1'], 'Bray-Curtis Dissimilarity', mapping_df, classes);
pcoa_figure(pcoa_dict['binary_jaccard'],      ax_[0], ['0','1'], 'Jaccard Distance', mapping_df, classes);

legend = ax_[3].legend(bbox_to_anchor=(1.1,0.55),bbox_transform=plt.gcf().transFigure,fontsize=20,frameon=True)
legend.get_frame().set_linewidth(3)

plt.suptitle('Grid of PCoA PLots',fontsize=30,fontweight='bold',y=1.02)

plt.gcf().subplots_adjust(left=0.075,right=0.925,top=0.925,bottom=0.075)

# plt.savefig('./figs/20170411_beta_diversity_panel.pdf',format='pdf')
```


![png](output_33_0.png)



```python

```
