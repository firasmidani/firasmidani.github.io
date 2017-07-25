---
layout: post
title: "Python functions for microbiome exploratory data analysis."
date: 2017-07-24
use_math: true
---

```python
def joinTables(list_of_tables):
    '''
    INPUT
    	list_of_tables : list of pandas.DataFrames
    OUTPUT
    	new_table : pandas.DataFrame of outer concatenation of tables across columns
    EXAMPLE
    	T1_D789 = JoinTables([M[1][df][1] for df in [7,8,9]])
    '''

    new_table = pd.concat(list_of_tables,axis=0,join='outer').fillna(0)
    
    return new_table 
```

```php
def joinTables(list_of_tables):
    '''
    INPUT
    	list_of_tables : list of pandas.DataFrames
    OUTPUT
    	new_table : pandas.DataFrame of outer concatenation of tables across columns
    EXAMPLE
    	T1_D789 = JoinTables([M[1][df][1] for df in [7,8,9]])
    '''

    new_table = pd.concat(list_of_tables,axis=0,join='outer').fillna(0)
    
    return new_table
```
