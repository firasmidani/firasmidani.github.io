---
layout: post
title: "Python functions for microbiome exploratory data analysis."
date: 2017-04-17
use_math: true
---

Over the past couple of years, I have been repeatedly munging through microbiome data and also repeatedly re-writing the same code to analyze similar datasets. On more occasions, I have been asked by peers and colleagues to share bits of the code that I have written. Eventually, I forced myself to create a central repository for code routinely used by me to analyze 16S amplicon sequencing data sets. In the David lab, we run our sequencing results through the [QIIME](http://qiime.org) pipeleine and accordingly most of the data analysis begins with a [BIOM](http://biom-format.org) file.

Here, I am sharing some of the primary functions that make exploration and analysis of microbiome data simpler.

|Function|Description|
|:---|:---|
|`read_biom`|reads a `BIOM` file and returns an OTU table and a seperate OTU to taxonomy mapping|
|`otu_taxa_dict`|breaks down the one-to-one OTU-to-taxonmy mapping to all taxonomic levels (phylum through species<a></a>)|
|`tss_norm`|converts the absolute abundances in an OTU table to relative abundances|
|`subsetTableBySampleIDs`|reutrns a susbet of the input OTU table based on the sample IDs of your choice|
|`subsetTableByMetadata`|returns a subset of the input OTU table based on metadata parameters of your choice|
|`summarizeTable`|returns a statistical summary of the distribution of OTU abundance across samples in the table|
|`phyloSummaryOtuTableSingle`|returns a taxonomy/phylogenic summary statistics of the OTUs comprising sample(s) in a table|
|`read_pcoa`|reads the `principal_coordinates.py` output by QIIME and extracts coordinates (eigenvectors/values<a></a>) and vairances explained|
|`pcoa_figure`|plots a PCoA figure with user-defined subset of data and figure aesthetics|


Here is a tutorial of these functions in a [jupyter notebook](/assets/ipynbs/2017_04_10_seq_analysis_post.html) and you can find the underyling functions in [firasmidani/microbiome-routines repository](https://github.com/firasmidani/microbiome_routines).
