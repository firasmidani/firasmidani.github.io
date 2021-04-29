---
layout: post
title: "New Software: AMiGA for the Analysis of Microbial Growth Assays"
date: 2021-04-23
use_math: true
comments: False
---

Check out my new software AMiGA (see <a href="https://github.com/firasmidani/amiga">code</a> and <a href="https://firasmidani.github.io/amiga">documentation</a>) for analyzing microbial growth curves. At the <a href="https://twitter.com/BrittonTML">@BrittonTML</a>, we had a pretty good system for analyzing microbial growth especially on Biolog PM plates, but we focused on a single metric which missed some intereting growth dynamics. 

I started with published tools and writing my own code to fit growth curves using classical models like Gompertz but quickly realized that clinical isoaltes including *Clostridioides difficile* are quite finicky and have varying shapes and growth dynamics. 

<br/>
{:refdef: style="text-align: center;"}
![amiga_1](/assets/img/amiga/amiga_1.png){:width="500px"}
{: refdef}
<br/>

To get around this issue, I adopted a fairly nascent approach of modelling growth curves with Gaussian Proces regression. It worked really well, and the rest of the team and collaborators became interested in this approach. 

<br/>
{:refdef: style="text-align: center;"}
![amiga_2](/assets/img/amiga/amiga_2.png){:width="400px"}
{: refdef}
<br/>

So I took my scripts and turned them into a user-friendly program **A**nalysis of **Mi**crobial **G**rowth **A**sssays (**AMiGA**). I spent lots of time on documentation and examples, see [https://firasmidani.github.io/amiga](https://firasmidani.github.io/amiga).

<br/>
{:refdef: style="text-align: center;"}
![amiga_3](/assets/img/amiga/amiga_3.png){:width="500px"}
{: refdef}
<br/>

AMiGA can fit individual curves or pool replicates (from different plates<a></a>) and give you confidence intervals too. It can estimate death rates; and detect & characterize diauxic shifts. I'm not aware of any existing software that automatically quantifies eitehr phenotypes, so feedback welcome. 

<br/>
<p style="text-align: center;">
  <img src="/assets/img/amiga/amiga_4a.png" width="200" />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img src="/assets/img/amiga/amiga_4b.png" width="200" /> 
</p>
<br/>

Users can interact with AMiGA via the command-lien (terminal). It requires minimal and often no pre-processing of your raw data.

<br/>
{:refdef: style="text-align: center;"}
![amiga_5](/assets/img/amiga/amiga_5.png){:width="500px"}
{: refdef}
<br/>

You can also throw in some meta-data to test hypotheses of differential growth. Feel free to report bugs or ask questions on email or GitHub repo. 

<br/>
<p style="text-align: center;">
  <img src="/assets/img/amiga/amiga_6a.png" width="400" />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <br/>
  <br/>
  <img src="/assets/img/amiga/amiga_6b.png" width="400" /> 
</p>
<br/>


# Paper Details

FS Midani, J Collins, and RA Britton.<br>
**AMiGA: software for the Analysis of Microbial Growth Assays**<br>
*bioRxiv*. 2020. <a href="https://www.biorxiv.org/content/10.1101/2020.11.04.369140v1">https://www.biorxiv.org/content/10.1101/2020.11.04.369140v1</a>

**Abstract**  <br>
The analysis of microbial growth is one of the central methods in the field of microbiology. Microbial growth dynamics can be characterized by growth parameters including carrying capacity, exponential growth rate, and growth lag. However, growth assays with clinical isolates, fastidious organisms, or microbes under stress often produce atypical growth shapes that do not follow the classical microbial growth pattern. Here, we introduce the Analysis of Microbial Growth Assays (AMiGA<a></a>) software which streamlines the analysis of growth curves without any assumptions about their shapes. AMiGA can pool replicates of growth curves and infer summary statistics for biologically meaningful growth parameters. In addition, AMiGA can quantify death phases and characterize diauxic shifts. It can also statistically test for differential growth under distinct experimental conditions. Altogether, AMiGA streamlines the organization, analysis, and visualization of microbial growth assays.

**Importance** <br> 
Our current understanding of microbial physiology relies on the simple method of measuring microbial populations’ size over time and under different conditions. Many advances have increased the throughput of those assays and enabled the study of non-lab adapted microbes under diverse conditions that widely affect their growth dynamics. Our software provides an all-in-one tool for estimating the growth parameters of microbial cultures and testing for differential growth in a high-throughput and user-friendly fashion without any underlying assumptions about how microbes respond to their growth conditions.

<br><br><br>
![footer_banner](/assets/img/mosaic_footer.png)
