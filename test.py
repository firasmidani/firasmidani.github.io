#ENA Upload Protocol

@(david.documents)[protocol, computing, data management]

| **Author**  Firas Said Midani
| **E-mail**  firas.midani@duke.edu
| **Date Created** 2017.01.13
| **Date Modified** 2017.01.13

id: 20160524-215821

[TOC]

##Required files

You need four files:

1. Index read file
2. Forward read file
3. Reverse read file
4. Mapping file


The sequencing core should provide you with the original index, forward, and reverse read files. You should have created a mapping file with proper format. The mapping file, at but at the bare minimum, it needs **\#SampleID**, **BarcodeSequence**, **LinkerPrimerSequence**, and **Description** columns.  It can also include additional columns and should look like the following

|\#SampleID|BarcodeSequence|LinkerPrimerSequence|Row|Column|Plate|Description|
|:---|:---|:---|:---|:---|:---|:---|
|SMIC.030.05.Day.02     |TCCCTTGTCTCC| CAAGCAGAAGACGGCATACGAGAT|A |1 |1 |806rcbc0|
|SMIC.040.01.Day.02     |TGCATACACTGG| CAAGCAGAAGACGGCATACGAGAT|B |1 |1 |806rcbc12|
|SMIC.014.01.Day.02     |GCGATATATCGC| CAAGCAGAAGACGGCATACGAGAT|C |1 |1 |806rcbc24|

##Desired output

For more details, see ENA <a href="https://www.ebi.ac.uk/ena/submit/read-file-formats"> Supported read file formats</a> page.

Briefly, the European Nucleotide Archives need sample-specific FASTQ files. Each FASTQ file must only include reads corresponding to a single sample. If paired end sequencing was performed, the forward and reverse reads must be in separate files.  Here is the ENA excerpt about the expected FASTQ format.

>**Fastq format**
>
>We recommend that read data is either submitted in BAM or CRAM format. However, single and paired reads are accepted as Fastq files that meet the following the requirements:
>* Quality scores must be in Phred scale. Both ASCII and space delimitered decimal encoding of quality scores are supported. We will automatically detect the Phred quality offset of either 33 or 64.
>* No technical reads (adapters, linkers, barcodes) are allowed.
>* Single reads must be submitted using a single Fastq file and can be submitted with or without read names.
>*Paired reads must split and submitted using either one or two Fastq files. The read names must have a suffix identifying the first and second read from the pair, for example '/1' and '/2' (regular expression for the reads: "^@([a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9-]+:[0-9-]+) ([12]):[YN]:[0-9]*[02468]:[ACGTN]+$").
>* The first line for each read must start with '@'.
>* The base calls and quality scores must be separated by a line starting with '+'.
>* The Fastq files must be compressed using gzip or bzip2.

By default (checked) for MiSeq Runs, the sequencing core provides us with FASTQ files with

- [x] Quality scores in Phred scale.
- [x] No technical reads (adapters, linkers, barcodes) included. Adapters and linkers were trimmed and barcode included seperately in the index read file.
- [ ] The reads names **DO NOT** have a suffix identifying the first and second read from the pair, for example '/1' and '/2'.
- [x] The first line for each starts with '@'.
- [x] The base calls and quality scores are sperated by a line starting with '+'.
- [ ] The fastq files compressed using gzip, **BUT** as you will see below,  we will have to decompress the files for proper formatting.

Accordingly, we will have to manipulate our raw FASTQ files to get them to the proper format.

## General outline of the process:

- [ ] For each run (in case, your project samples were processed over multiple runs):
     - [ ] Create a mapping file for the samples (in case, your sequencing run contains samples for other projects) to be used by qiime
     - [ ] For each forward and reverse reads separately
          - [ ] De-multiplex reads into one large FASTQ file
          - [ ] Distribute reads into sample-specific FASTQ files.
          - [ ] Format headers of sequences in FASTQ files to desired ENA format.
          - [ ] Compress files with Gzip.
          - [ ] Compute MD5 128-bit hashes for the compressed files with ```md5sum```
- [ ] Move all sample-specific compressed forward and reverse into a storage area that allows FTP transfer (e.g. isilon storage)
- [ ] Move MD5 hashes to local machine in order to incoroporate into ENA mapping files (needed for submission).

----------

###STEP 1 Create mapping file for your samples

Often, sequencing runs contain samples from a hodge-podge of experiments. To process the sequencing file in QIIME (for microbial community profiling), you should have already constructed mapping file. Here, you will simply need to remove rows for samples that you do not need. Only keep rows for samples that will be included in the ENA repository.  There are several ways to do it. You can use Excel or you can your favorite scripting language (```python```, ```R```,  or even ```bash```/```awk```  if you're adventurous).

----------

###STEP 2 De-multiplex reads into one large FASTQ file

The sequencing core provides you with three FASTQ files for the forward reads, reverse reads, and the index reads. The reads are organized in the same order across all three files. In other words, the 100th read in the forward file is paired with the 100th read in the reverse and index files. De-multiplexing takes the forward and reverse reads and modify their headers to include information about the barcode (and thus sample) corresponding to the reads.  

Here is what an index read looks like.
```tex
@M00830:212:000000000-A9BBF:1:1101:15315:1443 1:N:0:0
AACGAGAACTGA
+
A>A>>11AC1F@
```
</br>
Here is what the corresponding forward read looks.
```tex
@M00830:212:000000000-A9BBF:1:1101:15315:1443 1:N:0:0
TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGATTGGTCAGTCTGTCTTAAAAGTTCGGGGCTTAACCCCGTGATGGGATGGAAACTGCCAATCTAGAGTATCGGAGAGGAAAGTGGAATTCCTAGTG
+
>1>1>@B>B1>CA11FEFFEGGHF?00BGFHHHBEGEGGGHFDEGGAE/E?ECGC/>?G0112BGHHHHHHH>GFHH1>2?//>BCHFGG1/?/<B/0?F/?FFGGGFFFHG01?1<<1<1FFHFGGCA?-ACGGBC0=DBGDGCH00C;;
```
```tex
@M00830:212:000000000-A9BBF:1:1101:15315:1443 1:N:0:0
AACGAGAACTGA
+
A>A>>11AC1F@
```
</br>
De-multiplexing takes these two files and combines them into one with more informative header. For example, once the above read is de-multiplexed, this is what you get:

```tex
>SMIC.017.03.Day.07_0 M00830:212:000000000-A9BBF:1:1101:15315:1443 1:N:0:0 orig_bc=AACGA GAACTGA new_bc=AACGAGAACTGA bc_diffs=0
2. TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGATTGGTCAGTCTGTCTTAAAAGTTCGGGGCTTA ACCCCGTGATGGGATGGAAACTGCCGATCTAGAGTATCGGAGAGGAAAGTGGAATTCCTAGTG
3. +
4. >1>1>@B>B1>CA11FEFFEGGHF?00BGFHHHBEGEGGGHFDEGGAE/E?ECGC/>?G0112BGHHHHHHH>GFHH1>2?//>BCHF
   GG1/?/<B/0?F/?FFGGGFFFHG01?1<<1<1FFHFGGCA?-ACGGBC0=DBGDGCH00C;;
```
where
* **SMIC.07.03.Day.07_0** is the sample name as indicated from the mapping file and **_0** indicates that this is the first read in the FASTQ files (de-multiplexing with QIIME uses <a href="https://skillcrush.com/2013/01/17/why-programmers-start-counting-at-zero/">0-based indexing</a>).
* **orig_bc**  indicate the corresponding index sequenced in the run.
* **new_bc** indicates the corresponding barcode in the mapping file that matches the index sequenced. By default, QIIME pipeline allows for a single base error in the index/barcode. So, this new_bc might be different than orig_bc.
* **bc_diff** indicates the number of base differences between orig_bc and new_bc.
* the rest of the in-between are the default output by the Illumina MiSeq machine. Some of it is informative but you do not necessarily need it. If you are curious about what these numbers mean, see <a href="http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm">this document</a> by Illumina. Anyways, we will come back to this in a bit.

#### So, how do you de-multiplex with QIIME?

```bash
#!/bin/sh

#SBATCH --mem=2GB
#SBATCH --time-min=600

module load qiime

seq_path=/home/lad44/davidlab/users/fsm/seq_data/2014.07.07.seq

srun -o $seq_path/split.library.output.no.joining/slurm.log/split.library.forward.out \
     -e $seq_path/split.library.output.no.joining/slurm.log/split.library.forward.err \
     split_libraries_fastq.py \
     -m $seq_path/mapping/2014.07.07.mapping_corrected_cholera.only.txt \
     -b $seq_path/rawdata/decompressed/Undetermined_S0_L001_I1_001.fastq \
     -i $seq_path/rawdata/decompressed/Undetermined_S0_L001_R1_001.fastq \
     -o $seq_path/split.library.output.no.joining/seqs.forward.fna \
     --rev_comp_mapping_barcodes \
     --store_demultiplexed_fastq \
     --retain_unassigned_reads \
     -r 999 -p 0.0001 -n 999 -q 0 \
     --max_barcode_errors 1.5
```

As you have noticed, I am running my commands as jobs on a high-throughput computing (HPC) cluster. In our case, the cluster is HARDAC. I am basically\*
* requesting 8GB resources for 10 hours (I always ask for ten times the time I think I'll need)
* loading the QIIME module (into the environment so that my code can use QIIME functions),
* ```srun``` to request a job where
     *  we specify the path of the output and error log files with ```-o``` and ```-e``` arguments.
     *  we use the command ```split_libraries_fastq.py``` from QIIIME.
     *  ```-m``` corresponds to path of mapping file
     *  ```-b``` corresponds to path of index read file
     *  ```-i``` corresponds to path of forward (or reverse) read file
     *  ```-o``` specifies path of the output of de-multiplex FASTQ file and its name
     *  ```--rev_comp_mapping_barcodes``` forces lookup of the reverse complement barcode reads
     *  ```--store_demultiplexed_fastq``` forces output into FASTQ format otherwise it would be a FASTA file without quality scores.
\* Don't forget to to do this for reverse reads as well.

Time estimate: de-multiplexing takes about 20 minutes for each GB of data. If you process forward and reverse reads in series, obviously, the process would take the double time.

----------

### STEP 3 Distribute reads into sample-specific files

De-multiplexing should produce a single large FASTQ file for the forward reads and another large FASTQ files for the reverse. In these files, the header of each read indicates the sample name. The example above identified the read to correspond to *SMIC.017.03.Day.07*.  The distribution process basically reads each header and copies the read to a new file specific to its original sample. For example, the read in the example would be appended to the file specific to *SMIC.017.03.Day.07*.

#### So, how do you distribute samples with QIIME?
```bash
#!/bin/sh

#SBATCH --mem=2GB
#SBATCH --time-min=600

module load qiime

seq_path=/home/lad44/davidlab/users/fsm/seq_data/2014.07.07.seq

srun -o $seq_path/split.library.output.no.joining/slurm.log/distribute.library.forward.out \
     -e $seq_path/split.library.output.no.joining/slurm.log/distribute.library.forward.err \
     split_sequence_file_on_sample_ids.py \
     -i $seq_path/split.library.output.no.joining/seqs.forward.fna/seqs.fastq \
     -o $seq_path/split.library.output.no.joining/split_by_sample/forward/original --file_type fastq
```
Here, I am\*
* requesting 8GB resources for 10 hours (I always ask for ten times the time I think I'll need)
* loading the QIIME module (into the environment so that my code can use QIIME functions),
* ```srun``` to requet a job where
     *  we specify the path of the output and error log files with ```-o``` and ```-e``` arguments.
     *  we use the command ```split_sequence_file_on_sample_ids.py``` from QIIIME.
     *  ```-i``` corresponds to path of forward (or reverse) read file
     *  ```-o``` specifies path of the output of de-multiplex FASTQ file and its name
     *  ```--file_type fastq``` specifies the output to remain in FASTQ format

\* Don't forget to to do this for reverse reads as well.

Time estimate: de-multiplexing takes about 2.5 minutes for each GB of data. If you process forward and reverse reads in series, obviously, the process would take the double time.

----------

### STEP 4 Formatting headers to proper ENA format

As I discussed in the [Desired output](#Desired output) section, we need to modify the headers of all of the reads in the FASTQ files. In particular, we need to do the following:

* For each read, add suffix identifying whether it is the forward or reverse from a pair, for example '/1' or '/2'
* Remove any unnecessary detail in the header that is included by the defaults of the MiSeq (optional)

In regards to the latter, let me explain briefly the different information included in the header. I am basing the following on <a href="http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm">this document</a> by Illumina.  For reference, let's revisit the example of the forward read I showed earlier.

```tex
@M00830:212:000000000-A9BBF:1:1101:15315:1443 1:N:0:0
TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGATTGGTCAGTCTGTCTTAAAAGTTCGGGGCTTAACCCCGTGATGGGATGGAAACTGCCAATCTAGAGTATCGGAGAGGAAAGTGGAATTCCTAGTG
+
>1>1>@B>B1>CA11FEFFEGGHF?00BGFHHHBEGEGGGHFDEGGAE/E?ECGC/>?G0112BGHHHHHHH>GFHH1>2?//>BCHFGG1/?/<B/0?F/?FFGGGFFFHG01?1<<1<1FFHFGGCA?-ACGGBC0=DBGDGCH00C;;
```

The elements in the header are the following.

|Element|Description|Example from case above|
|:---|:---|
|`@`|Each sequence identifier line starts with @|
|`<instrument>`|Instrument ID|M00830|
|`<run number>`|Run number on instrument|212|
|`flowcell ID>`|ID of flow cell|000000000-A9BBF|
|`<lane>`|Flow cell lane number|1
|`<tile>`|Lame tile number|1101
|`<x_pos>`|X coordinate of cluster|15315|
|`<y_pos>`|Y coordinate of cluster|1443|
|`<read>`|Read number, 1 can be single read or read 2 of paired end|1|
|`<is filtered>`|Y if the read is filtered, N otherwise|N|
|`<control number>`|0 when none of the control bits are on, otherwise it is an even number|0|
|`<index sequence>`|Index sequence|0 (because index reads are sparate in our case)|

Most of this information is likely never to be used by you or anyone else. However, it is important to maintain some of it for provenance. For the example below, I only dismiss the run number and the read number, read filter status, and control number. I also appended to the end of the header a hashtag `#` followed by the index read, then a slash `\` followed by the read number (`#AACGAGAACTGA/1`). See the below example

Before formatting header

```tex
@SMIC.017.03.Day.07_0 M00830:212:000000000-A9BBF:1:1101:15315:1443 1:N:0:0 orig_bc=AACGAGAACTGA new_bc=AACGAGAACTGA bc_diffs=0
TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGATTGGTCAGTCTGTCTTAAAAGTTCGGGGCTTAACCCCGTGATGGGATGGAAACTGCCAATCTAGAGTATCGGAGAGGAAAGTGGAATTCCTAGTG
+
>1>1>@B>B1>CA11FEFFEGGHF?00BGFHHHBEGEGGGHFDEGGAE/E?ECGC/>?G0112BGHHHHHHH>GFHH1>2?//>BCHFGG1/?/<B/0?F/?FFGGGFFFHG01?1<<1<1FFHFGGCA?-ACGGBC0=DBGDGCH00C;;
```

After formatting headers

```tex
@M00830:000000000-A9BBF:1:1101:15315:1443#AACGAGAACTGA/1
TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGATTGGTCAGTCTGTCTTAAAAGTTCGGGGCTTAACCCCGTGATGGGATGGAAACTGCCAATCTAGAGTATCGGAGAGGAAAGTGGAATTCCTAGTG
+
>1>1>@B>B1>CA11FEFFEGGHF?00BGFHHHBEGEGGGHFDEGGAE/E?ECGC/>?G0112BGHHHHHHH>GFHH1>2?//>BCHFGG1/?/<B/0?F/?FFGGGFFFHG01?1<<1<1FFHFGGCA?-ACGGBC0=DBGDGCH00C;;
```

#### So, how do you actually format the headers efficiently?

We won't be using QIIME here. I wrote my own script in `python` (see below). It acts two arguments:
* absolute path to FASTQ file*
* read number (1 or 2) indicating whether file contains forward or reverse reads

\* If you give the script a FASTA file, it will give you back a FASTQ file which would not make sense. FASTA files do not have any quality scores for the bases. So, to properly return a FASTQ file, the script by default uses the base calls as the quality scores which is wrong. Ideally, this script would include error checking to either not accept FASTA files or return FASTA files in the same format. If I have time, I will do this <sup>yea ... right</sup>
<br/>

**ENA_format_fastq.py**
```python
#!/usr/bin/env python

#Author #Firas Said Midani
#E_mail #firas.midani@duke.edu
#Date Created #2016.11.27

from Bio import SeqIO
import sys

input_seq_file = sys.argv[1] ## file name with absolute path
seq_direction = sys.argv[2] ## 1 or 2 for forward or reverse respectively

print 'input_seq_file\t',input_seq_file

output_seq_file = input_seq_file.split('.fastq')[0]
output_seq_file = output_seq_file.replace('original','ENA_formatted') # make sure ENA_Formatted directory exists
output_seq_file = '%s.fastq' % output_seq_file

print '~~~ Formatting ~~~'
print '   input %s' % input_seq_file
print '   output %s' % output_seq_file

input_seq_handle = open(input_seq_file,"rU")
output_seq_handle = open(output_seq_file,"w+")

for seq_record in SeqIO.parse(input_seq_handle,"fastq"):

    old_description = seq_record.description;

    id_h, info_r, info_f, bc_orig, bc_new, bc_diff = seq_record.description.split(' ');

    new_description = '%s:%s#%s/%s' % (info_r.split(':')[0],
                                       ":".join(info_r.split(':')[2:]),
                                       bc_orig.split('=')[-1],
                                       seq_direction)

    seq_record.description = "";
    seq_record.id = new_description;

    SeqIO.write(seq_record,output_seq_handle,"fastq")

print '~~~~~~~~~~~~~~~~~~'
print

input_seq_handle.close();
output_seq_handle.close();
```
The above script only formats one FASTQ file at a time. So, I wrote up a shell script to loop through all of the files of interest. I also submits this script as the job on the cluster.
<br/>

**ENA_format_fastq_library.sh**
```bash
#!/bin/sh

#SBATCH --time-min=120
#SBATCH --mem=2GB

source /home/lad44/davidlab/users/fsm/smic_study/vpython/bin/activate

seq_path=/home/lad44/davidlab/users/fsm/seq_data/2014.07.07.seq
seq_library=$seq_path/split.library.output.no.joining

for i in $seq_library/split_by_sample/forward/original/*;
    do python $seq_library/ENA_format_fastq.py $i 1;
done;

for i in $seq_library/split_by_sample/reverse/original/*;
    do python $seq_library/ENA_format_fastq.py $i 2;
done;
```
<br/>
**ENA_format_fastq_library.slurm**
```bash
#!/bin/sh

#SBATCH --time-min=120
#SBATCH --mem=2GB

seq_path=/home/lad44/davidlab/users/fsm/seq_data/2014.07.07.seq
seq_library=$seq_path/split.library.output.no.joining

srun -o $seq_library/slurm.log/ENA_format_fastq_library.out\
     -e $seq_library/slurm.log/ENA_format_fastq_library.err \
     sh $seq_library/ENA_format_fastq_library.sh
```

<br/>
Time estimate: de-multiplexing takes abou 2 minutes for each GB of data. If you process forward and reverse reads in series, as shown in the example above, the process would take the double time.

----------
### STEP 5 Compressing file with Gzip and computing their MD5 hashes

This is pretty straightforward except that `gzip` by default over-writes the input file.  To keep both the raw and compressed version of the FASTQ files, we need to be extra verbose. Basically, we loop through the FASTQ files of interest in our directory. Define the new file name and its file path. Then, compress each file and pipe its compressed to the new filename. FInally, compute the MD5 hash for each file and save it in the same directory.

**gzip_fastq_library.sh**
```bash
#!/bin/sh

#SBATCH --time-min=120
#SBATCH --mem=2GB

seq_path=/home/lad44/davidlab/users/fsm/seq_data/2014.07.07.seq
seq_library=$seq_path/split.library.output.no.joining

for i in $seq_library/split_by_sample/forward/ENA_formatted/*;
    do
        out_file=$(echo $i | sed 's/ENA_formatted/ENA_formatted_gz/');
        echo $out_file;
        gzip -c $i > $out_file.gz;
        md5sum $out_file.gz > $out_file.gz.md5sum
done;


for i in $seq_library/split_by_sample/reverse/ENA_formatted/*;
    do
        out_file=$(echo $i | sed 's/ENA_formatted/ENA_formatted_gz/');
        echo $out_file;
        gzip -c $i > $out_file.gz;
        md5sum $out_file.gz > $out_file.gz.md5sum
    done
```


**gzip_fastq_library.slurm**
```bash
#!/bin/sh

#SBATCH --time-min=120
#SBATCH --mem=2GB

seq_path=/home/lad44/davidlab/users/fsm/seq_data/2014.07.07.seq
seq_library=$seq_path/split.library.output.no.joining

srun -o $seq_library/slurm.log/gzip_fastq_library.out\
     -e $seq_library/slurm.log/gzip_fastq_library.err \
     sh $seq_library/gzip_fastq_library.sh                                         
```
<br/>
Time estimate: de-multiplexing takes about 1 minute for each GB of data. If you process forward and reverse reads in series, as shown in the example above, the process would take the double time.

#### Computing MD5 hashes of gzippped FASTQ files into a single file
In some cases, it is convenient to have all of the MD5 hashes in the same file. For example, when registering samples in the European Nucleotide Archive, you might have to include the MD5 hashes in its own columns along with other meta-data regarding each sample.

To do so, I submit the ```md5sum_library.sh``` to the cluster with the only argument as the absolute path leading to the FASTQ library.
```bash
sbatch md5sum_library.sh data/davidlab/users/fsm/seq_data/2014.07.07.seq/split.library.output.no.joining/split_by_sample/
```
**md5sum_library.sh**
```bash
#!/bin/sh

#SBATCH --time-min=120
#SBATCH --mem=2GB

md5sum $1/forward/ENA_formatted_gz/*fastq.gz > $1/forward/ENA_formatted_gz/forward.md5sum
md5sum $1/reverse/ENA_formatted_gz/*fastq.gz > $1/reverse/ENA_formatted_gz/reverse.md5sum
```

The output are two files (one for forward reads and one for reverse reads) that look like this. I shortened the full path name only for easier visualization.
|MD5 hash|file|
|:---|:---|
|ac4040592106e5359f5d0c07af7a599f|/home/lad44/.../forward/ENA_formatted_gz/SMIC.001.01.Day.02.fastq.gz|
|9f9aa5018335741dcd159f423ec48cf7|/home/lad44/.../forward/ENA_formatted_gz/SMIC.001.02.Day.02.fastq.gz|
|d8a64e2920bb461f7113c9c69b184375|/home/lad44/.../forward/ENA_formatted_gz/SMIC.001.02.Day.10.fastq.gz|

Time estimate: This is very fast. It might only take a minute for a MiSeq run sequencing data.

----------
The protocols should discuss briefly the metadata needed by ENA for all samples. There are two types of metadata: (1) metadata about the biological samples, and (2) metadata about the processing/sequencing of the samples. Below is an example of the metadata I include with the sequences from the cholera project uploaded to ENA.

**(1) Metadata about the biological samples (cholera project example)**
|key|value|example|
|:---|:---|:---|
|sample_alias|---|SMIC.001.01.Day.02|
|sample_title|---|1|
|investigation type|---|mimarks-survey|
|project name|---|cholera susceptibility|
|sequencing method|---|Illumina|
|collection date|---|2013|
|geographic location (country and/or sea)|---|Bangladesh|
|geographic location (latitude)|---|23.7|
|geographic location (longitude)|---|90.3667|
|human gut environmental package     |---|human-gut|
|environment (biome)|---|ENVO:urban biome|
|environment (feature)|---|ENVO:human-associated habitat|
|environment (material)|---|ENVO:feces|

**(2) Metadata about prcessing/sequencing of biological samples (cholera project example)**
|key|value|example|
|:---|:---|:---|
|project_accession|this is generated by ENA|PRJEB17860|
|sample_alias|---|SMIC.001.01.Day.02|ena-EXPERIMENT-DUKE UNIV-27-11-2016-23:51:21:382-1|
|run_alias|this is generated by ENA|ena-RUN-DUKE UNIV-27-11-2016-23:51:21:383-1|
|library_name|---|unspecified|
|library_source|---|METAGENOMIC|
|library_selection|---|PCR|
|library_strategy|---|AMPLICON|
|instrument_model|---|Illumina MiSeq|
|file_type|---|fastq|
|library_layout|---|PAIRED|
|insert_size|---|151|
|forward_file_name|---|SMIC.001.01.Day.02_1.fastq|
|reverse_file_md5|---|714fc74ed2f5a2b0bf3fa7944541bfd6|
|reverse_file_name|---|SMIC.001.01.Day.02_2.fastq|
|reverse_file_md5|---|4272a174c217dd193d2f26b18080200c|
