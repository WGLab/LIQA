# Instruction about how to use LIQA

The inputs of LIQA are aligned long-read RNA-seq data in BAM format and a reference isoform annotation file (Ensembl/Refseq). User needs to specify `-task` to perform in each step:
```
LIQA.py -task <task>:

        refgene:   preprocess reference file

        group:   group alternative splicing exon

        count:   count informative reads from indexed BAM file

        gene:    estimate mean gene expression for each single cell condition

        das:    detect differential alternative splicing (DAS) for each exon group between conditions

        sum:    summarize DAS test results

```

## Step 1: Transform isoforms to compatible matrix based on reference annotation file
LIQA requires a reference annotation file `example.refFile` in following format:
```
749     NM_001397       chr1    -       21543739        21616982        21546447        21616907        19      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,21616562,21616856,  21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21605825,21616649,21616982,     0       ECE1cmpl     cmpl    0,0,1,2,2,0,0,0,0,2,0,0,0,0,1,1,0,0,0,
93      NM_001113348    chr1    -       21543739        21672034        21546447        21671871        19      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,21616562,21671868,  21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21605825,21616649,21672034,     0       ECE1cmpl     cmpl    0,0,1,2,2,0,0,0,0,2,0,0,0,0,1,1,0,0,0,
749     NM_001113349    chr1    -       21543739        21616766        21546447        21616691        18      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,21616562,   21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21605825,21616766,      0       ECE1    cmpl    cmpl0,0,1,2,2,0,0,0,0,2,0,0,0,0,1,1,0,0,
749     NM_001113347    chr1    -       21543739        21606183        21546447        21605927        17      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,    21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21606183,       0       ECE1    cmpl    cmpl    0,0,1,2,2,0,0
,0,0,2,0,0,0,0,1,1,0,
```
Reference file in this format can be downloaded at [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?command=start) by selecting "all fields from selected table" in output format.

We preprocess `example.refFile` by using `python LIQA.py -task refgene`. An example is given below.
```
python LIQA.py -task refgene -ref example.refFile -out example.refgene
```

## Step 2: Extract informative read count for each exon group from alignment file
LIQA requires a headerless `metafile` in this step to tell LIQA that how and where to find the aligment BAM files to extract cell-specific informative read count. BAM files have to be indexed. Here is an example of `metafile`
```
AACACGTCACATAACC-1      A       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
GGACAAGTCTCCCTGA-1      A       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
CACAGGCAGATCCCGC-1      B       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
ATCTGCCGTCATCGGC-1      B       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
GGAAAGCGTTGCTCCT-1      C       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
CGAGCACGTGTTCTTT-1      C       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
CCTATTACAATGGATA-1      D       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
AAGGAGCAGCGTCAAG-1      D       ~/1dot1/outs/possorted_genome_bam.bam  UB      CB
```
where <strong>1st</strong> column contains cell <strong>barcode/cell name</strong>, <strong>2nd</strong> column represents <strong>condition group</strong>, <strong>3rd</strong> column represents the <strong>location of BAM file</strong>. <strong>4th</strong> and <strong>5th</strong> columns represent the <strong>tag names of UMI barcode and cell barcode</strong> in BAM file. For example
```
NS500497:57:H27CKBGX2:3:12506:1885:16376        272     1       3014861 1       98M     *       0       0       TGGCGTTCCCCTGTACTGGGGCTTATAAAGTTTGCAAGTCCAATGGGCCTCTCTTTGCAGTGATGGCCGACTAGGCCATCTTTTGATACATATGCAGC      //A/A/A/EEE<A66EA/EAE//</EAEE//E/AEEAE/EEEAEE//AEEEE/AAAEEEEAEEEEE6EEEEEEEEEEEEAEE/EE6EEEEAEEAAAAA   NH:i:3  HI:i:3  AS:i:94 nM:i:1  NM:i:1  CR:Z:GAGGTGAAGTGACATA   CY:Z:AAAAAEEEEEEEEEEE   CB:Z:GAGGTGAAGTGACATA-1 UR:Z:CCATACATGA UY:Z:EEEEEEEEEE UB:Z:CCATACATGA      BC:Z:GGTTTACT   QT:Z:AAAAAEEE   RG:Z:CellRangerCount-1dot1_combined:MissingLibrary:1:H27CKBGX2:3
```
where UMI barcode and cell barcode tag names are indicated by "UB" and "CB" in this BAM file.

The specific usage details are given below.
```
        USAGE:
        python LIQA.py -task count [count options] -meta <metafile> -refgene <refgene_file> -gpinfo <gpinfo_file>

        [count options]    type 'python LIQA.py -task count' to check two important count options.
        
          -umi  <yes/no> collect UMI count or not
        
          -onebam  <yes/no> whether all aligned reads are merged in one BAM files

        OUTPUT:
        
          count_*.sh    script files will be generated under directory `./tmp/count_script`.
```
where '-umi' and '-onebam' are two important options:
* `-umi yes -onebam yes`: UMI and cell barcode tag names have to be specified in the <strong>4th</strong> and <strong>5th</strong> columns of `metafile`.
* `-umi yes -onebam no`: only UMI barcode tag name is needed. It has to be specified in the <strong>4th</strong> column of `metafile`.
* `-umi no -onebam yes`: only cell barcode tag name is needed. It has to be specified in the <strong>4th</strong> column of `metafile`.
* `-umi no -onebam no`: no tag name is needed.

Outputs of `python LIQA.py -task count` are `count_*.sh` script files located at `./tmp/count_script`. User needs to run all of them to obtain informative read count for each single cell.

## Step 3: Quantify gene-level expression accounting for technical noises
In this step, user needs to give  `metafile` to LIQA and specify the number of cores to use for each pairwise comparison between conditions:
```
python LIQA.py -task gene -ncore 20 -meta metafile
```
Outputs of `python LIQA.py -task gene` are `gene_*.sh` script files located at `./tmp/gene_script`. User needs to run all of them to obtain accurate gene expression estimations for each cell condition group.

## Step 4: Detect differential alternative splicing (DAS) across cell conditions accounting for technical noises
In this step, user needs to give  `metafile` and `example.gpinfo` to LIQA and specify the number of cores to use for each pairwise comparison between conditions:
```
python LIQA.py -task das -ncore 20 -meta metafile -gpinfo example.gpinfo
```
Outputs of `python LIQA.py -task das` are `das_*.sh` script files located at `./tmp/das_script`. User needs to run all of them to obtain differential alternative splicing even at exon group level across cell conditions.

## Step 5: Summarize DAS test results
```
python LIQA.py -task sum -gpinfo example.gpinfo
```
